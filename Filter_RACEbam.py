import pysam
import re
import argparse
from Bio.Seq import Seq
import gzip
from Levenshtein import distance

parser=argparse.ArgumentParser()

parser.add_argument("-p","--primers_file",help="Input file with primers", required=True)
parser.add_argument("-f","--bam_file",help="Input bam file", required=True)
parser.add_argument("-b","--skip_bp",help="Skip this number of bp from read start", default=0, type=int)
parser.add_argument("-e","--ed_primer",help="Edit distance with primer", default=5, type=int)
parser.add_argument("-l","--tail_seq",help="5'-tail to remove from each primer", default="")
parser.add_argument("-d","--max_dist",help="The max tolerated distance from primer 5'-end to the first aligned coordinate in read", default=20, type=int)
parser.add_argument("-o","--output_prefix",help="Prefix for output files", required=True)

args=parser.parse_args()

#find primer sequence and its reference coordinates
primers_file=open(args.primers_file,'rt')
input_header=primers_file.readline().strip("\n").split("\t")
primer_dict={}
primer_res_dict={}
for line in primers_file:
    line_dict={}
    line_vals=line.strip("\n").split("\t")
    for i in range(0,len(input_header)):
        line_dict[input_header[i]]=line_vals[i]
    primer_res_dict[line_dict["name"]]=[]
    primer_seq=(line_dict["primer_seq"])[len(args.tail_seq):]
    if not re.match("chr",line_dict["Chr"]):
        line_dict["Chr"]="chr"+line_dict["Chr"]
    if not (line_dict["Chr"]) in primer_dict.keys():
        primer_dict[line_dict["Chr"]]=[[line_dict["name"]],[primer_seq,],[int(line_dict["primer_start"]),],[int(line_dict["primer_end"]),]]
    else:
        primer_dict[line_dict["Chr"]][0].append(line_dict["name"])
        primer_dict[line_dict["Chr"]][1].append(primer_seq)
        primer_dict[line_dict["Chr"]][2].append(int(line_dict["primer_start"]))
        primer_dict[line_dict["Chr"]][3].append(int(line_dict["primer_end"]))
primers_file.close()
  
def process_reads_with_same_name(reads_with_same_name):
    no_mapped_r1=1
    for same_read in reads_with_same_name:
        if same_read.is_read1 and same_read.is_mapped and (not same_read.is_secondary) and (not same_read.is_supplementary):
            no_mapped_r1=0
            #check chr
            if not re.match('chr', same_read.reference_name):
                same_read_refname="chr"+same_read.reference_name
            else:
                same_read_refname=same_read.reference_name
            #find closest primer
            if same_read_refname in primer_dict.keys():
                primer_idx=0
                max_dist=args.max_dist+1
                for primer_name,primer_seq,primer_start,primer_end in zip(primer_dict[same_read_refname][0],
                                                                          primer_dict[same_read_refname][1],
                                                                          primer_dict[same_read_refname][2],
                                                                          primer_dict[same_read_refname][3]):
                    if primer_start<primer_end:
                        current_dist=abs(int(same_read.reference_start)-primer_start)
                    else:
                        current_dist=abs(int(same_read.reference_end)-primer_start)
                    if current_dist<max_dist:
                        max_dist=current_dist
                        primer_idx=primer_dict[same_read_refname][0].index(primer_name)
                if max_dist<=args.max_dist:
                    levenstain_dist=distance(same_read.get_forward_sequence()[args.skip_bp:args.skip_bp+len(primer_dict[same_read_refname][1][primer_idx])],
                                         primer_dict[same_read_refname][1][primer_idx])
                    if levenstain_dist<=args.ed_primer:
                        primer_res_dict[primer_dict[same_read_refname][0][primer_idx]].append(same_read.query_name)
                        return(1)
                    else:
                        return(0)
                else:
                    return(0)
            else:
                return(0)
    if no_mapped_r1:
        return(0)

samfile=pysam.AlignmentFile(args.bam_file, "rb")
results_file=pysam.AlignmentFile(args.output_prefix+"_filtered.bam","wb", template=samfile)
number_of_queries=0
queries_left=0
current_read_name=""
reads_with_same_name=[]

for read in samfile.fetch(until_eof=True):
    if not current_read_name:
        current_read_name=read.query_name
    else:
        if read.query_name==current_read_name:
            reads_with_same_name.append(read)
        else:
            number_of_queries+=1
            res=process_reads_with_same_name(reads_with_same_name)
            queries_left+=res
            if res:
                for selected_read in reads_with_same_name:
                    results_file.write(selected_read)
                
            current_read_name=read.query_name
            reads_with_same_name=[read,]
else:
    number_of_queries+=1
    res=process_reads_with_same_name(reads_with_same_name)
    queries_left+=res
    if res:
        for selected_read in reads_with_same_name:
            results_file.write(selected_read)
samfile.close()        
results_file.close()

print("Total number of templates:",number_of_queries)
print("Queries mapped correctly:",queries_left)

filter_file=open(args.output_prefix+"_filter.txt",'wt')

for each_primer in primer_res_dict.keys():
    print(each_primer,len(primer_res_dict[each_primer]))
    print(each_primer,(",").join(primer_res_dict[each_primer]),file=filter_file)

filter_file.close()