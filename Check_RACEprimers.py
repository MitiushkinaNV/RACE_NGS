import pysam
import argparse
from Bio.Seq import Seq
import gzip
from Levenshtein import distance

parser=argparse.ArgumentParser()

parser.add_argument("-r","--ref_file",help="Genomic reference sequence", required=True)
parser.add_argument("-p","--primers_file",help="Input file with primers", required=True)
parser.add_argument("-f","--fastq_file",help="Input R1 fastq file", required=True)
parser.add_argument("-b","--skip_bp",help="Skip this number of bp from read start", default=0, type=int)
parser.add_argument("-e","--ed_primer",help="Max edit distance with primer", default=2, type=int)
parser.add_argument("-d","--ed_ref",help="Max edit distance with reference", default=2, type=int)
parser.add_argument("-i","--compare_min",help="Min number of bp to compare with reference", default=10, type=int)
parser.add_argument("-a","--compare_max",help="Max number of bp to compare with reference", default=20, type=int)
parser.add_argument("-l","--tail_seq",help="5'-tail to remove from each primer", default="")
parser.add_argument("-o","--output_file",help="Output file name", default="specificity_expression_res.txt")
parser.add_argument("-s","--primer_strand", help="Were primers selected on forward or on reverse strand", required=True)
parser.add_argument("-n","--sample_name",help="Sample name", required=True)

args=parser.parse_args()

if args.primer_strand!="forward" and args.primer_strand!="reverse":
    print("The primer_strand argument should be either forward (for 3' RACE) or reverse (for 5' RACE)")
    exit()

ref_seq=pysam.FastaFile(filename=args.ref_file)
primers_file=open(args.primers_file,'rt')
results_file=open(args.output_file,'at')

input_header=primers_file.readline().strip("\n").split("\t")
header=["Sample"]+input_header+["ref_to_compare","all_reads","spec_reads","nonspec_reads",
                                "primer_eff","primer_spec"]
print("\t".join(header), file=results_file)

primer_line_dict={}
primer_res_dict={}
primer_ref_dict={}

primer_line_dict["unknown"]=[""]*len(input_header)
primer_line_dict["unknown"][0]="unknown"
primer_res_dict["unknown"]=[0,0,0]
primer_ref_dict["unknown"]=""

for line in primers_file:
    line_dict={}
    line_vals=line.strip("\n").split("\t")
    for i in range(0,len(input_header)):
        line_dict[input_header[i]]=line_vals[i]
    primer_seq=(line_dict["primer_seq"])[len(args.tail_seq):]
    primer_line_dict[primer_seq]=line_vals
    primer_res_dict[primer_seq]=[0,0,0]
    if (args.primer_strand=="forward" and \
        abs(int(line_dict["End"])-int(line_dict["primer_end"]))>=args.compare_min) or\
       (args.primer_strand=="reverse" and \
        abs(int(line_dict["Start"])-int(line_dict["primer_end"]))>=args.compare_min):
        chrom=line_dict["Chr"]
        #chrom="chr"+line_dict["Chr"]
        if args.primer_strand=="forward":
            ref_start=min(int(line_dict["End"]),int(line_dict["primer_end"]))
            ref_end=max(int(line_dict["End"]),int(line_dict["primer_end"]))
        elif args.primer_strand=="reverse":
            ref_start=min(int(line_dict["Start"]),int(line_dict["primer_end"]))
            ref_end=max(int(line_dict["Start"]),int(line_dict["primer_end"]))
        if ((args.primer_strand=="reverse" and int(line_dict["Start"])<int(line_dict["End"])) or\
            (args.primer_strand=="forward" and int(line_dict["Start"])>int(line_dict["End"]))):
            ref_start-=1
            ref_end-=1
        if ref_end-ref_start>args.compare_max:
            if ((args.primer_strand=="reverse" and int(line_dict["Start"])<int(line_dict["End"])) or\
                (args.primer_strand=="forward" and int(line_dict["Start"])>int(line_dict["End"]))):
                ref_start=ref_end-args.compare_max
            else:
                ref_end=ref_start+args.compare_max
        ref_to_compare=ref_seq.fetch(reference=chrom,\
            start=ref_start,
            end=ref_end)
        if (args.primer_strand=="reverse" and int(line_dict["Start"])<int(line_dict["End"])) or\
            (args.primer_strand=="forward" and int(line_dict["Start"])>int(line_dict["End"])):
            ref_to_compare=Seq(ref_to_compare.upper())
            ref_to_compare=ref_to_compare.reverse_complement()
            ref_to_compare=str(ref_to_compare)
        primer_ref_dict[primer_seq]=ref_to_compare
    else:
        primer_ref_dict[primer_seq]=""
        
primers_file.close()

fastq_file=gzip.open(args.fastq_file,'rt')

n=0

for line in fastq_file:
    if n==1:
        line=line[args.skip_bp:]
        line.strip()
        line_primer="unknown"
        for primer in primer_ref_dict.keys():
            if len(line)>=len(primer):
                if primer==line[:len(primer)]:
                    line_primer=primer
                    break
        if line_primer=="unknown":
            #check with edit distance
            for primer in primer_ref_dict.keys():
                if len(line)>=len(primer):
                    if distance(primer,line[:len(primer)])<=args.ed_primer:
                        line_primer=primer
                        break
        primer_res_dict[line_primer][0]+=1
        #check specificity
        if line_primer=="unknown" or not(primer_ref_dict[line_primer]):
            pass
        elif len(line)>=len(line_primer)+len(primer_ref_dict[line_primer]):
            if line[len(line_primer):len(line_primer)+len(primer_ref_dict[line_primer])]==\
                primer_ref_dict[line_primer]:
                primer_res_dict[line_primer][1]+=1
            elif distance(line[len(line_primer):len(line_primer)+len(primer_ref_dict[line_primer])],
                primer_ref_dict[line_primer])<=args.ed_ref:
                primer_res_dict[line_primer][1]+=1
            else:
                primer_res_dict[line_primer][2]+=1
        elif len(line)>=len(line_primer)+args.compare_min:
            if line[len(line_primer):len(line_primer)+args.compare_min]==\
                primer_ref_dict[line_primer][:args.compare_min]:
                primer_res_dict[line_primer][1]+=1
            elif distance(line[len(line_primer):len(line_primer)+args.compare_min],
                primer_ref_dict[line_primer][:args.compare_min])<=args.ed_ref:
                primer_res_dict[line_primer][1]+=1
            else:
                primer_res_dict[line_primer][2]+=1
    n+=1
    if n>3:
        n=0


print("Primer sequence was not found in ", primer_res_dict["unknown"][0], "reads")
no_primer=primer_res_dict["unknown"][0]
del primer_line_dict["unknown"]
del primer_res_dict["unknown"]
del primer_ref_dict["unknown"]

total_reads=0
total_reads_specific=0
total_reads_unspecific=0

gene_index=input_header.index("Gene")
gene_reads_dict={}
for each_primer in primer_line_dict.keys():
    gene_name=primer_line_dict[each_primer][gene_index]
    if gene_name in gene_reads_dict:
        gene_reads_dict[gene_name].append(primer_res_dict[each_primer][1])
    else:
        gene_reads_dict[gene_name]=[primer_res_dict[each_primer][1]]
    


for each_primer in primer_line_dict.keys():
    gene_name=primer_line_dict[each_primer][gene_index]
    if primer_res_dict[each_primer][1]+primer_res_dict[each_primer][2]>0:
        primer_spec=round(100*(primer_res_dict[each_primer][1])/(primer_res_dict[each_primer][1]\
                                                             +primer_res_dict[each_primer][2]),1)
    else:
        primer_spec=""
    if len(gene_reads_dict)>1 and max(gene_reads_dict[gene_name])>0:
        primer_eff_max=round(100*(primer_res_dict[each_primer][1])/(max(gene_reads_dict[gene_name])),1)
    else:
        primer_eff_max=""
    
    result_line=[args.sample_name]+primer_line_dict[each_primer]+[primer_ref_dict[each_primer]]\
        +list(map(str, primer_res_dict[each_primer]))+list(map(str, [primer_eff_max,primer_spec]))
    print("\t".join(result_line), file=results_file)
    total_reads+=primer_res_dict[each_primer][0]
    total_reads_specific+=primer_res_dict[each_primer][1]
    total_reads_unspecific+=primer_res_dict[each_primer][2]
    
print("Reads with primer mapped: ", total_reads, "(",round(total_reads*100/(total_reads+no_primer),1),"%)")
if total_reads_specific+total_reads_unspecific>0:
    print("Specific annealing: ", total_reads_specific, "(",\
      round(total_reads_specific*100/(total_reads_specific+total_reads_unspecific),1),"%)")
print("Unspecific reads: ", total_reads_unspecific)
    
results_file.close()
ref_seq.close()

    
