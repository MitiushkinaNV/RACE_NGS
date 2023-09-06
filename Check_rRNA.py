import pysam
import re
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import parse
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


parser=argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="Input file with primer_seq", required=True)
parser.add_argument("-s","--primer_strand", help="Select either 'forward' for 3'RACE or 'reverse' for 5'RACE method", required=True)
parser.add_argument("-l","--tail_seq",help="5'-tail of each primer, if present", default="")
parser.add_argument("-o","--output_file",help="Output file name")
parser.add_argument("-a","--spec_fasta",help="File with reference transcripts in FASTA format", required=True)
parser.add_argument("-m","--mis_num",help="Minimum number of mismatches considered safe", type=int, default=5)

args=parser.parse_args()

if args.primer_strand!="forward" and args.primer_strand!="reverse":
    print("The primer_strand argument should be either forward (for 3' RACE) or reverse (for 5' RACE)")
    exit()

spec_file=open(args.spec_fasta)
spec_records=parse(spec_file,"fasta")

ribo_mt_dict={}
for spec_record in spec_records:
    if re.search("PREDICTED:",spec_record.description):
        continue
    if re.search("ribosomal RNA$",spec_record.description):
        ribo_mt_dict[spec_record.id]=((spec_record.description)[len(spec_record.id)+1:],str(spec_record.seq))

input_file=open(args.input_file,'rt')
if args.output_file:
    output_file=open(args.output_file,'wt')
input_header=input_file.readline().strip("\n").split("\t")

for line in input_file:
    line_list=line.strip("\n").split("\t")
    line_dict={}
    for i in range(0,len(input_header)):
        line_dict[input_header[i]]=line_list[i]
        if args.output_file:
            print(input_header[i],line_dict[input_header[i]], file=output_file)
        else:
            print(input_header[i],line_dict[input_header[i]])
    primer_seq=Seq((line_dict["primer_seq"])[len(args.tail_seq):])
    if args.output_file:
        print("Primer seq without tail",str(primer_seq), file=output_file)
    else:
        print("Primer seq without tail",str(primer_seq))
        
    for ribo_mt_id in ribo_mt_dict.keys():
        ribo_mt_name,ribo_mt_seq=ribo_mt_dict[ribo_mt_id]
        if args.output_file:
            print("checking primer vs",ribo_mt_name, file=output_file)
        else:
            print("checking primer vs",ribo_mt_name)
        
        ribo_mt_seq=Seq(ribo_mt_seq)
        if args.primer_strand=="reverse":
            ribo_mt_seq=ribo_mt_seq.reverse_complement()
        a_found=0
        for a in pairwise2.align.globalms(primer_seq,ribo_mt_seq,1,0,-10, 0, penalize_end_gaps=(False,False)):
            
            if a.score>(len(str(primer_seq))-args.mis_num):
                if args.output_file:
                    print(format_alignment(*a),file=output_file)
                else:
                    print(format_alignment(*a))
                a_found=1
        if a_found==0:
            if args.output_file:
                print("OK",file=output_file)
            else:
                print("OK")
spec_file.close()
if args.output_file:
    output_file.close()
input_file.close()
