import pysam
import re
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import parse
import ahocorasick
import math

parser=argparse.ArgumentParser()
parser.add_argument("-c","--coordinates_file",help="Input file with coordinates", required=True)
parser.add_argument("-r","--ref_file",help="Genomic reference sequence", required=True)
parser.add_argument("-v","--vcf_file",help="Reference SNP database file", required=True)
parser.add_argument("-f","--vcf_freq",help="Threshold for polymorphism frequency", type=float, default=0.01)
parser.add_argument("-s","--primer_strand", help="Select primers on forward or on reverse strand", required=True)
parser.add_argument("-t","--melt_temp",help="Melting temperature", required=True, type=float)
parser.add_argument("-l","--tail_seq",help="5'-tail to add to each primer", default="")
parser.add_argument("-o","--output_file",help="Output file name", default="primer_list.txt")
parser.add_argument("-a","--spec_fasta",help="File to be used in check of primer specificity", required=True)
parser.add_argument("-p","--conc_primer",help="Primer concentration, nM", required=True, type=float)
parser.add_argument("-m","--conc_mg",help="Mg concentration, mM", required=True, type=float)
parser.add_argument("-b","--conc_salt",help="Salt concentration, mM", required=True, type=float)

parser.add_argument("-i","--limit_genes",help="Limit number of outputted gene names", type=int, default=10)
parser.add_argument("-n","--min_length",help="Min primer length", type=int, default=15)

args=parser.parse_args()

if args.primer_strand!="forward" and args.primer_strand!="reverse":
    print("The primer_strand argument should be either forward (for 3' RACE) or reverse (for 5' RACE)")
    exit()

#from table at http://www.ncbi.nlm.nih.gov/pmc/articles/PMC19045/table/T2/ (SantaLucia, 1998)
#enthalpy values
enthalpy_dict={}
enthalpy_dict["AA"]= -7.9
enthalpy_dict["AC"]= -8.4
enthalpy_dict["AG"]= -7.8
enthalpy_dict["AT"]= -7.2
enthalpy_dict["CA"]= -8.5
enthalpy_dict["CC"]= -8.0
enthalpy_dict["CG"]=-10.6
enthalpy_dict["CT"]= -7.8
enthalpy_dict["GA"]= -8.2
enthalpy_dict["GC"]= -9.8
enthalpy_dict["GG"]= -8.0
enthalpy_dict["GT"]= -8.4
enthalpy_dict["TA"]= -7.2
enthalpy_dict["TC"]= -8.2
enthalpy_dict["TG"]= -8.5
enthalpy_dict["TT"]= -7.9
entropy_dict={}
entropy_dict["AA"]=-22.2
entropy_dict["AC"]=-22.4
entropy_dict["AG"]=-21.0
entropy_dict["AT"]=-20.4
entropy_dict["CA"]=-22.7
entropy_dict["CC"]=-19.9
entropy_dict["CG"]=-27.2
entropy_dict["CT"]=-21.0
entropy_dict["GA"]=-22.2
entropy_dict["GC"]=-24.4
entropy_dict["GG"]=-19.9
entropy_dict["GT"]=-22.4
entropy_dict["TA"]=-21.3
entropy_dict["TC"]=-22.2
entropy_dict["TG"]=-22.7
entropy_dict["TT"]=-22.2

def define_Tm(primer_seq, args, enthalpy_dict, entropy_dict):
    h=0
    s=0
    
    #effect on entropy by salt correction; von Ahsen et al 1999
    #Increase of stability due to presence of Mg;
    salt_effect= (args.conc_salt/1000)+((args.conc_mg/1000)*140)
    #effect on entropy
    s+=0.368*(len(primer_seq)-1)*math.log(salt_effect)
    
    #terminal corrections. Santalucia 1998
    firstnucleotide=primer_seq[0]
    if (firstnucleotide=="G" or firstnucleotide=="C"):
        h+=0.1
        s+=(-2.8)
    elif (firstnucleotide=="A" or firstnucleotide=="T"):
        h+=2.3
        s+=4.1

    lastnucleotide=primer_seq[len(primer_seq)-1]
    if (lastnucleotide=="G" or lastnucleotide=="C"):
        h+=0.1
        s+=(-2.8)
    if (lastnucleotide=="A" or lastnucleotide=="T"):
        h+=2.3
        s+=4.1
    
    #compute new H and s based on sequence. Santalucia 1998
    for i in range(0,len(primer_seq)-1):
        subc=primer_seq[i:i+2]
        h+=enthalpy_dict[subc]
        s+=entropy_dict[subc]
    primer_tm=((1000*h)/(s+(1.987*math.log(args.conc_primer/2000000000))))-273.15
    return(primer_tm)


ref_seq=pysam.FastaFile(filename=args.ref_file)
ref_chr=False
if "chr1" in ref_seq.references:
    ref_chr=True
vcf_file=pysam.VariantFile(filename=args.vcf_file)
coordinates_file=open(args.coordinates_file,'rt')
results_file=open(args.output_file,'wt')
input_header=coordinates_file.readline().strip("\n\t\s").split("\t")
header=["FILTER"]+input_header+["primer_start","primer_end","primer_seq","primer_CG","primer_Tm","primer_len"]
header+=["SNP","SNP_3position","SNP_info"]
header+=["Ngenes_0mis","Genes_0mis","Ngenes_1mis","Genes_1mis","Ngenes2mis","Genes_2mis","Ngenes_3mis",\
         "Genes_3mis"]
print("\t".join(header), file=results_file)

spec_file=open(args.spec_fasta)
spec_records=parse(spec_file,"fasta")
spec_dict={}
for spec_record in spec_records:
    if re.search("PREDICTED:",spec_record.description):
        continue
    spec_gene=(re.split(", transcript variant",spec_record.description))[0]
    spec_dict[spec_record.id]=((spec_gene)[len(spec_record.id)+1:],str(spec_record.seq))
vcf_file = pysam.VariantFile(args.vcf_file)
if "chr1" in vcf_file.header.contigs:
    vcf_chr=True

for line in coordinates_file:
    print(line)
    line_dict={}
    line_vals=line.strip("\n\t\s").split("\t")
    for i in range(0,len(input_header)):
        line_dict[input_header[i]]=line_vals[i]
    if ref_chr and not(re.match("chr",line_dict["Chr"])):
        chrom="chr"+line_dict["Chr"]
    elif not ref_chr and re.match("chr",line_dict["Chr"]):
        chrom=line_dict["Chr"][3:]
    else:
        chrom=line_dict["Chr"]
    if int(line_dict["Start"])<int(line_dict["End"]):
        dna_strand="plus"
        pos_start=int(line_dict["Start"])
        pos_end=int(line_dict["End"])
    else:
        dna_strand="minus"
        pos_start=int(line_dict["End"])
        pos_end=int(line_dict["Start"])
    
    reference_seq=ref_seq.fetch(reference=chrom,start=pos_start-1,end=pos_end)
    reference_seq=Seq(reference_seq.upper())
    if (dna_strand=="plus" and args.primer_strand=="reverse") or\
        (dna_strand=="minus" and args.primer_strand=="forward"):
        reference_seq=reference_seq.reverse_complement()
        position_vector=[i for i in range(pos_end,pos_start-1,-1)]
    else:
        position_vector=[i for i in range(pos_start,pos_end+1)]
    reference_seq=str(reference_seq)
    print(reference_seq)
    if (vcf_chr and ref_chr) or (not(vcf_chr) and not(ref_chr)):
        vcf_records=vcf_file.fetch(contig=chrom,start=pos_start-1, stop=pos_end)
    elif vcf_chr and not(ref_chr):
        vcf_records=vcf_file.fetch(contig="chr"+chrom,start=pos_start-1, stop=pos_end)
    elif not(vcf_chr) and ref_chr:
        vcf_records=vcf_file.fetch(contig=chrom[3:],start=pos_start-1, stop=pos_end)
    vcf_alleles={}
    vcf_id={}
    vcf_AF={}
    vcf_maxAF={}
    vcf_maxPOP={}
    vcf_type={}
    for vcf_record in vcf_records:
        #find frequencies
        freq=[]
        freq_names=[]
        for info in vcf_record.info.keys():
            if re.match('AF',info):
                freq.append(vcf_record.info[info][0])
                freq_names.append(info)
        maxAF=max(freq)
        maxPOP_indexes=[i for i, e in enumerate(freq) if e == maxAF]
        maxPOP=""
        for maxPOP_index in maxPOP_indexes:
            maxPOP+=freq_names[maxPOP_index]
        if vcf_record.id:
            vcf_rec_id=vcf_record.id
        else:
            vcf_rec_id='.'
        if len(vcf_record.alleles[0])==len(vcf_record.alleles[1]):
            var_type="SNP"
            position_list=[vcf_record.pos]
        elif len(vcf_record.alleles[0])>len(vcf_record.alleles[1]):
            var_type="DEL"
            position_list=[i for i in range(vcf_record.pos+1,vcf_record.pos+1+(len(vcf_record.alleles[0])-len(vcf_record.alleles[1])))]
        elif len(vcf_record.alleles[0])<len(vcf_record.alleles[1]):
            var_type="INS"
            position_list=[vcf_record.pos]
        for position in position_list:
            if position in vcf_alleles:
                vcf_alleles[position].append(">".join(vcf_record.alleles))
                vcf_id[position].append(vcf_rec_id)
                vcf_AF[position].append(freq[0])
                vcf_maxAF[position].append(maxAF)
                vcf_maxPOP[position].append(maxPOP)
                vcf_type[position].append(var_type)
            else:
                vcf_alleles[position]=[">".join(vcf_record.alleles)]
                vcf_id[position]=[vcf_rec_id]
                vcf_AF[position]=[freq[0]]
                vcf_maxAF[position]=[maxAF]
                vcf_maxPOP[position]=[maxPOP]
                vcf_type[position]=[var_type]
    SNP_type={}
    SNP_info={}
    for position in vcf_alleles.keys():
        sumAF=sum(vcf_AF[position])
        if max(sumAF,max(vcf_maxAF[position]))>=args.vcf_freq:
            if len(vcf_type[position])>1:
                SNP_type[position]="MULT"
            else:
                if vcf_type[position][0]=="SNP":
                    SNP_type[position]="SNP"
                else:
                    SNP_type[position]="INDEL"
            info=[]
            for i in range(0,len(vcf_alleles[position])):
                info.append(":".join([vcf_id[position][i],vcf_alleles[position][i],str(round(vcf_AF[position][i],4)),
                                  vcf_maxPOP[position][i],str(round(vcf_maxAF[position][i],4))]))
            SNP_info[position]="/".join(info)
    
    #create primers
    
    print_res=[]
    for i in range(0,len(input_header)):
        print_res.append(line_dict[input_header[i]])
    primer_end=len(reference_seq)-1
    while True:
        primer_seq=""
        primer_start=primer_end
        primer_CG=0
        primer_Tm=0
        primer_len=0
        #collect nucleotides
        primer_SNP_3positions=[]
        primer_SNP_types=[]
        primer_SNP_infos=[]
        while primer_Tm<args.melt_temp and primer_start>=0:
            primer_seq=reference_seq[primer_start]+primer_seq
            primer_len+=1
            if(primer_len>=args.min_length):
                primer_Tm=define_Tm(primer_seq, args, enthalpy_dict, entropy_dict)
            if reference_seq[primer_start]=="C" or reference_seq[primer_start]=="G":
                primer_CG+=1
            if position_vector[primer_start] in SNP_type:
                primer_SNP_3positions.append(len(primer_seq))
                primer_SNP_types.append(SNP_type[position_vector[primer_start]])
                primer_SNP_infos.append(":".join([str(len(primer_seq)),SNP_info[position_vector[primer_start]]]))
            primer_start-=1
        if primer_Tm>=args.melt_temp:
            #check if any snp in primers
            filt_res=[]
            if primer_CG/primer_len<0.4 or primer_CG/primer_len>0.6:
                filt_res.append("CG")
            if primer_SNP_3positions:
                filt_res.append("SNP")
                primer_SNP_3position=str(min(primer_SNP_3positions))
                primer_SNP_info=";".join(primer_SNP_infos)
                if len(primer_SNP_types)>1:
                    primer_SNP_type="MULT"
                else:
                    primer_SNP_type=primer_SNP_types[0]
            else:
                primer_SNP_3position=""
                primer_SNP_info=""
                primer_SNP_type=""
            
            #mutate primer
            
            print("Primer: ", primer_seq)
            
            A = ahocorasick.Automaton()
            A.add_word(primer_seq, 0)
            
            for i in range(0,len(primer_seq)-1):
                for nucl1 in ["A","T","G","C"]:
                    if primer_seq[i]!=nucl1:
                        mut_seq1=primer_seq[0:i]+nucl1+primer_seq[(i+1):len(primer_seq)]
                        
                        A.add_word(mut_seq1, 1)
                        
                        for j in range(0,len(primer_seq)-1):
                            if (i!=j) and not(i>(len(primer_seq)-5) and j>(len(primer_seq)-5)):
                                for nucl2 in ["A","T","G","C"]:
                                    if primer_seq[j]!=nucl2:
                                        mut_seq2=mut_seq1[0:j]+nucl2+mut_seq1[(j+1):len(primer_seq)]
                                        
                                        A.add_word(mut_seq2, 2)
                                        
                                        for k in range(0,len(primer_seq)-1):
                                            if (k!=i and k!=j) and\
                                                not(i>(len(primer_seq)-6) and j>(len(primer_seq)-6)\
                                                and k>(len(primer_seq)-6) and \
                                                not(k>(len(primer_seq)-5) and j>(len(primer_seq)-5)) and \
                                                not(i>(len(primer_seq)-5) and k>(len(primer_seq)-5))):
                                                for nucl3 in ["A","T","G","C"]:
                                                    if primer_seq[k]!=nucl3:
                                                        mut_seq3=mut_seq2[0:k]+nucl3+mut_seq2[(k+1):len(primer_seq)]
                                                        
                                                        A.add_word(mut_seq3, 3)
                                                        
            A.make_automaton()
            
            spec_mut_dict={}
            spec_mut_dict[0]=[]
            spec_mut_dict[1]=[]
            spec_mut_dict[2]=[]
            spec_mut_dict[3]=[]
            for spec_id in spec_dict.keys():
                spec_description,spec_seq=spec_dict[spec_id]
                if args.primer_strand=="forward":
                    Ares=A.iter(str(spec_seq))
                elif args.primer_strand=="reverse":
                    Ares=A.iter("NNN"+str(Seq(spec_seq[3:]).reverse_complement()))
                for each_Ares in Ares:
                    #print("Ares",spec_id,spec_description,each_Ares)
                    spec_mut_dict[each_Ares[1]].append(spec_description)
            n_target=0
            for i in range(0,4):
                spec_mut_dict[i]=list(set(spec_mut_dict[i]))
                n_target+=len(spec_mut_dict[i])
                if len(spec_mut_dict[i])>args.limit_genes:
                    spec_mut_dict[i]=["MORE than "+str(args.limit_genes)+" genes",]
            if n_target>1:
                filt_res.append("SPEC")
                
            primer_seq=args.tail_seq+primer_seq
            to_print=print_res+[str(position_vector[primer_start+1]),str(position_vector[primer_end]),\
                                primer_seq,str(round(primer_CG/primer_len,2)),str(primer_Tm),\
                                str(primer_len),primer_SNP_type,primer_SNP_3position,primer_SNP_info,\
                                str(len(spec_mut_dict[0])),";".join(spec_mut_dict[0]),\
                                str(len(spec_mut_dict[1])),";".join(spec_mut_dict[1]),\
                                str(len(spec_mut_dict[2])),";".join(spec_mut_dict[2]),\
                                str(len(spec_mut_dict[3])),";".join(spec_mut_dict[3])
                                ]
            print(";".join(filt_res)+"\t"+"\t".join(to_print),file=results_file)
        if primer_start+1==0:
            break
        primer_end-=1
        
ref_seq.close()
vcf_file.close()
coordinates_file.close()
results_file.close()