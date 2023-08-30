import pysam
import argparse
import scipy.stats as stats
import re

parser=argparse.ArgumentParser()

parser.add_argument("-f","--file_in",help="Input bam file name", required=True)
parser.add_argument("-b","--regions_bed",help="BED file with regions", required=True)
parser.add_argument("-o","--file_out",help="Output results file name", required=True)
parser.add_argument("-r","--reference_file",help="Reference fasta file", required=True)
parser.add_argument("-q","--filter_file",help="File with query names to filter at primer coordinates")
parser.add_argument("-x","--primers_file",help="Input file with primers")

#variant-level filtering parameters
parser.add_argument("-n","--min_depth",help="Min depth to make a call", type=int, default=5)
parser.add_argument("-a","--fract_alt",help="Min fraction of alternate allele to keep variant", type=float, default=0.03)
parser.add_argument("-m","--min_alt",help="Min reads with alternate allele to keep variant", type=int, default=5)
parser.add_argument("-c","--or_cutoff",help="Filter variants with significant strand bias at this OR", type=int, default=10)

#parameters to pass to pileup
parser.add_argument("-p","--min_map_qual",help="Min required read mappinq quality", type=int, default=10)
parser.add_argument("-s","--min_base_qual",help="Min required base quality", type=int, default=13)
parser.add_argument("-d","--max_depth",help="Max read depth permitted", type=int, default=8000)

args=parser.parse_args()

bam_chr=False
in_bamfile=pysam.AlignmentFile(args.file_in, "rb")
#print(in_bamfile.references)
if "chr1" in in_bamfile.references:
    bam_chr=True

ref_chr=False
ref_seq=pysam.FastaFile(filename=args.reference_file)
#print(ref_seq.references)
if "chr1" in ref_seq.references:
    ref_chr=True
    
filter_reads=False
if args.filter_file and args.primers_file:
    filter_file=open(args.filter_file, "rt")
    filter_reads=True
    filter_dict={}
    for line in filter_file:
        line_vals=line.strip("\n").split(" ")
        filter_dict[line_vals[0]]=line_vals[1].split(",")
    filter_file.close()
    primers_file=open(args.primers_file,'rt')
    input_header=primers_file.readline().strip("\n").split("\t")
    primer_dict={}
    for line in primers_file:
        line_dict={}
        line_vals=line.strip("\n").split("\t")
        for i in range(0,len(input_header)):
            line_dict[input_header[i]]=line_vals[i]
        if re.match("chr",line_dict["Chr"]):
            line_dict["Chr"]=line_dict["Chr"][3:]
        if int(line_dict["primer_start"])<int(line_dict["primer_end"]):
            primer_coordinates=range(int(line_dict["primer_start"])-1,int(line_dict["primer_end"]))
        else:
            primer_coordinates=range(int(line_dict["primer_end"])-1,int(line_dict["primer_start"]))
        for each_coordinate in primer_coordinates:
            primer_dict[(line_dict["Chr"],each_coordinate)]=line_dict["name"]
    primers_file.close()
    '''
    for each_key in primer_dict.keys():
        print(each_key,primer_dict[each_key])
    '''
elif args.filter_file or args.primers_file:
    print("!!!")
    print("To switch on filtering of certain bases in primer positions BOTH filter-file and primers-file are required")
    print("!!!")

def var_finder(f_pileup, coordinates, ref_dict):
    res_dict={}
    for column in f_pileup:
        coord=column.reference_pos
        ref_name=column.reference_name
        if re.match("chr",ref_name):
            ref_name=ref_name[3:]
        all_bases=column.get_query_sequences(add_indels=True)
        
        all_queries=column.get_query_names()
        
        if tuple([ref_name,column.reference_pos]) in primer_dict:
            
            queries_to_filter=filter_dict[primer_dict[(ref_name,column.reference_pos)]]
            
            for i in range(0,len(all_bases)):
                if all_queries[i] in queries_to_filter:
                    all_bases[i]="."
        
        base_dict={}
        total_bases=0
        for each_base in all_bases:
            if each_base=="*" or each_base=="<" or each_base==">" or each_base==".":
                continue
            total_bases+=1
            each_base=each_base.upper()
            if each_base in base_dict:
                base_dict[each_base]+=1  
            else:
                base_dict[each_base]=1
                
        if total_bases>=args.min_depth:
            if ref_dict[str(coord)] in base_dict.keys():
                ref_depth=base_dict[ref_dict[str(coord)]]
                del base_dict[ref_dict[str(coord)]]
            else:
                ref_depth=0
            
            for each_allele in base_dict.keys():
                if base_dict[each_allele]/total_bases>=args.fract_alt and\
                    base_dict[each_allele]>=args.min_alt:
                    ref_indices  = [index for (index, item) in enumerate(all_bases) if item.upper() == ref_dict[str(coord)]]
                    alt_indices = [index for (index, item) in enumerate(all_bases) if item.upper() == each_allele]
                    all_positions = column.get_query_positions()
                    if len(ref_indices)>=1:
                        ref_positions = [all_positions[index] for index in ref_indices]
                        alt_positions = [all_positions[index] for index in alt_indices]
                        p_value_position = stats.mannwhitneyu(ref_positions,alt_positions).pvalue
                        all_qualities = column.get_query_qualities()
                        ref_qualities = [all_qualities[index] for index in ref_indices]
                        alt_qualities = [all_qualities[index] for index in alt_indices]
                        p_value_qualities = stats.mannwhitneyu(ref_qualities,alt_qualities).pvalue
                    else:
                        p_value_position = 1
                        p_value_qualities = 1
                    #mind the strand bias
                    ADF=[str(all_bases.count(ref_dict[str(coord)])),str(all_bases.count(each_allele))]
                    ADR=[str(all_bases.count(ref_dict[str(coord)].lower())),str(all_bases.count(each_allele.lower()))]
                    odd_ratio, p_value = stats.fisher_exact([ADF,ADR])
                    if not (p_value<0.05 and (odd_ratio<=1/args.or_cutoff or odd_ratio>=args.or_cutoff)):
                        res_dict[str(coord)+","+each_allele]=[total_bases,ref_depth,base_dict[each_allele],ADF,ADR,odd_ratio,p_value,
                                                          p_value_position,p_value_qualities]
                    #print(all_bases)
    return(res_dict)

res_file=open(args.file_out,'wt')
res_file_header=["#Chr","Start","End","Ref","Alt","DP","AD_Ref","AD_Alt","ALTfraction","ADF","ADR","STRbias_OR","STRbias_p","POSbias_p","QUALbias_p"]
print("\t".join(res_file_header),file=res_file)

regions_file=open(args.regions_bed,'rt')
for line in regions_file:
    chrom,start,end=line.strip("\n").split("\t")
    if re.match("chr",chrom):
        chrom=chrom[3:]
    if ref_chr:
        reference_seq=ref_seq.fetch(reference="chr"+chrom,start=int(start), end=int(end))
    else:
        reference_seq=ref_seq.fetch(reference=chrom,start=int(start), end=int(end))
    reference_seq=reference_seq.upper()
    coordinates=list(range(int(start),int(end)))
    ref_dict={}
    for idx_coord in range(0,len(coordinates)):
        ref_dict[str(coordinates[idx_coord])]=reference_seq[idx_coord]
    if bam_chr:
        f1_pileup=in_bamfile.pileup(contig="chr"+chrom, start=int(start), end=int(end), stepper="all", truncate=True, max_depth=args.max_depth,\
                              compute_baq=False, min_mapping_quality=args.min_map_qual, min_base_quality=args.min_base_qual)
    else:
        f1_pileup=in_bamfile.pileup(contig=chrom, start=int(start), end=int(end), stepper="all", truncate=True, max_depth=args.max_depth,\
                              compute_baq=False, min_mapping_quality=args.min_map_qual, min_base_quality=args.min_base_qual)
    res_dict=var_finder(f1_pileup, coordinates, ref_dict)
    
    for res_key in res_dict.keys():
        print(chrom,":",res_key,":",res_dict[res_key])
        res_file_header=["#Chr","Start","End","Ref","Alt","DP","AD_Ref","AD_Alt","ALTfraction","ADF","ADR","STRbias_OR","STRbias_p","POSbias_p","QUALbias_p"]
        res_line_dict={}
        res_line_dict["#Chr"]=chrom
        res_line_dict["Start"],res_line_dict["Alt"]=res_key.split(",")
        res_line_dict["Ref"]=ref_dict[res_line_dict["Start"]]
        res_line_dict["Start"]=str(int(res_line_dict["Start"])+1)
        res_line_dict["End"]=res_line_dict["Start"]
        res_line_dict["DP"]=str(res_dict[res_key][0])
        res_line_dict["AD_Ref"]=str(res_dict[res_key][1])
        res_line_dict["AD_Alt"]=str(res_dict[res_key][2])
        res_line_dict["ALTfraction"]=str(round(res_dict[res_key][2]/(res_dict[res_key][1]+res_dict[res_key][2]),2))
        res_line_dict["ADF"]=",".join(res_dict[res_key][3])
        res_line_dict["ADR"]=",".join(res_dict[res_key][4])
        res_line_dict["STRbias_OR"]=str(round(res_dict[res_key][5],2))
        res_line_dict["STRbias_p"]=str(round(res_dict[res_key][6],5))
        res_line_dict["POSbias_p"]=str(round(res_dict[res_key][7],5))
        res_line_dict["QUALbias_p"]=str(round(res_dict[res_key][8],5))
        
        #find ref
        #format insertions and deletions
        '''
        for i in range(0,len(res_file_header)):
            print(res_file_header[i]+": "+res_line_dict[res_file_header[i]])
        '''    
        if re.search('\-\d+N+',res_line_dict["Alt"]):
            #print("Deletion found")
            #print("Formatting...")
            num_del=int((re.search('\d+',res_line_dict["Alt"])).group())
            res_line_dict["Start"]=int(res_line_dict["Start"])+1
            res_line_dict["Ref"]=""
            for k in range(res_line_dict["Start"],res_line_dict["Start"]+num_del):
                res_line_dict["Ref"]+=ref_dict[str(k-1)]
            res_line_dict["Alt"]="-"
            res_line_dict["End"]=str(res_line_dict["Start"]+num_del-1)
            res_line_dict["Start"]=str(res_line_dict["Start"])
            '''
            for i in range(0,len(res_file_header)):
                print(res_file_header[i]+": "+res_line_dict[res_file_header[i]])
            '''
        elif re.search('\+\d+',res_line_dict["Alt"]):
            #print("Insertion found")
            #print("Formatting...")
            res_line_dict["Ref"]="-"
            res_line_dict["Alt"]=re.sub('\D\+\d+','',res_line_dict["Alt"])
            
        res_to_print=[]
        for i in range(0,len(res_file_header)):
            res_to_print.append(res_line_dict[res_file_header[i]])
        print("\t".join(res_to_print),file=res_file)
res_file.close()    
in_bamfile.close()
ref_seq.close()