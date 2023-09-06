# RACE_NGS
This collection of python scripts was developed to help in the design of the targeted RNA panels on the basis of anhored multiplex PCR-based methods (primarily, 3' and 5' rapid anplification of cDNA ends, or RACE), and for the bioinformatic analysis of the obtained data. While the analysis of rearrangements from such data can effectively be done with the available software (e.g., the STAR-fusion https://github.com/STAR-Fusion/STAR-Fusion), analysis of somatic mutations and gene expression can be more challenging in such datasets. Thus, additional scripts were created for these two purposes.

FindRACEprimers.py 

This is the main script used for the primer design for 3'RACE or 5'RACE-based approaches. The input coordinates file in .tsv format should contain at least 3 necessary fields: "Chr"	"Start"	"End" (1-based coordinates are required; if the gene is located on a 'minus' strand, the first coordinate should be smaller than the second one!!!). It is recomended to provide the coordinates for the whole exons, if not too long, for the convenience of the following data analysis. For the 3'RACE select --primer_strand as "forward" and for the 5'RACE select --primer_strand as "reverse". Set the desired primers' melting temperature with --melt_temp parameter. The following arguments are required for the calculation of melting temperatures of potential primers with the SantaLucia method (SantaLucia, 1998), with the corrections, suggested by von Ahsen et al. (1999): 

  --conc_primer, Primer concentration in PCR, nM;
  
  --conc_mg, Mg concentration in PCR;
  
  --conc_salt, Salt concentration in PCR, mM;
  
The code for the calculation of melting temperature is adapted from http://www.biophp.org/minitools/melting_temperature. 

Other required parameters:

  --ref_file, Genomic reference sequence in FASTA format;
  
  --vcf_file, Reference SNP database file (the gnomAD v.3 database https://gnomad.broadinstitute.org/ can be used, but only as a single vcf file);
  
  --spec_fasta, File with reference transcripts in FASTA format (this database is required to check for the primers' specificity, the RefSeq transcripts database https://www.ncbi.nlm.nih.gov/refseq/ is recommended);

All optional parameters can be viewed using --help option.

The output file with default name 'primer_list.tsv' contains the following fields:

FILTER: this field contains the characteristics of primers, which are not optimal (e.g. SPEC means specificity issues, CG - too low or too high CG-content, etc.);

all the fields, which were present in the input file with genomic coordinates;

primer_start, primer_end: primer's coordinates;	

primer_seq: primer's sequence;

primer_CG: primer's CG-content;

primer_Tm: primer's melting temperature;

primer_len: primer's length;

SNP: SNPs, with the population frequencies above the threshold (default 0.01), within primer coordinates; 

SNP_3position: the position of SNP within primer's sequence, starting from the 3' end; 

SNP_info: the database SNP information;

Ngenes_0mis, Ngenes_1mis, Ngenes2mis, Ngenes_3mis: number of genes, whose transcripts contain regions, where the primer can hybridize with 0, 1, 2 or 3 mismatches, respectively;

Genes_0mis, Genes_1mis, Genes_2mis, Genes_3mis: the respective genes.

Regarding the check of specificity, it is supposed that the Taq polymerase lacking 3'-error correction activity is used in reaction; thus, the targets with 3'end mismatches with the primer are ignored.



Check_rRNA.py

This script can be used to check the selected primers' possibility to hybridize with ribosomal RNAs. The input file in .tsv format should contain at least the 'primer_seq' field. File with reference transcripts in FASTA format, specified by --spec_fasta, is required and should contain non-abbreviated gene names. The minimum number of mismatches, which is considered safe, is specified with --mis_num parameter. The alignments with rRNA species, which have fewer number of mismatches will be outputted. 


Check_RACEprimers.py

This is a script written for the analysis of the results obtained from pilot experiments. It requires R1.fastq file only (provided via --fastq_file option). The .tsv file with information about primers used for the enrichment and the genomic coordinates of the exons, where the primers are located, also should be provided. Normally, the file generated with FindRACEprimers.py, is supposed to be used. The required columns are "Start", "End" (which are exon coordinates), "primer_start", "primer_end" (which are primer coordinates), "primer_seq" (primer sequence) and "Gene" (gene name). Additional columns, e.g. "primer_name" can be added, and all input columns will be outputted. The program searches for the primer sequence in the beginning of a Read 1, then, if there is enough distance left till the exon end, the sequence in the Read 1, which immediately follows the primer sequence, is compared to the reference sequence (the reference fasta file should be provided with --ref_file argument). The other required arguments are --primer_strand ("forward" for 3' RACE and "reverse" for 5' RACE methods) and the --sample_name. The program counts reads, where the expected reference sequence follows primer sequence, as "specific". In the outputted file with default name "Count_reads_per_primer.tsv", the following new columns will appear:

Sample: the sample name, which is useful then the file containing multiple samples is created;

ref_to_compare: the reference sequence, within the exonic coordinates, which was compared to the actual Read 1 sequence (if the primer's location is too close to the exon border, this field can be empty);

all_reads: the total number of reads, starting from the primer sequence;

spec_reads: the number of reads, in which the exprected reference sequence follows the primer sequence;

nonspec_reads: the number of reads, in which the exprected reference sequence does not follow the primer sequence;

primer_eff: primer efficiency, calculated relative to the primer with maximum number of reads, among all primers located in the same gene;

primer_spec: primer specificity, calculated as the percentage of "specific" reads relative to the total number of reads.

Additional arguments can be accessed with the --help option.

Filter_RACEbam.py

The script prepares deduplicated and aligned (presumably, with the STAR RNA-seq aligner https://github.com/alexdobin/STAR) .bam files for the following mutation calling and expression analysis. The input .bam file should be sorted by queryname!!! Input .tsv file with primers has to be provided using --primers_file option. This is normally the file generated by the FindRACEprimers.py, where only rows with primers included in the panel are left, but with at least one additional column "name" with the user-defined names for the primers. At least "Chr", "primer_start", "primer_end" and "primer_seq" columns should be kept from the original file generated by FindRACEprimers.py. The last required option is the prefix for the output files (e.g., sample name).

The program outputs the filtered .bam file, where only those templates are kept, which were mapped to the coordinates, predefined by the primer, whose sequence is found in the Read 1. The "_filter.txt" contains the lists of the respective reads for each primer and is used in the downstream analysis.

Additional arguments can be accessed with the --help option.

RACE_Caller.py

This script was created for the variant calling. The sensitivity is defined by the minimun fraction of the alternate allele (--fract_alt option with default value 0.03) and the minimum number of reads suppoting the alternate allele (--min_alt option with default value 5). The sorted by coordinate and indexed .bam file, filtered with Filter_RACEbam.py program, together with the "_filter.txt" file should be provided using options --file_in and --filter_file, respectively. Also, a file in BED format, specifying regions, where variant calling should be performed, is required, as well as the file with primers used for enrichment. These files need to be provided using --regions_bed and --primers_file options, respectively. The other necessary options are --reference_file (reference fasta file) and --file_out (the output file name). Additional arguments can be accessed with the --help option.  

RACE_Caller.py outputs results in the table format, compatible with ANNOVAR (https://annovar.openbioinformatics.org). The following columns are included:

#Chr,Start,End,Ref,Alt: as defined in ANNOVAR documentation (https://annovar.openbioinformatics.org/en/latest/#annovar-documentation);

DP: depth of coverage;

AD_Ref: allelic depth (reference allele);

AD_Alt: allelic depth (alternate allele);

ALTfraction: alternate allele fraction;

ADF: allelic depths (forward strand);

ADR: allelic depths (reverse strand);

STRbias_OR: odds ratio for the strand bias;

STRbias_p: Fisher's exact test p-value for the strand bias;

POSbias_p: Mann-Whitney U-test p-value for the position bias;

QUALbias_p: Mann-Whitney U-test p-value for the base quality bias.

Expression_count.py

This script was written for the analysis of gene expression at the specified genomic regions. The regions should be provided in a .tsv file, formatted as the example file Expression_count_example_input.tsv, specified using the --coordinates option. The necessary fields are "Primer" (the user-defined names for the primers, the same as were used with Filter_RACEbam.py), "Chr", "Start", "End", and "Name" (the designations for selected regions). In the calculation of coverage within specified regions only reads, generated by the predefined primers are included. The sorted by coordinate and indexed .bam file, filtered with Filter_RACEbam.py program, together with the "_filter.txt" file should be provided using options --file_in and --filter_file, respectively. The --sample_name argument is also required.

If the program is called multiple times (one time for each sample), a single file with the default name "Expression_results.tsv" will be created. This file will contain the calculated coverage (expression) values for each sample in each of the selected regions. 

In case of any problems or questions regarding the programs in this repository, please contact Natalia V. Mitiushkina by e-mail nmmail@inbox.ru.


![image](https://github.com/MitiushkinaNV/RACE_NGS/assets/96590759/372639ee-f20a-48fa-bc00-1783ea5a2a9f)

