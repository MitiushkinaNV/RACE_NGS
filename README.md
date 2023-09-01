# RACE_NGS
This collection of python scripts was developed to help in the design of the targeted RNA panels on the basis of anhored multiplex PCR-based methods (primarily, 3' and 5' rapid anplification of cDNA ends or RACE), and for the following analysis of the obtained data. While the analysis of rearrangements can effectively be done with the available software (e.g., the STAR-fusion https://github.com/STAR-Fusion/STAR-Fusion), analysis of somatic mutations and gene expression can be more challenging in such datasets.

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

This script can be used to check the selected primers' possibility to hybridize with ribosomal RNAs. The 



![image](https://github.com/MitiushkinaNV/RACE_NGS/assets/96590759/372639ee-f20a-48fa-bc00-1783ea5a2a9f)

