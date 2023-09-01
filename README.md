# RACE_NGS
This collection of python scripts was developed to help in the design of the targeted RNA panels on the basis of anhored multiplex PCR-based methods (primarily, 3' and 5' rapid anplification of cDNA ends or RACE), and for the following analysis of the obtained data. While the analysis of rearrangements can effectively be done with the available software (e.g., the STAR-fusion https://github.com/STAR-Fusion/STAR-Fusion), analysis of somatic mutations and gene expression can be more challenging in such datasets.

FindRACEprimers.py 

This is the main script used for the primer design for 3'RACE or 5'RACE-based approaches. The input coordinates file in .tsv format should contain at least 3 necessary fields: "Chr"	"Start"	"End" (1-based coordinates are required). It is recomended to provide the coordinates for the whole exons, if not too long, for the convenience of the following data analysis. For the 3'RACE select --primer_strand as "forward" and for the 5'RACE select --primer_strand as "reverse". Set the desired primers' melting temperature with --melt_temp parameter. The following arguments are required for the calculation of melting temperatures of potential primers with the SantaLucia method (SantaLucia, 1998), with the corrections, suggested by von Ahsen et al. (1999): 

  --conc_primer, Primer concentration in PCR, nM
  
  --conc_mg, Mg concentration in PCR
  
  --conc_salt, Salt concentration in PCR, mM
  
The code for the calculation of melting temperature is adapted from http://www.biophp.org/minitools/melting_temperature. 
