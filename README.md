# PRIME
Prime - Profile Motif Search based on TSS


**PRIME Parameters:**
```
-f --fasta  
                  Genome file to analyze in FASTA format <FILE> (default="")
-p --pssm
                  PSSM (position-specific scoring matrix) file, otherwise it will be computed <FILE> [OPTIONAL, if parameter -m is used] (default="")
-m --motif
                  Hand selected motif patterns to build up a PSSM <FILE> [OPTIONAL, if parameter -p is used] (default="")
-s --tss
                  TSS (transcriptional start sites) file coming from ReadXplorer "https://www.uni-giessen.de/fbz/fb08/Inst/bioinformatik/software/ReadXplorer"
                  <FILE> (default="")
-x --pssm_s1      
                  PSSM start min.: minimum upstream distance <INT value> (default=6)
-y --pssm_s2
                  PSSM start max.: maximum upstream distance <INT value> (default=8)
-c --pseudo_count
                  Value for pseudocounts (PSSM calculation). Value is added to zero as well as to non-zero values <FLOAT value> (default=0.7)
-n --num_bkg_runs
                  Number of randomly picked PSSM scores from a given genome to set up a suitable background model <INT value> (default=100000)
-a --alpha        
                  Adjusting significance level alpha [0.0, 1.0] <FLOAT value> (default=0.05)
-o --out_folder
                  Folder is used to store the final results <FILE> (default="")
-d --description
                  Feature type; only used for the GFF file <STRING value> (default="TATA")
-k --distribution
                  PSSM calculated values are normal distributed or follow a generalized extreme value distribution [norm, genextreme] (default="norm")
```
