# CHARM-SexDiff-Splicing

A collection of jupyter notebooks for analyzing sex differences in blood RNA-seq alternative splicing in the CHARM cohort.

## Content

```{bash}
├── 00-preprocess.ipynb          # Preprocessing, clean-up meta data annotation, analysis for batch effects
├── 01-exclusion.ipynb           # Cohort definition; Exclusion of seropositive subjects from downstream analyses
├── 02_SE.ipynb                  # Visualization of Skipped Exons analysis
├── 03-A5SS.ipynb                # Alternative 5' Splice Site
├── 04-A3SS.ipynb                # Alternative 3' Splice Site
├── 05-RI.ipynb                  # Retained Introns
├── 06-ISG_compilation.ipynb     # Use projection matrix to compute ISG LVs for all samples
├── 07-mediation.ipynb           # Causal mediation analysis using R package (Figure 4)
├── 07-total-mediation.ipynb     # Total mediation effects calculated in two ways (Figure 4)
├── 08-summarize.ipynb           # Get summary numbers; Example plot (Figure 2)
├── 09-viral_load.ipynb          # Viral Load analysis (Figure 1)
├── 10-network.ipynb             # Prepare gene sets for Network analysis (Figure 2)
├── md5sum_matching.ipynb        # Check MD5sums to match fastq files between versions 
└── run_jem.py                   # Script for generating results in 02-05
```
