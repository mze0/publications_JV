## CNV part:
Rscript GBS_CNV_part.R Dir_with_bam_file Result_folder Bin_size 'Disease_loci'

## Haplotyping part of converting sequencing to snparray format:
Rscript convertGBStoSNParray.R PGD_family_ID Phasing_seed Working_Dir