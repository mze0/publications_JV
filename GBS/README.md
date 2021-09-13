## step1 CNV part:
Rscript GBS_CNV_part.R Dir_with_bam_file Result_folder Bin_size 'Disease_loci'

## step2 Haplotyping part of converting sequencing to snparray format:
Rscript convertGBStoSNParray.R PGD_family_ID Phasing_seed Working_Dir

## step3 siCHILD part:
Rscript siCHILD_core.R PGD_family_ID SNParray_format_input_data parameter_file interval_file out_Dir script_Dir 
