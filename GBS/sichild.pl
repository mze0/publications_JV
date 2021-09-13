package sichild;
use strict;
use warnings 'all';
use File::Basename;
use IO::File;


our $java_pref;
my $UsedHOST;
my %external_versions = (
    samtools    =>      '0.1.19',
    sichild	=>	'v_2018_02_28',	#'v_2016_04_15',
);


my $PGD_SERVER = "/uz/data/hydra/genomicscore/";
my $RSCRIPT = "/uz/data/hydra/shared_app/apps/R/3.4.0/bin/";


# Default settings, use key = value pair in sets_variants_freebayes.request file to modify
sub init{
    my %default_values = (
    species                 =>      'homo_sapiens',
    genome_build            =>      'hg19',
    mismatches		    =>      'false'     #true if a graph of the mismatches must be made
    );

    &::addIfMissing(\%default_values);

    my @required_in_request = (
    );

    my $failed  = &::failIfMissing(\@required_in_request);

    $failed;
}



sub sichild() {

    &::myprt(0,"Running sichild version $::versions{'sichild'}");
    my $failed = &init;
    my @shFileList = ();
    my $QUEUE = "lowprio";  #"lowmem"; 
#    my $JOBS_DIR = "$::cf{task_dir}/jobs"; &mysys(1,"mkdir -p $JOBS_DIR");

    my $RESULT_DIR = &::makeResultDir;
    my $RESULT_TMP_DIR = &::makeResultTmpDir;

#    my $RESULT_DIR = "$::cf{task_dir}/result";
#    my $RESULT_TMP_DIR = "$RESULT_DIR/tmp";  &mysys(1,"mkdir -p $RESULT_TMP_DIR");
    my $SICHILD_TO_GC_DIR = "$::cf{script_dir}/sichild/$::versions{'sichild'}";

    my $PGD_EXPORTED_DIR = "$PGD_SERVER/cme_genome_raw/PGD_Exported/$::cf{nameSampleOrSet}";
    my $PGD_SICHILD_DIR = "$PGD_SERVER/cme_genome_calculated/siCHILD";
    my $PGD_PARAMETERS_FIXED = "$SICHILD_TO_GC_DIR/Fixed_Parameters.tsv";
    my $PGD_DATA = "$PGD_EXPORTED_DIR/$::cf{nameSampleOrSet}.txt";

    my $PGD_DATA_ADJ = "$RESULT_TMP_DIR/$::cf{nameSampleOrSet}.adj";
    my $PGD_PARAMETERS = "$PGD_EXPORTED_DIR/$::cf{nameSampleOrSet}_Parameters.txt";
    my $PGD_PARAMETERS_ADJ = "$RESULT_TMP_DIR/$::cf{nameSampleOrSet}_Parameters.adj";
    my $PGD_INTERVALS = "$PGD_EXPORTED_DIR/$::cf{nameSampleOrSet}_Intervals.txt";
    my $PGD_DATA_INTERVALS = "$RESULT_TMP_DIR/$::cf{nameSampleOrSet}_Intervals.txt";
    my $PGD_RESULT_DIR = "$PGD_SICHILD_DIR/$::cf{nameSampleOrSet}/current";
#    my $failed = 0;

    &mymail(1,"jia.ding\@uzleuven.be,cindy.melotte\@uzleuven.be,eftychia.dimitriadou\@uzleuven.be,luc.dehaspe\@uzleuven.be","genomicscore\@uzleuven.be","Started task siCHILD $::cf{nameSampleOrSet}","in: $PGD_EXPORTED_DIR\nprogress: $::cf{task_dir}/$::cf{task}.log");

    if (-e $PGD_DATA) {
	&myprt(1,"Found data file $PGD_DATA");
    } else {
	&myprt(1,"Could not find data file $PGD_DATA");
	$failed = 1;
    }
    if (-e $PGD_PARAMETERS) {
	&myprt(1,"Found parameters file $PGD_PARAMETERS");
    } else {
	&myprt(1,"Could not find parameters file $PGD_PARAMETERS");
	$failed = 1;
    }
    if (-e $PGD_INTERVALS) {
	&myprt(1,"Found intervals file $PGD_INTERVALS");
    } else {
	&myprt(1,"Could not find intervals file $PGD_INTERVALS");
	$failed = 1;
    }

    unless ($failed) {
    	&mysys(1,"sed -e 's/,/\\./g' $PGD_DATA > $PGD_DATA_ADJ");
	&mysys(1,"cp $PGD_PARAMETERS_FIXED  $PGD_PARAMETERS_ADJ");
	&mysys(1,"cp $PGD_INTERVALS $PGD_DATA_INTERVALS");	
	&mysys(1,"echo 'GC_File\t$SICHILD_TO_GC_DIR/SeqStat_CytoSNP12_Window10000.txt' >> $PGD_PARAMETERS_ADJ");
	&mysys(1,"echo 'siCHILD_DIR\t$SICHILD_TO_GC_DIR' >> $PGD_PARAMETERS_ADJ");
	&mysys(1,"sed -e 's/ *= */\\t/' $PGD_PARAMETERS | awk -F'\t' 'BEGIN{OFS=FS}{for (i=1; i<=NF; ++i) {\$i=toupper(substr(\$i,1,1)) tolower(substr(\$i,2))} print }' >> $PGD_PARAMETERS_ADJ");
	&mysys(1,"$RSCRIPT $SICHILD_TO_GC_DIR/siCHILD_core.R $::cf{nameSampleOrSet} RESULT_TMP_DIR $RESULT_DIR"); 	

    	$failed = &execAndWait(1,"siC.plot.$::cf{nameSampleOrSet}",\@shFileList,$JOBS_DIR,"lowmem");
    }

    unless($failed) {
	&mysys(1,"mkdir -p $PGD_RESULT_DIR");
	&mysys(1,"chmod 500 -R $RESULT_DIR/*; cp -r $RESULT_DIR/* $PGD_RESULT_DIR/.; chmod 555 -R $RESULT_DIR/*");
	&mysys(1,"chmod 500 -R $RESULT_DIR/*; cp -r $RESULT_DIR/* $PGD_RESULT_DIR/.; chmod 555 -R $RESULT_DIR/*");
	&mymail(1,"jia.ding\@uzleuven.be,cindy.melotte\@uzleuven.be,eftychia.dimitriadou\@uzleuven.be,luc.dehaspe\@uzleuven.be","genomicscore\@uzleuven.be","siCHILD $::cf{nameSampleOrSet} finished","in: $PGD_EXPORTED_DIR\nout: $PGD_RESULT_DIR");
    }

    &grantPermissions($::cf{task_dir});

    $failed; 
}

1;
