#!/usr/bin/env perl

if(@ARGV != 3){
    print "%perl " . $0 . " [scaleFactor.txt] [Absolute input_MACS2_directory] [Absolute output_directory]\n";
    exit;
}

### Configuration
require("./ATACprofWS.conf");

$scale_file = shift(@ARGV);
$input_dir  = shift(@ARGV);
$output_dir = shift(@ARGV);

if(!-e $scale_file){
    print "scaleFactor file " . $scale_file . " was not found\n";
    print "Aborted\n";
    exit;
}
if(substr($input_dir, -1) ne "/"){
    $input_dir .= "/";
}
if(substr($output_dir, -1) ne "/"){
    $output_dir .= "/";
}

print "## Reading scale factors from " . $scale_file . "\n";
undef(%SCALES);
undef(@REP_IDs);
undef(%SAMPLE_IDs);

open(IN, $scale_file);
while($line = <IN>){
    if($line =~ /^\#/){ next; }
    if($line =~ /_Norm/){ next; }

    chop($line);
    undef(@data);
    @data = split(/\t/, $line);

    # F4
    $SCALES{ $data[0] } = $data[ scalar(@data) -1];
    push(@REP_IDs, $data[0]);
    
    # sample ID is like "replicate_number"
    $data[0] =~ /(.*)\_(\d+)$/;
    $SID = $1;
    $RID = $2;
    $SAMPLE_IDs{ $SID }{$RID} = 1;
}
close(IN);

$howmany = keys %SCALES;
print "#\t" . $howmany . " SFs\n";
if($howmany eq "" || $howmany == 0){
    print "Error in reading scale factors\n";
    print "Aborted\n";
    exit;
}

$PV_PROG = $ATACprofWS_ROOT . "script/add_Pvalue.pl";
if(!-e $PV_PROG){
    print $PV_PROG . " was not found\n";
    print "Aborted" . "\n";
    exit;
}


$howmany_samples = keys %SAMPLE_IDs;
print "#\t" . $howmany_samples . " SFs: ";
print join(" ", keys %SAMPLE_IDs) . "\n";
print "\n";

# This is for common peaks in KO and Ctrl samples. FDR filtering is not applied.

$MACS2    = "postMACS2_v4-1";

## Control 
$this_nc_peak = $input_dir . $CTRL . "/". $CTRL . "_final.narrowPeak.bed.gz";
if(!-e $this_nc_peak){
    print "Not NC exist " . $this_nc_peak . "\n";
    exit;
}
## This file: chr start end ID rep1 rep2 ....
## Here, we assumed the number of replicates is two.
## If you have more two replicates,, you have to modify this part

$command1 = $command2 = "";
$this_SID = $CTRL;
$this_SF1 = $SCALES{ $this_SID . "_1" };
$this_SF2 = $SCALES{ $this_SID . "_2" };

$this_nc_peak_cut = $input_dir . $this_SID . "/". $this_SID . "_final.narrowPeak.cut.bed.gz";

$REP1_nc_target_bed_sw = $input_dir . $this_SID . "/" . $this_SID . "_1_150.bed.gz";
$REP2_nc_target_bed_sw = $input_dir . $this_SID . "/" . $this_SID . "_2_150.bed.gz";

if(!-e $REP1_nc_target_bed_sw){
    print $REP1_nc_target_bed_sw . " was not found\n";
    exit;
}
if(!-e $REP2_nc_target_bed_sw){
    print $REP2_nc_target_bed_sw . " was not found\n";
    exit;
}


if($this_SF1 eq "" || $this_SF2 eq "" ){
    print "## Scale factor for " . $this_SID . " was empty\n";
    print "Aborted" . "\n";
    exit;
}

print "\n";
print "#" . $this_SID . " is scaled and averaged and filtered" . "\n";

# Remove ""SCALED"" peaks <= cut
$this_nc_peak_cut = $input_dir . $CTRL . "/". $CTRL . "_final.narrowPeak.cut.bed.gz";
unlink($this_nc_peak_cut);
#print "## Create " . $this_nc_peak_cut . "\n";

$this_nc_peak_cut_tmp = $input_dir . $CTRL . "/". $CTRL . "_final.narrowPeak.cut.bed";
open(OUT, ">$this_nc_peak_cut_tmp");
print OUT "# Scale factor for Rep1 : " . $this_SF1 . "\n";
print OUT "# Scale factor for Rep2 : " . $this_SF2 . "\n";
print OUT "# Cut if 'scaled'_avg <= "    . $MIN_SCALED_READ_CUT . "\n";
print OUT "#Chr Start  End ID #.reads_Rep1 #.reads_Rep2 Scaled_Rep1 Scaled_Rep2 Scaled_Avg" . "\n";
close(OUT);

#chr start end ID raw_rep1 raw_rep2 scaled_rep1 scaled_rep2 scaled_ave
$command1 = " zcat " . $this_nc_peak . " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$5 * $this_SF1; s2=\$6 * $this_SF2; ss = (s1+s2)/2; if ( ss > $MIN_SCALED_READ_CUT ){ print \$0 \"\\t\" s1 \"\\t\" s2 \"\\t\" ss}}\' >> " . $this_nc_peak_cut_tmp;
$command2 = "gzip -f " . $this_nc_peak_cut_tmp;
system($command1);
system($command2);
#===========================================#

if(!-e $this_nc_peak_cut || -z $this_nc_peak_cut){
    print "Error in generating " . $this_nc_peak_cut . "\n";
    print "Aborted." . "\n";
    exit;
}
print "#" . $this_nc_peak_cut . "\n";
$REP1_SF_nc = $this_SF1;
$REP2_SF_nc = $this_SF2;
#========= End of Ctrl   


# Samples
$command3 = $command4 = "";
foreach $key (keys %SAMPLE_IDs){
    $this_SID = $key;
    if($this_SID eq $CTRL){ next; }

    $REP1_SF = $SCALES{ $this_SID . "_1" };
    $REP2_SF = $SCALES{ $this_SID . "_2" };
    $this_peak = $input_dir . $this_SID . "/". $this_SID . "_final.narrowPeak.bed.gz";
    
    if($REP1_SF eq "" || $REP2_SF eq "" ){
	print "## Scale factor for " . $this_SID . " was empty\n";
	print "Aborted" . "\n";
	exit;
    }

    $REP1_target_bed_sw = $input_dir . $this_SID . "/" . $this_SID . "_1_150.bed.gz";
    $REP2_target_bed_sw = $input_dir . $this_SID . "/" . $this_SID . "_2_150.bed.gz";
    
    if(!-e $REP1_target_bed_sw){
	print $REP1_target_bed_sw . " was not found\n";
	exit;
    }
    if(!-e $REP2_target_bed_sw){
	print $REP2_target_bed_sw . " was not found\n";
	exit;
    }

    print "#" . $this_SID . " is scaled and averaged and filtered" . "\n";
    $this_output = $output_dir . $this_SID . "/";
    print $this_output . "\n";

    if(!-e $this_output){
	system("mkdir -p $this_output");
    }

    # Remove ""SCALED"" peaks <= cut
    $this_peak_cut = $this_output . $this_SID . "_final.narrowPeak.cut.bed.gz";
    unlink($this_peak_cut);
    $this_peak_cut_tmp  = $this_output . $this_SID . "_final.narrowPeak.cut.bed";
    $this_peak_cut_tmp2 = $this_output . $this_SID . "_final.narrowPeak.cut.tmp.bed";


    # STEP 1: filtering PEAKS using $MIN_SCALED_READ_CUT (nc sample was done already)
    $command3  = "if [ -e " . $this_peak_cut_tmp2 . " ]; then " . "\n";
    $command3 .= " rm " . $this_peak_cut_tmp2 . "\n";
    $command3 .= "fi" . "\n";

    $command3 .= "if [ -e " . $this_peak_cut_tmp . " ]; then " . "\n";
    $command3 .= " rm " . $this_peak_cut_tmp . "\n";
    $command3 .= "fi" . "\n";
    
    $command3 .= "echo \"# Sample: $this_SID \" > " . $this_peak_cut_tmp . "\n";
    $command3 .= "echo \"# Ctrl: $CTRL \" >> " . $this_peak_cut_tmp . "\n";
    $command3 .= "echo \"# Scale factor for Rep1 : $REP1_SF \" >> " . $this_peak_cut_tmp . "\n";
    $command3 .= "echo \"# Scale factor for Rep2 : $REP2_SF \" >> " . $this_peak_cut_tmp . "\n";
    $command3 .= "echo \"# Scale factor for ctrl_Rep1 : $REP1_SF_nc \" >> " . $this_peak_cut_tmp . "\n";
    $command3 .= "echo \"# Scale factor for ctrl_Rep2 : $REP2_SF_nc \" >> " . $this_peak_cut_tmp . "\n";
    $command3 .= "echo \"# Cut if scaled_avg <=  $MIN_SCALED_READ_CUT \" >> " .  $this_peak_cut_tmp . "\n";
    $command3 .= "echo \"#Chr Start  End ID #.reads_Rep1 #.reads_Rep2 Scaled_Rep1 Scaled_Rep2 Scaled_Avg Ctrl_Rep1 Ctrl_Rep2 ctrl_sR1 ctrl_sR2 Ctrl_sR_Avg\" >> " . $this_peak_cut_tmp . "\n\n";

    $command3 .= "zcat " . $this_peak . " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$5 * $REP1_SF; s2=\$6 * $REP2_SF; ss = (s1+s2)/2; if ( ss > $MIN_SCALED_READ_CUT ){ print \$0 \"\\t\" s1 \"\\t\" s2 \"\\t\" ss}}\' > " . $this_peak_cut_tmp2;

    ## ADD scaled reads of counter part (=CTRL)
    $command3 .= "\n" . $COVB . " -counts -a " . $this_peak_cut_tmp2 . " -b " . $REP1_nc_target_bed_sw;
    $command3 .= " | " . $COVB . " -counts -a stdin -b " . $REP2_nc_target_bed_sw;

    $command3 .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$10 * $REP1_SF_nc; s2=\$11 * $REP2_SF_nc; ss = (s1+s2)/2;  print \$0 \"\\t\" s1 \"\\t\" s2 \"\\t\" ss}\' >> " . $this_peak_cut_tmp . "\n\n";

    $command4  = "if [ -e " . $this_peak_cut . " ]; then " . "\n";
    $command4 .= " rm " . $this_peak_cut . "\n";
    $command4 .= "fi" . "\n";
    $command4 .= "gzip " . $this_peak_cut_tmp . "\n";
    $command4 .= "rm " . $this_peak_cut_tmp2 . "\n\n";;


    ##$this_nc_peak_cut = $temp . "/MACS2/NC/NC_final.narrowPeak.cut.bed.gz";
    $this_nc_peak_cut2 = $this_output . $this_SID . "_nc_final.narrowPeak.cut.bed";
    $command4 .= "\n\n" . "if [ -e " . $this_nc_peak_cut2 . " ]; then " . "\n";
    $command4 .= " rm " . $this_nc_peak_cut2 . "\n";
    $command4 .= "fi" . "\n";

    $command4 .= "echo \"# Sample: $CTRL \" > "    . $this_nc_peak_cut2 . "\n";
    $command4 .= "echo \"# Ctrl: $this_SID \" >> " . $this_nc_peak_cut2 . "\n";

    $command4 .= "echo \"# Scale factor for Rep1 : $REP1_SF_nc \" >> " . $this_nc_peak_cut2 . "\n";
    $command4 .= "echo \"# Scale factor for Rep2 : $REP2_SF_nc \" >> " . $this_nc_peak_cut2 . "\n";
    $command4 .= "echo \"# Scale factor for ctrl_Rep1 : $REP1_SF \" >> " . $this_nc_peak_cut2 . "\n";
    $command4 .= "echo \"# Scale factor for ctrl_Rep2 : $REP2_SF \" >> " . $this_nc_peak_cut2 . "\n";
    $command4 .= "echo \"# Cut if scaled_avg <=  $MIN_SCALED_READ_CUT \" >> " .  $this_nc_peak_cut2 . "\n";
    $command4 .= "echo \"#Chr Start  End ID #.reads_Rep1 #.reads_Rep2 Scaled_Rep1 Scaled_Rep2 Scaled_Avg Ctrl_Rep1 Ctrl_Rep2 ctrl_sR1 ctrl_sR2 ctrl_sAvg\" >> " . $this_nc_peak_cut2 . "\n\n";
		
    $command4 .= "\n" . $COVB . " -counts -a " . $this_nc_peak_cut . " -b " . $REP1_target_bed_sw;
    $command4 .= " | " . $COVB . " -counts -a stdin -b " . $REP2_target_bed_sw;
    $command4 .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$10 * $REP1_SF; s2=\$11 * $REP2_SF; ss = (s1+s2)/2;  print \$0 \"\\t\" s1 \"\\t\" s2 \"\\t\" ss}\' | sort -V -k1,1 -k 2,2 >> " . $this_nc_peak_cut2 . "\n\n";


    # STEP 2: genome wide Fold change 10K-bp (howmany: 308,839) using cut_PEAKs
    $this_peak_cut_10K = $this_output . $this_SID . "_final.narrowPeak.cut.10K.bed.gz";
    
    $outbed = $this_peak_cut_10K;
    $winbed = $GL_10K;

    $commond_line  = "";
    &GET_Command_Win();
    $command4_1 = $command_line;


    #------------------ All Peaks in KO
    # 1     2      3    4    5                   6                     7                  8                 9
    # chr start end ID #.reads_Rep1 #.reads_Rep2 scaled_Rep1 scaled_Rep2 scaled_avg
    # CUP: common and unique peaks
    # Chr Start  End ID r1 r2 sr1 sr2 savg #.reads_Rep1 #.reads_Rep2 Scaled_Rep1 Scaled_Rep2 Scaled_Avg Ctrl_Rep1 Ctrl_Rep2 ctrl_sR1 ctrl_sR2 Ctrl_sR_Avg
    
    # STEP 1: for calc. FC, PV, FDR with NC reads (not considering NC peaks)
    #Chr Start  End ID #.reads_Rep1 #.reads_Rep2 Scaled_Rep1 Scaled_Rep2 Scaled_Avg
    # skip 5,6
    $this_common     =  $this_output . $this_SID . "_common.bed";
    $this_common_tmp =  $this_output . $this_SID . "_common.tmp.bed";
    $command5  = "if [ -e " . $this_common_tmp . " ]; then " . "\n";
    $command5 .= " rm " . $this_common_tmp . "\n";
    $command5 .= "fi" . "\n";
    $command5 .= "\n";
    # chr start end ID r1 r2 s1 s2 avg r1 r2 c1 c2 avg (c=NC)
    # chr start end ID r1 r2 s1 s2 avg r1 r2 c1 c2 avg (c=KO)
    $command5 .= $INTB . " -wo -a " . $this_peak_cut . " -b " . $this_nc_peak_cut2;
    $command5 .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$3-\$2; s2=\$17-\$16; if ((\$29/s1 >= 0.5) || (\$29/s2 >= 0.5)){ print \$0}}\' | cut -f 1-28 | sort -V -k1,1 -k2,2 | uniq > " . $this_common_tmp;

    ## one KO peak can include multiple NC peaks
    ## So, get only the list of KO peaks, and calc. NC scaled reads in these (??? NO). Take average!!
    $command5 .= "\n\n" . "if [ -e " . $this_common . " ]; then " . "\n";
    $command5 .= " rm " . $this_common . "\n";
    $command5 .= "fi" . "\n";
    $command5 .= "\n";

    #Chr Start  End ID Scaled_Rep1 Scaled_Rep2 Scaled_Avg Scaled_NC1 Scaled_NC2 Scaled_Avg_NC
    $command5 .= $MERGEB . " -i " . $this_common_tmp . " -c 4,7,8,9,21,22,23 -o distinct,sum,sum,sum,mean,mean,mean -d -1 | sort -V -k1,1 -k2,2 > " . $this_common . "\n";
    $command5 .= "\n" . "rm " . $this_common_tmp;
    $command5 .= "\n";
    #------------------ End of All Peaks in KO

    #STEP 2: add NC unique info. to this_CUP
    $this_nc_common_id   =  $this_output . $this_SID . "_nc_common.id.bed";
    $this_nc_denovo_id   =  $this_output . $this_SID . "_nc_uniq.id.bed";
    $this_nc_denovo      =  $this_output . $this_SID . "_nc_uniq.bed";

    # NC unique
    # output NC peak list (split NC parts)
    $command6_1  = "if [ -e " . $this_nc_common_id . " ]; then " . "\n";
    $command6_1 .= " rm " . $this_nc_common_id . "\n";
    $command6_1 .= "fi"   . "\n";
    # chr start end ID r1 r2 s1 s2 avg r1 r2 c1 c2 avg (c=NC)14
    # chr start end ID r1 r2 s1 s2 avg r1 r2 c1 c2 avg (c=KO)14
    $command6_1 .= $INTB . " -wo -a " . $this_peak_cut . " -b " . $this_nc_peak_cut2;
    $command6_1 .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$3-\$2; s2=\$17-\$16; if ((\$29/s1 >= 0.5) || (\$29/s2 >= 0.5)){ print \$0}}\' | cut -f 15-18 | sort -V -k1,1 -k2,2 | uniq > " . $this_nc_common_id;
    # this is common
    
    # and now NC denovo with reads info.
    $command6_2  = "if [ -e " . $this_nc_denovo . " ]; then " . "\n";
    $command6_2 .= " rm " . $this_nc_denovo . "\n";
    $command6_2 .= "fi" . "\n";
    $command6_2 .= $INTB . " -v -a " . $this_nc_peak_cut2 . " -b " . $this_nc_common_id;
    # chr start end ID r1 r2 s1 s2 avg r1 r2 c1 c2 avg (c=KO)14
    # whether denovo set to zero or take counter-part reads?
    # Flip KO -> NC
    
    # this one is taking KO mapped reads
    $command6_2 .= " | cut -f 1-4,7-9,12-14 | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ print \$1 \"\\t\" \$2 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$8 \"\\t\" \$9 \"\\t\" \$10 \"\\t\" \$5 \"\\t\" \$6 \"\\t\" \$7 }\' | sort -V -k1,1 -k2,2 | uniq > " . $this_nc_denovo;
    
    # this one is set to zero
    #$command6_2 .= " | cut -f 1-4,7-9,12-14 | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ \$8=0; \$9=0; \$10=0; print \$1 \"\\t\" \$2 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$8 \"\\t\" \$9 \"\\t\" \$10 \"\\t\" \$5 \"\\t\" \$6 \"\\t\" \$7 }\' | sort -V -k1,1 -k2,2 | uniq > " . $this_nc_denovo;
    
    $command6_2 .= "\n\n" . "cut -f 1-4 " . $this_nc_denovo . " > " . $this_nc_denovo_id;
    $command6_2 .= "\n\n" . "rm " . $this_nc_common_id;
    
    #STEP 3: add KO unique info. to this_CUP
    $this_common_id     = $this_output . $this_SID . "_common.id.bed";
    $this_ko_denovo_id  = $this_output . $this_SID . "_ko_uniq.id.bed";
    $this_ko_denovo     = $this_output . $this_SID . "_ko_uniq.bed";
    
    $command7  = "if [ -e " . $this_common_id . " ]; then " . "\n";
    $command7 .= " rm " . $this_common_id . "\n";
    $command7 .= "fi" . "\n";
    $command7 .= $INTB . " -wo -a " . $this_peak_cut . " -b " . $this_nc_peak_cut2;
    $command7 .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$3-\$2; s2=\$17-\$16; if ((\$29/s1 >= 0.5) || (\$29/s2 >= 0.5)){ print \$0}}\' | cut -f 1-4 | sort -V -k1,1 -k2,2 | uniq > " . $this_common_id;
    
    # chr start end ID r1 r2 s1 s2 avg r1 r2 c1 c2 avg (c=NC)14
    $command7 .= "\n" . $INTB . " -v -a " . $this_peak_cut . " -b " . $this_common_id;
    
    # take counter-part mapped reads
    $command7 .= " | cut -f 1-4,7-9,12-14 | sort -V -k1,1 -k2,2 | uniq > " . $this_ko_denovo;
    
    # set to all zero
    #$command7 .= " | cut -f 1-4,7-9,12-14 | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ \$8=0; \$9=0; \$10=0; print \$1 \"\\t\" \$2 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$5 \"\\t\" \$6 \"\\t\" \$7 \"\\t\" \$8 \"\\t\" \$9 \"\\t\" \$10 }\' | sort -V -k1,1 -k2,2 | uniq > " . $this_ko_denovo;
    
    $command7 .= "\n\n";
    $command7 .= "cut -f 1-4 " . $this_ko_denovo . " > " . $this_ko_denovo_id;
    $command7 .= "\n\n";
		
    ## Calc. PV and FDR
    # now, this_common this_ko_denovo this_nc_denovo
		
    #All in one -> Calc. FDR -> split
    $this_CUP     = $this_output . $this_SID . "_CUP.bed";
    $this_CUP_fdr = $this_output . $this_SID . "_CUP.fdr.bed";

    #Chr Start  End ID Scaled_Rep1 Scaled_Rep2 Scaled_Avg Scaled_NC1 Scaled_NC2 Scaled_Avg_NC logFC PV FDR
    $command7 .= "\n\n" . "cat " . $this_common . " " . $this_ko_denovo . " " . $this_nc_denovo . " | sort -V -k1,1 -k2,2 > " . $this_CUP;

    #-------------------------------------------------------
    #Split
    $this_common_fdr     = $this_output . $this_SID . "_common.fdr.bed";
    $this_ko_denovo_fdr  = $this_output . $this_SID . "_ko_uniq.fdr.bed";
    $this_nc_denovo_fdr  = $this_output . $this_SID . "_nc_uniq.fdr.bed";
    #&Split_Each();

    # norm not used
    $command7 .= "\n\n" . "perl " . $PV_PROG . " " . $this_common . " " . $this_common_fdr . " \"5,6\" \"8,9\" 0";

    # norm used at EdgeR
    $command7 .= "\n\n" . "perl " . $PV_PROG . " " . $this_ko_denovo . " " . $this_ko_denovo_fdr . " \"5,6\" \"8,9\" 1";
    $command7 .= "\n\n" . "perl " . $PV_PROG . " " . $this_nc_denovo . " " . $this_nc_denovo_fdr . " \"5,6\" \"8,9\" 1";

    # all in one
    $command7 .= "\n\n" . "cat " . $this_common_fdr . " " . $this_ko_denovo_fdr . " " . $this_nc_denovo_fdr . " | sort -V -k1,1 -k2,2 > " . $this_CUP_fdr;
    $command7 .= "\n\n";

    # show numbers
    $command10_1 = "echo \"ALL KO peaks (common+unique): \" `zcat $this_peak_cut | grep -v \"#\" | wc -l`";

    $command10_2   = "echo \"Common peaks: \" `cat $this_common_fdr | grep -v \"#\" | wc -l`";
    $command10_2_1 = "echo \"Common logFC>1 & FDR<0.05: \" `grep -v \"\#\" $this_common_fdr | awk '{if( \$11 > 1 && \$13 < 0.05){ print \$0 } }' | wc -l`";
    $command10_2_2 = "echo \"Common logFC< -1 & FDR<0.05: \" `grep -v \"\#\" $this_common_fdr | awk '{if( \$11 < -1 && \$13 < 0.05){ print \$0 } }' | wc -l`";

    $command10_3   = "echo \"KO unique peaks: \" `cat $this_ko_denovo_fdr | grep -v \"#\" | wc -l`";
    $command10_3_1 = "echo \"KO logFC>1 & FDR<0.05: \" `grep -v \"\#\" $this_ko_denovo_fdr | awk '{if( \$11 > 1 && \$13 < 0.05){ print \$0 } }' | wc -l`";
    $command10_3_2 = "echo \"KO logFC< -1 & FDR<0.05: \" `grep -v \"\#\" $this_ko_denovo_fdr | awk '{if( \$11 < -1 && \$13 < 0.05){ print \$0 } }' | wc -l`";

    $command10_4 = "echo \"NC unique peaks: \" `cat $this_nc_denovo_fdr | grep -v \"#\" | wc -l`";
    $command10_4_1 = "echo \"NC logFC>1 & FDR<0.05: \" `grep -v \"\#\" $this_nc_denovo_fdr | awk '{if( \$11 > 1 && \$13 < 0.05){ print \$0 } }' | wc -l`";
    $command10_4_2 = "echo \"NC logFC< -1 & FDR<0.05: \" `grep -v \"\#\" $this_nc_denovo_fdr | awk '{if( \$11 < -1 && \$13 < 0.05){ print \$0 } }' | wc -l`";

    
    $AGE = $this_output . $this_SID . "_pstMACS2.age";
    $MEM = "32G";
    $JID = $this_SID . "_pstMACS2";
    
    &THIS_QSUB();
    #system("qsub -v PATH=\$PATH $AGE");
    print "###########\n";
    print "%>sh " . $AGE . "\n";
    print "or\n";
    print "%>qsub -v $AGE" . "\n";
    print "###########\n";
    print "\n\n";

}

exit;



sub THIS_QSUB(){

    open(OUT, ">$AGE");
print OUT <<EOF;
#!/bin/bash
#\$ -S /bin/bash
#\$ -o $this_output
#\$ -e $this_output
#\$ -N $JID
#\$ -v LD_LIBRARY_PATH=""
#\$ -cwd -j y -l s_vmem=$MEM -l mem_req=$MEM

####cd $this_output

### pre-processing
$command3

$command4

### window blocks
$command4_1

$command4_2

$command4_3


### peaks (common, KO unique, NC uniue)
$command5

$command6_1

$command6_2

$command6_3

$command7

$command8

$command9


## Numbers
$command10_1

$command10_2
$command10_2_1
$command10_2_2

$command10_3
$command10_3_1
$command10_3_2

$command10_4
$command10_4_1
$command10_4_2

echo "DONE"
EOF
    close(OUT);

}


#$this_peak_cut, $this_nc_peak_cut2 are used
sub GET_Command_Win{
    
    # chr star end ID chr start end ID raw_rep1 raw_rep2 scaled_rep1 scaled_rep2 scaled_ave cov_bp
    $command_line  = "if [ -e " . $outbed . " ]; then " . "\n";
    $command_line .= " rm " . $outbed . "\n";
    $command_line .= "fi" . "\n";
    $command_line .= $INTB . " -wao -a " . $winbed . " -b " . $this_peak_cut;
    
    # 50% overlap (modified 2019/03)
    # chr start end ID r1 r2 s1 s2 avg r1 r2 c1 c2 avg
    $command_line .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ if(\$19>0){ s1=\$3-\$2; s2=\$7-\$6; if ((\$19/s1 >= 0.5) || (\$19/s2 >= 0.5)){ print \$0}else{ \$11=0; \$12=0; \$13=0; print \$0;} }else{ print \$0} }\'";
    
    # output: chr star end ID s1 s1 savg c1 c2 cavg
    $command_line .= " | " . $MERGEB . " -i stdin -c 4,11,12,13,16,17,18 -o distinct,sum,sum,sum,sum,sum,sum -d -1";
    # + chr srt end ID r1 r2 s1 s2 savg r1 r2 sample_c1 sample_c2 sample_avg
    $command_line .= " | " . $INTB . " -wao -a stdin " . " -b " . $this_nc_peak_cut2;
    $command_line .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ if(\$25>0){ s1=\$3-\$2; s2=\$13-\$12; if ((\$25/s1 >= 0.5) || (\$25/s2 >= 0.5)){ print \$0}else{ \$17=0; \$18=0; \$19=0; print \$0;} }else{ print \$0} }\'";
    
    # chr srt end ID s1 s2 savg c1 c2 avgc nc1 nc2 avgnc s1 s2 avgs
    $command_line .= " | " . $MERGEB . " -i stdin -c 4,5,6,7,8,9,10,17,18,19,22,23,24 -o distinct,max,max,max,max,max,max,sum,sum,sum,sum,sum,sum -d -1 | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ if(\$5 <1){ \$5=0 };if(\$6 <1){ \$6=0};if(\$7 <1){ \$7=0};if(\$8 <1){ \$8=0};if(\$9 <1){ \$9=0};if(\$10 <1){ \$10=0}; print \$0 }\' > " . $outbed . ".tmp";
    
    $command_line .=  "\n\n". "awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ s1=\$5+\$6; s2=\$11+\$12; if(s1>0 && s2>0){ print \$0 } }\' " . $outbed . ".tmp" . " | cut -f 1-4,5-7,11-13 > " . $outbed . ".common.bed";

    #get to zeros!!!
    #$command_line .= "\n\n" . "awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ s1=\$5+\$6; s2=\$11+\$12; if(s1>0 && s2==0){ print \$0 } }\' " . $outbed . ".tmp" . " | cut -f 1-4,5-7,11-13 > " . $outbed . ".ko_denovo.bed";

    #get counter-part mapped reads
    $command_line .= "\n\n" . "awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ s1=\$5+\$6; s2=\$11+\$12; if(s1>0 && s2==0){ print \$0 } }\' " . $outbed . ".tmp" . " | cut -f 1-4,5-7,8-10 > " . $outbed . ".ko_denovo.bed";

    #get counter-part mapped reads and flip ko and nc
    $command_line .= "\n\n" . "awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ s1=\$5+\$6; s2=\$11+\$12; if(s1==0 && s2>0){ print \$0 } }\' " . $outbed . ".tmp" . " | cut -f 1-4,11-13,14-16 | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{ print \$1 \"\\t\" \$2 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$8 \"\\t\" \$9 \"\\t\" \$10 \"\\t\" \$5 \"\\t\" \$6 \"\\t\" \$7 }\' > " . $outbed . ".nc_denovo.bed";

    # all in one
    #$command_line .= "\n\n" . "cat " . $outbed . ".common.bed" . " " . $outbed . ".ko_denovo.bed" . " " . $outbed . ".nc_denovo.bed | sort -V -k1,1 -k2,2 | uniq | gzip -cf > " . $outbed . "\n\n";

    # split
    ## Calc. PV and FDR
    $command_line .="\n\n". "perl " . $PV_PROG . " " . $outbed.".common.bed" . " " . $outbed.".common.fdr.bed" . " \"5,6\" \"8,9\" 0";

    ## Calc. PV and FDR
    $command_line .="\n\n". "perl " . $PV_PROG . " " . $outbed.".ko_denovo.bed" . " " . $outbed.".ko_denovo.fdr.bed" . " \"5,6\" \"8,9\" 1";

    ## Calc. PV and FDR
    $command_line .="\n\n". "perl " . $PV_PROG . " " . $outbed.".nc_denovo.bed" . " " . $outbed.".nc_denovo.fdr.bed" . " \"5,6\" \"8,9\" 1";

    ## Merge
    $outbed_fdr = $outbed;
    $outbed_fdr =~ s/\.bed.gz$//;
    $outbed_fdr .= ".fdr.bed";
    $command_line .= "\n\n" . "cat " . $outbed . ".common.fdr.bed" . " " . $outbed . ".ko_denovo.fdr.bed" . " " . $outbed . ".nc_denovo.fdr.bed | sort -V -k1,1 -k2,2 | uniq > " . $outbed_fdr . "\n\n";

    $command_line .= "\n\n" . "cat " . $outbed . ".common.bed" . " " . $outbed . ".ko_denovo.bed" . " " . $outbed . ".nc_denovo.bed | sort -V -k1,1 -k2,2 | uniq | gzip -cf > " . $outbed . "\n\n";

    $command_line .= "\n" . "rm ". $outbed . ".tmp";
    
}


sub Split_Each{

    $command8  = "if [ -e " . $this_ko_denovo_fdr . " ]; then " . "\n";
    $command8 .= " rm " . $this_ko_denovo_fdr . "\n";
    $command8 .= "fi" . "\n";
    $command8 .= $INTB . " -f 1 -a " . $this_CUP_fdr . " -b " . $this_ko_denovo_id . " | sort -V -k1,1 -k2,2 | uniq > " . $this_ko_denovo_fdr;
    $command8 .= "\n" . "rm " . $this_ko_denovo_id;
    
    $command8 .= "\n\n" . "if [ -e " . $this_common_fdr . " ]; then " . "\n";
    $command8 .= " rm " . $this_common_fdr . "\n";
    $command8 .= "fi" . "\n";
    $command8 .= "\n" . $INTB . " -f 1 -a " . $this_CUP_fdr  . " -b " . $this_common_id . " | sort -V -k1,1 -k2,2 | uniq > " . $this_common_fdr;
    $command8 .= "\n" . "rm " . $this_common_id;
    
    $command9  = "if [ -e " . $this_nc_denovo_fdr . " ]; then " . "\n";
    $command9 .= " rm " . $this_nc_denovo_fdr . "\n";
    $command9 .= "fi" . "\n";
    $command9 .= $INTB . " -f 1 -a " . $this_CUP_fdr  . " -b " . $this_nc_denovo_id . " | sort -V -k1,1 -k2,2 | uniq > " . $this_nc_denovo_fdr;
    $command9 .= "\n" . "rm " . $this_nc_denovo_id;
    #-------------------------------------------------------
}





