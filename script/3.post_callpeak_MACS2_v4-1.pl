#!/usr/local/bin/perl

require("/home/park/PERL/get_numberBowtie.pl");
use File::Basename 'basename', 'dirname';

if(@ARGV != 1){
    print "%perl " . $0 . " [Abs. Path: conf. file]\n";
    exit;
}

require shift(@ARGV);
#-------------------

$qc_number_file = $STAT_DIR . "before_after_QC.txt";
undef(@BAMS);

@SUB=("HYBRID");
$TARGET = "_mapped_sorted.bam";

#--------------------------------------------------
open(IN, $qc_number_file);
while($line = <IN>){
    if($line =~ /^Group/){ next; }

    chop($line);
    undef(@data);
    @data = split(/\t/, $line);

    $GRP = $data[0];
    $SAMPLE = $data[1];

    $this_DIR = $ROOT . $GRP . "/Bowtie2/";
    $BAM       = $SAMPLE .  $TARGET;

    for($i=0; $i < scalar(@SUB); $i++){
        $this_BAM = $this_DIR . $SUB[$i] . "/" . $SAMPLE . "/" . $BAM;
        if(!-e $this_BAM){
            print $this_BAM . " was not found\n";
            exit;
        }
        push(@BAMS, $this_BAM);
    }

}
close(IN);
#--------------------------------------------------
print "#" . scalar(@BAMS) . " were found\n";


#Get Scale Factor
if(!-e $SCALE_FILE){
    print $SCALE_FILE  ." file was not found\n";
    exit;
}

undef(%SCALES);
open(IN, $SCALE_FILE);
while($line=<IN>){
    if($line =~ /^#/){ next; }
    if($line =~ /^Method/){ next; }

    chop($line);
    undef(@data);
    @data = split(/\t/, $line);

    $method = $data[0];
    $grp = $data[1];
    $sample = $data[2];

    $ref = $data[3]; # if ref is HG38, scale factors from DM6 should be used
    $ID = $method . "_" . $grp . "_" . $sample;

    #which scale factor you want to use?
    $sf  = $data[scalar(@data) -1]; # F4
    $SCALES{ $ID } = $sf;
}
close(IN);


$MACS2 = "MACS2_v4-1";
$FDR_mode = 0; # scale factor

for($cnt = 0; $cnt < scalar(@G); $cnt++){
    #group

    for($cnt2 = 0; $cnt2 < scalar(@T); $cnt2++){
	#tool

	$pref= $P[ $cnt2 ];

	for($cnt3 = 0; $cnt3 < scalar(@R); $cnt3++){
	    #ref

	    # DO NC data first
	    $SF_ID1 = $T[ $cnt2] . "_" . $G[$cnt] . "_NC_1";
	    $SF_ID2 = $T[ $cnt2] . "_" . $G[$cnt] . "_NC_2";
	    $REP1_SF_nc = $SCALES{ $SF_ID1 };
	    $REP2_SF_nc = $SCALES{ $SF_ID2 };
		
	    if($REP1_SF_nc eq ""){
		print "SF not found " . $SF_ID1 . "\n";
		exit;
	    }
	    if($REP2_SF_nc eq ""){
		print "SF not found " . $SF_ID2 . "\n";
		exit;
	    }

	    #===========================================#
	    # filtering PEAKS using $MIN_SCALED_READ_CUT
	    $temp  = $ROOT . $G[$cnt] . "/Bowtie2/" . $R[$cnt3] .  "/";
	    $this_nc_peak = $temp . "/" . $MACS2 . "/NC/NC_final.narrowPeak.bed.gz";
	    if(!-e $this_nc_peak){
		print "Not NC exist " . $this_nc_peak . "\n";
		exit;
	    }

	    $this_nc_peak_cut = $temp . "/" . $MACS2 . "/NC/NC_final.narrowPeak.cut.bed.gz";
	    print "Create " . $this_nc_peak_cut . "\n";

	    unlink($this_nc_peak_cut);
	    $this_nc_peak_cut_tmp = $temp . "/" . $MACS2 . "/NC/NC_final.narrowPeak.cut.bed";
	    open(OUT, ">$this_nc_peak_cut_tmp");
	    print OUT "# Scale factor for Rep1 : " . $REP1_SF_nc . "\n";
	    print OUT "# Scale factor for Rep2 : " . $REP2_SF_nc . "\n";
	    print OUT "# Cut if scaled_avg <= " . $MIN_SCALED_READ_CUT . "\n";
	    print OUT "#Chr Start  End ID #.reads_Rep1 #.reads_Rep2 Scaled_Rep1 Scaled_Rep2 Scaled_Avg" . "\n";
	    close(OUT);

	    # IF average scaled reads > MIN_SCALED_READ_CUT, 
	    #chr start end ID raw_rep1 raw_rep2 scaled_rep1 scaled_rep2 scaled_ave
	    $command1 = " zcat " . $this_nc_peak . " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$5 * $REP1_SF_nc; s2=\$6 * $REP2_SF_nc; ss = (s1+s2)/2; if ( ss > $MIN_SCALED_READ_CUT ){ print \$0 \"\\t\" s1 \"\\t\" s2 \"\\t\" ss}}\' >> " . $this_nc_peak_cut_tmp;
	    $command2 = "gzip " . $this_nc_peak_cut_tmp;
	    system($command1);
	    system($command2);
	    #===========================================#


	    # for each SAMPLE (Rep1 and Rep2 at once)
	    foreach $BAM (@BAMS){
		#Get ID of this BAM
		$this_ID = basename($BAM);
		$this_ID =~ s/$TARGET//;

		$this_ID =~ /(.*)\_(\d+)$/;
		$this_ID = $1;
		$this_REP = $2;
		if($this_REP == 2){ next; } ## skip for rep 2 
		if($SKIP_KO{ $this_ID } eq "yes"){ next; }

		#Get the fullpath of this BAM
		$this_DIR = dirname($BAM) . "/";
		
		# pick up MODE
		$temp = dirname($this_DIR);
		$up_mode_dir = dirname($temp);
		
		$MODE = "";
		if( $temp =~ /HYBRID$/){
		    $MODE = "HYBRID";
		}
		# get group ID
		$this_DIR  =~ /.*\/(G\d+)\/Bowtie2.*/;
		$this_GRP = $1;

		if($this_GRP ne $G[$cnt] ){ next; } # skip not target group
		if($MODE ne $R[$cnt3] ){ next; }   # skip other reference mode (take only HYBRID)
		if($this_ID =~ /^NC/){ next; }       # skip NC

		#-----------------------
		# NC sample for this KO sample
		$this_nc_peak_cut = $temp . "/" . $MACS2 . "/NC/NC_final.narrowPeak.cut.bed.gz";
		if(!-e $this_nc_peak_cut){
		    print "Not NC exist " . $this_nc_peak_cut . "\n";
		    exit;
		}
		$SF_ID1 = $T[ $cnt2] . "_" . $this_GRP . "_NC_1";
		$SF_ID2 = $T[ $cnt2] . "_" . $this_GRP . "_NC_2";
		$REP1_SF_nc = $SCALES{ $SF_ID1 };
		$REP2_SF_nc = $SCALES{ $SF_ID2 };
		
		if($REP1_SF_nc eq ""){
		    print "SF not found " . $SF_ID1 . "\n";
		    exit;
		}
		if($REP2_SF_nc eq ""){
		    print "SF not found " . $SF_ID2 . "\n";
		    exit;
		}
		#-----------------------
		
		#======== PEAKs in KO sample
		$target_dir  = $temp . "/" . $MACS2 . "/" . $this_ID . "/";
		if(!-e $target_dir){
		    print "Not exist " . $target_dir . "\n";
		    exit;
		}

		$this_peak = $target_dir . $this_ID . "_final.narrowPeak.bed.gz";
		if(!-e $this_peak){
		    print "Not PEAK exist " . $this_peak . "\n";
		    exit;
		}

		$SF_ID1 = $T[ $cnt2] . "_" . $this_GRP . "_" . $this_ID . "_1";
		$SF_ID2 = $T[ $cnt2] . "_" . $this_GRP . "_" . $this_ID . "_2";
		$REP1_SF = $SCALES{ $SF_ID1 };
		$REP2_SF = $SCALES{ $SF_ID2 };
		
		if($REP1_SF eq ""){
		    print "SF not found " . $SF_ID1 . "\n";
		    exit;
		}
		if($REP2_SF eq ""){
		    print "SF not found " . $SF_ID2 . "\n";
		    exit;
		}


		# read files 
		$REP1_nc_target_bed_sw = $temp . "/" . $MACS2 . "/NC/REP1_150.bed.gz";
		$REP2_nc_target_bed_sw = $temp . "/" . $MACS2 . "/NC/REP2_150.bed.gz";
		if(!-e $REP1_nc_target_bed_sw){
		    print $REP1_nc_target_bed_sw . " was not found\n";
		    exit;
		}
		if(!-e $REP2_nc_target_bed_sw){
		    print $REP2_nc_target_bed_sw . " was not found\n";
		    exit;
		}

		$REP1_target_bed_sw =  $target_dir . "REP1_150.bed.gz";
		$REP2_target_bed_sw =  $target_dir . "REP2_150.bed.gz";
		if(!-e $REP1_target_bed_sw){
		    print $REP1_target_bed_sw . " was not found\n";
		    exit;
		}
		if(!-e $REP2_target_bed_sw){
		    print $REP2_target_bed_sw . " was not found\n";
		    exit;
		}
		#---------------


		#=========
		# STEP 1: filtering PEAKS using $MIN_SCALED_READ_CUT (nc sample was done already)
		$this_peak_cut = $target_dir . $this_ID . "_final.narrowPeak.cut.bed.gz";
		
		$this_peak_cut_tmp2 = $target_dir . $this_ID . "_final.narrowPeak.cut.tmp.bed";
		$this_peak_cut_tmp = $target_dir . $this_ID . "_final.narrowPeak.cut.bed";
		$command3  = "if [ -e " . $this_peak_cut_tmp2 . " ]; then " . "\n";
		$command3 .= " rm " . $this_peak_cut_tmp2 . "\n";
		$command3 .= "fi" . "\n";

		$command3 .= "if [ -e " . $this_peak_cut_tmp . " ]; then " . "\n";
		$command3 .= " rm " . $this_peak_cut_tmp . "\n";
		$command3 .= "fi" . "\n";
		$command3 .= "echo \"# Scale factor for Rep1 : $REP1_SF \" > " . $this_peak_cut_tmp . "\n";
		$command3 .= "echo \"# Scale factor for Rep2 : $REP2_SF \" >> " . $this_peak_cut_tmp . "\n";
		$command3 .= "echo \"# Scale factor for ctrl_Rep1 : $REP1_SF_nc \" >> " . $this_peak_cut_tmp . "\n";
		$command3 .= "echo \"# Scale factor for ctrl_Rep2 : $REP2_SF_nc \" >> " . $this_peak_cut_tmp . "\n";
		$command3 .= "echo \"# Cut if scaled_avg <=  $MIN_SCALED_READ_CUT \" >> " .  $this_peak_cut_tmp . "\n";
		$command3 .= "echo \"#Chr Start  End ID #.reads_Rep1 #.reads_Rep2 Scaled_Rep1 Scaled_Rep2 Scaled_Avg Ctrl_Rep1 Ctrl_Rep2 ctrl_sR1 ctrl_sR2 Ctrl_sR_Avg\" >> " . $this_peak_cut_tmp . "\n\n";


		$command3 .= "zcat " . $this_peak . " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$5 * $REP1_SF; s2=\$6 * $REP2_SF; ss = (s1+s2)/2; if ( ss > $MIN_SCALED_READ_CUT ){ print \$0 \"\\t\" s1 \"\\t\" s2 \"\\t\" ss}}\' > " . $this_peak_cut_tmp2;
		## ADD scaled reads of counter part
		$command3 .= "\n" . "/home/park/TOOL/BEDTools/bin/coverageBed -counts -a " . $this_peak_cut_tmp2 . " -b " . $REP1_nc_target_bed_sw;
                $command3 .= " | /home/park/TOOL/BEDTools/bin/coverageBed -counts -a stdin -b " . $REP2_nc_target_bed_sw;
		$command3 .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$10 * $REP1_SF_nc; s2=\$11 * $REP2_SF_nc; ss = (s1+s2)/2;  print \$0 \"\\t\" s1 \"\\t\" s2 \"\\t\" ss}\' >> " . $this_peak_cut_tmp . "\n\n";

		$command4  = "if [ -e " . $this_peak_cut . " ]; then " . "\n";
		$command4 .= " rm " . $this_peak_cut . "\n";
		$command4 .= "fi" . "\n";
		$command4 .= "gzip " . $this_peak_cut_tmp . "\n";
		$command4 .= "rm " . $this_peak_cut_tmp2 . "\n\n";;


		##$this_nc_peak_cut = $temp . "/MACS2/NC/NC_final.narrowPeak.cut.bed.gz";
		$this_nc_peak_cut2 = $target_dir . $this_ID . "_nc_final.narrowPeak.cut.bed";
		$command4 .= "\n\n" . "if [ -e " . $this_nc_peak_cut2 . " ]; then " . "\n";
		$command4 .= " rm " . $this_nc_peak_cut2 . "\n";
		$command4 .= "fi" . "\n";
		$command4 .= "echo \"# Scale factor for Rep1 : $REP1_SF_nc \" > " . $this_nc_peak_cut2 . "\n";
		$command4 .= "echo \"# Scale factor for Rep2 : $REP2_SF_nc \" >> " . $this_nc_peak_cut2 . "\n";
		$command4 .= "echo \"# Scale factor for ctrl_Rep1 : $REP1_SF \" >> " . $this_nc_peak_cut2 . "\n";
		$command4 .= "echo \"# Scale factor for ctrl_Rep2 : $REP2_SF \" >> " . $this_nc_peak_cut2 . "\n";
		$command4 .= "echo \"# Ctrl is $this_ID \" >> " . $this_nc_peak_cut2 . "\n";
		$command4 .= "echo \"# Cut if scaled_avg <=  $MIN_SCALED_READ_CUT \" >> " .  $this_nc_peak_cut2 . "\n";
		$command4 .= "echo \"#Chr Start  End ID #.reads_Rep1 #.reads_Rep2 Scaled_Rep1 Scaled_Rep2 Scaled_Avg Ctrl_Rep1 Ctrl_Rep2 ctrl_sR1 ctrl_sR2 ctrl_sAvg\" >> " . $this_nc_peak_cut2 . "\n\n";
		
		$command4 .= "\n" . "/home/park/TOOL/BEDTools/bin/coverageBed -counts -a " . $this_nc_peak_cut . " -b " . $REP1_target_bed_sw;
                $command4 .= " | /home/park/TOOL/BEDTools/bin/coverageBed -counts -a stdin -b " . $REP2_target_bed_sw;
		$command4 .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$10 * $REP1_SF; s2=\$11 * $REP2_SF; ss = (s1+s2)/2;  print \$0 \"\\t\" s1 \"\\t\" s2 \"\\t\" ss}\' | sort -V -k1,1 -k 2,2 >> " . $this_nc_peak_cut2 . "\n\n";

		#=========
		# STEP 2: genome wide Fold change 10K-bp (howmany: 308,839) using cut_PEAKs
		$this_peak_cut_10K = $target_dir . $this_ID . "_final.narrowPeak.cut.10K.bed.gz";

		$outbed = $this_peak_cut_10K;
		$winbed = $GL_10K;
		$commond_line  = "";
		&GET_Command_Win();
		$command4_1 = $command_line;


		#=========
		$this_peak_cut_20K = $target_dir . $this_ID . "_final.narrowPeak.cut.20K.bed.gz";
		$outbed = $this_peak_cut_20K;
		$winbed = $GL_20K;
		$commond_line  = "";
		&GET_Command_Win();
		$command4_2 = $command_line;

		#=========
		$this_peak_cut_50K = $target_dir . $this_ID . "_final.narrowPeak.cut.50K.bed.gz";
		$outbed = $this_peak_cut_50K;
		$winbed = $GL_50K;
		$commond_line  = "";
		&GET_Command_Win();
		$command4_3 = $command_line;

		#------------------ All Peaks in KO
		# 1     2      3    4    5                   6                     7                  8                 9
		# chr start end ID #.reads_Rep1 #.reads_Rep2 scaled_Rep1 scaled_Rep2 scaled_avg
		# CUP: common and unique peaks
		# Chr Start  End ID r1 r2 sr1 sr2 savg #.reads_Rep1 #.reads_Rep2 Scaled_Rep1 Scaled_Rep2 Scaled_Avg Ctrl_Rep1 Ctrl_Rep2 ctrl_sR1 ctrl_sR2 Ctrl_sR_Avg

 		# STEP 1: for calc. FC, PV, FDR with NC reads (not considering NC peaks)
		#Chr Start  End ID #.reads_Rep1 #.reads_Rep2 Scaled_Rep1 Scaled_Rep2 Scaled_Avg
		# skip 5,6
		$this_common =  $target_dir . $this_ID . "_common.bed";
		$this_common_tmp =  $target_dir . $this_ID . "_common.tmp.bed";
		$command5  = "if [ -e " . $this_common_tmp . " ]; then " . "\n";
		$command5 .= " rm " . $this_common_tmp . "\n";
		$command5 .= "fi" . "\n";
		# chr start end ID r1 r2 s1 s2 avg r1 r2 c1 c2 avg (c=NC)
		# chr start end ID r1 r2 s1 s2 avg r1 r2 c1 c2 avg (c=KO)
 		$command5 .= $INTB . " -wo -a " . $this_peak_cut . " -b " . $this_nc_peak_cut2;
                $command5 .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$3-\$2; s2=\$17-\$16; if ((\$29/s1 >= 0.5) || (\$29/s2 >= 0.5)){ print \$0}}\' | cut -f 1-28 | sort -V -k1,1 -k2,2 | uniq > " . $this_common_tmp;

                ## one KO peak can include multiple NC peaks
                ## So, get only the list of KO peaks, and calc. NC scaled reads in these (??? NO). Take average!!
                $command5 .= "\n\n" . "if [ -e " . $this_common . " ]; then " . "\n";
                $command5 .= " rm " . $this_common . "\n";
                $command5 .= "fi" . "\n";

		#Chr Start  End ID Scaled_Rep1 Scaled_Rep2 Scaled_Avg Scaled_NC1 Scaled_NC2 Scaled_Avg_NC
                $command5 .= $MERGEB . " -i " . $this_common_tmp . " -c 4,7,8,9,21,22,23 -o distinct,sum,sum,sum,mean,mean,mean -d -1 | sort -V -k1,1 -k2,2 > " . $this_common . "\n";
		$command5 .= "\n" . "rm " . $this_common_tmp;
		#------------------ End of All Peaks in KO

		#STEP 2: add NC unique info. to this_CUP
		$this_nc_common_id   = $target_dir . $this_ID . "_nc_common.id.bed";
		$this_nc_denovo_id     =  $target_dir . $this_ID . "_nc_uniq.id.bed";
		$this_nc_denovo         =  $target_dir . $this_ID . "_nc_uniq.bed";

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
		$this_common_id     = $target_dir . $this_ID . "_common.id.bed";
		$this_ko_denovo_id  = $target_dir . $this_ID . "_ko_uniq.id.bed";
		$this_ko_denovo     =  $target_dir . $this_ID . "_ko_uniq.bed";

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

		$command7 .= "\n\n" . "cut -f 1-4 " . $this_ko_denovo . " > " . $this_ko_denovo_id;
		
		## Calc. PV and FDR
		# now, this_common this_ko_denovo this_nc_denovo
		
		#All in one -> Calc. FDR -> split
		$this_CUP     = $target_dir . $this_ID . "_CUP.bed";
		$this_CUP_fdr = $target_dir . $this_ID . "_CUP.fdr.bed";

		#Chr Start  End ID Scaled_Rep1 Scaled_Rep2 Scaled_Avg Scaled_NC1 Scaled_NC2 Scaled_Avg_NC logFC PV FDR
		$command7 .= "\n\n" . "cat " . $this_common . " " . $this_ko_denovo . " " . $this_nc_denovo . " | sort -V -k1,1 -k2,2 > " . $this_CUP;
		#$command7 .= "\n\n" . "perl " . $ROOT . "script/add_Pvalue.pl " . $this_CUP . " " . $this_CUP_fdr . " \"5,6\" \"8,9\" " . $FDR_mode;

		#-------------------------------------------------------
		#Split
		$this_common_fdr  = $target_dir . $this_ID . "_common.fdr.bed";
		$this_ko_denovo_fdr  = $target_dir . $this_ID . "_ko_uniq.fdr.bed";
		$this_nc_denovo_fdr  = $target_dir . $this_ID . "_nc_uniq.fdr.bed";
		#&Split_Each();

		$command7 .= "\n\n" . "perl " . $ROOT . "script/add_Pvalue.pl " . $this_common . " " . $this_common_fdr . " \"5,6\" \"8,9\" " . $FDR_mode;
		$command7 .= "\n\n" . "perl " . $ROOT . "script/add_Pvalue.pl " . $this_ko_denovo . " " . $this_ko_denovo_fdr . " \"5,6\" \"8,9\" " . $FDR_mode;
		$command7 .= "\n\n" . "perl " . $ROOT . "script/add_Pvalue.pl " . $this_nc_denovo . " " . $this_nc_denovo_fdr . " \"5,6\" \"8,9\" " . $FDR_mode;

		$command7 .= "\n\n" . "cat " . $this_common_fdr . " " . $this_ko_denovo_fdr . " " . $this_nc_denovo_fdr . " | sort -V -k1,1 -k2,2 > " . $this_CUP_fdr;
		

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

		
		#sge
		$SGE = $target_dir . $this_ID . "_pstMACS2.sge";
		$MEM = "32G";
		$JID = $this_ID . "_pstMACS2";
		
		print $target_dir . "\n";
		&THIS_QSUB();
		system("qsub -v PATH=\$PATH $SGE");

	    } # each bams

	} # each Ref
	
    } # each Tool

} # each Group
exit;



sub THIS_QSUB(){

    open(OUT, ">$SGE");
print OUT <<EOF;
#!/bin/bash
#\$ -S /bin/bash
#\$ -o $target_dir
#\$ -e $target_dir
#\$ -N $JID
#\$ -v LD_LIBRARY_PATH=""
#\$ -cwd -j y -l s_vmem=$MEM -l mem_req=$MEM

cd $target_dir

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




#$outbed = $this_peak_cut_20K;
#$winbed = $GL_20K;
#$FDR_mode = 0; #no scale factor
#$commond_line  = "";
#&GET_Command_Win();

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
    $command_line .="\n\n". "perl " . $ROOT . "script/add_Pvalue.pl " . $outbed.".common.bed" . " " . $outbed.".common.fdr.bed" . " \"5,6\" \"8,9\" " . $FDR_mode;

    ## Calc. PV and FDR
    $command_line .="\n\n". "perl " . $ROOT . "script/add_Pvalue.pl " . $outbed.".ko_denovo.bed" . " " . $outbed.".ko_denovo.fdr.bed" . " \"5,6\" \"8,9\" " . $FDR_mode;

    ## Calc. PV and FDR
    $command_line .="\n\n". "perl " . $ROOT . "script/add_Pvalue.pl " . $outbed.".nc_denovo.bed" . " " . $outbed.".nc_denovo.fdr.bed" . " \"5,6\" \"8,9\" " . $FDR_mode;

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





