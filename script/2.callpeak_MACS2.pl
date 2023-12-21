#!/usr/bin/env perl

if(@ARGV != 3){
    print "%perl " . $0 . " [scaleFactor.txt] [Absolute inputBED_directory] [Absolute outputPeak_directory]\n";
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

$howmany_samples = keys %SAMPLE_IDs;
print "#\t" . $howmany_samples . " SFs: ";
print join(" ", keys %SAMPLE_IDs) . "\n";
print "\n";

### Check input BED files
undef(%BEDS);
for($i=0; $i<scalar(@REP_IDs); $i++){
    $this_ID = $REP_IDs[ $i ];

    #This BED file is the conversion of a BAM file
    #%>bamToBed -i bamfile_of_Bowtie2_mapping | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > output_bed.tmp
    #%>awk -F $'\t' 'BEGIN {OFS = FS}{ if($6=="+"){$2=$2+4}elseif($6=="-"){$3=$3-5} print $0}' output_bed.tmp > this_bedfile

    $bed_file = $input_dir.  $this_ID . "_HG38.bed.gz";
    if(!-e $bed_file){
	print "Error: not found BED file " . $bed_file . "\n";
	print "Aborted\n";
	exit;
    }
    $BEDS{ $this_ID } = $bed_file;
}


$move = int( $SMOOTH_WIN / 2.0);

foreach $key (keys %SAMPLE_IDs){
    $this_SID = $key;
    $this_output_dir = $output_dir . $this_SID . "/";
    if(!-e $this_output_dir){
	system("mkdir -p $this_output_dir");
    }
    
    #make command lines
    $command = "";
    undef(@rep_bedFiles);

    #for each replicate
    foreach $rep (keys %{$SAMPLE_IDs{$key}}){
	$this_ID  = $this_SID . "_" . $rep;
	$this_BED = $BEDS{ $this_ID };
	$this_BED_sw = $this_output_dir . $this_ID . "_" . $SMOOTH_WIN . ".bed.gz";

	if(!-e $this_BED){
	    print "File : " . $this_BED . " was not found\n";
	    print "Aborted\n";
	    exit;
	}

	# Peak call
	$command .= $MACS2 . " " . $MACS2_OPT . " -t " . $this_BED . " -n " . $this_ID . " --outdir " . $this_output_dir;
	$command .= "\n\n";
	push(@rep_bedFiles, $this_BED);

	#shift reads with $SMOOTH_WIN
	$command .= "zcat " . $this_BED . " | awk -F \$'\\t' 'BEGIN {OFS = FS}{ if(\$6 == \"+\"){ t = \$2; \$2 = t - $move; \$3= t + $move }else if(\$6 == \"-\"){ t = \$3; \$2 = t - $move; \$3 = t + $move} print \$0 }' | sort -V -k1,1 -k2,2 | gzip -cf > " . $this_BED_sw;
	$command .= "\n\n";
    }
    # end of each replicate

    # pooled
    $pooled_bed = $this_output_dir . "pooled.bed.gz";
    $command .= "zcat " . join(" ", @rep_bedFiles) . " | sort -V -k1,1 -k2,2 | gzip -cf > " . $pooled_bed;
    $command .= "\n\n";
    # Peak call
    $command .= $MACS2 . " " . $MACS2_OPT . " -t " . $pooled_bed . " -n " . $this_SID . "_pooled" . " --outdir " . $this_output_dir;
    $command .= "\n\n";


    # BED (BED6 + 4 format) narrowPeak
    #1st: chromosome name
    #2nd: start position of peak
    #3rd: end position of peak
    #4th: name of peak
    #5th: integer score for display in genome browser (e.g. UCSC)
    #6th: strand, either "." (=no strand) or "+" or "-"
    #7th: fold-change
    #8th: -log10pvalue
    #9th: -log10qvalue
    #10th: relative summit position to peak start
    
    # Filter and Merge
    undef(@Final_Peaks);
    foreach $rep (keys %{$SAMPLE_IDs{$key}}){
	$this_ID  = $this_SID . "_" . $rep;
	$narrow   = $this_output_dir . $this_ID . "_peaks.narrowPeak";
	$narrow_tmp = $narrow . ".tmp";
	$narrow_flt = $narrow . ".filt.bed.gz";

	$command .= "\n";
	$command .= "if [[ ! -s ".$narrow." ]] || [[ ! -e " . $narrow . " ]]; then" . "\n";
	$command .= " echo 'Aborted'" . "\n";
	$command .= "fi" . "\n";
    
	$command .= "sort -V -k1,1 -k2,2 ". $narrow . " | " . $MERGEB . " -i stdin -c 4 -o distinct";
	$command .= " | awk 'BEGIN{ OFS= \"\\t\"}{\$4=\"Peak_\"NR; print \$0}' > " . $narrow_tmp;
	$command .= "\n\n";

	if(-e $BLACK_LIST){
	    $command .= $INTB . " -v -a " . $narrow_tmp . " -b " . $BLACK_LIST . " | sort -V -k1,1 -k2,2 | gzip -c > " . $narrow_flt;
	    $command .= "\n\n";

	}else{
	    #print "## Black_List file was not found\n";
	    $command .= "sort -V -k1,1 -k2,2 " . $narrow_tmp . " | gzip -c > " . $narrow_flt;
	    $command .= "\n\n";
	}

	push(@Final_Peaks, $narrow_flt);
    }
    # end or each replicaiton

    ## for POOLED
    $temp_this_ID = $this_SID . "_pooled";
    $narrow       = $this_output_dir . $temp_this_ID . "_peaks.narrowPeak";
    $narrow_tmp   = $narrow . ".tmp";
    $narrow_flt   = $narrow . ".filt.bed.gz";
    
    $command .= "sort -V -k1,1 -k2,2 ". $narrow . " | " . $MERGEB . " -i stdin -c 4 -o distinct";
    $command .= " | awk 'BEGIN{ OFS= \"\\t\"}{\$4=\"Peak_\"NR; print \$0}' > " . $narrow_tmp;
    $command .= "\n\n";

    if(-e $BLACK_LIST){
	$command .= $INTB . " -v -a " . $narrow_tmp . " -b " . $BLACK_LIST . " | sort -V -k1,1 -k2,2 | gzip -c > " . $narrow_flt;
	$command .= "\n\n";
    }else{
	$command .= "sort -V -k1,1 -k2,2 " . $narrow_tmp . " | gzip -c > " . $narrow_flt;
	$command .= "\n\n";
    }
    
    
    ######
    $peak_ovl  = $this_output_dir . $this_SID . "_ovl.narrowPeak.bed.gz";
    $peak_final= $this_output_dir . $this_SID . "_final.narrowPeak.bed.gz";
    $peak_pool = $this_output_dir . $this_SID . "_pooled_peaks.narrowPeak.filt.bed.gz";
    
    $peak_rep = $Final_Peaks[0];
    $command .= "\n";
    $command .= "if [[ ! -s ".$peak_rep." ]] || [[ ! -e " . $peak_rep . " ]]; then" . "\n";
    $command .= " echo 'Aborted'" . "\n";
    $command .= "fi" . "\n";

    $command .= "\n";
    $command .= "if [[ ! -s ".$peak_pool." ]] || [[ ! -e " . $peak_pool . " ]]; then" . "\n";
    $command .= " echo 'Aborted'" . "\n";
    $command .= "fi" . "\n";

    $command .= $INTB . " -wo -a " . $peak_pool . " -b " . $peak_rep;
    $command .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$3-\$2; s2=\$7-\$6; if ((\$9/s1 >= 0.5) || (\$9/s2 >= 0.5)){ print \$0}}\' | cut -f 1-4 | sort -V -k1,1 -k2,2 | uniq";

    for($i=1; $i<scalar(@Final_Peaks); $i++){
	$peak_rep = $Final_Peaks[$i];

	$command .= " | " . $INTB . " -wo -a stdin -b " . $peak_rep;
	$command .= " | awk \'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{s1=\$3-\$2; s2=\$7-\$6; if ((\$9/s1 >= 0.5) || (\$9/s2 >= 0.5)){ print \$0}}\' | cut -f 1-4 | sort -V -k1,1 -k2,2 | uniq";

    }
    $command .= " | gzip -c > " . $peak_ovl;
    $command .= "\n\n";


    $command .= "\n";
    $command .= "if [[ ! -s ".$peak_ovl." ]] || [[ ! -e " . $peak_ovl . " ]]; then" . "\n";
    $command .= " echo 'Aborted'" . "\n";
    $command .= "fi" . "\n";

    # the final peak file: chr start end ID #.reads_REP1 #.reads_REP2 ....
    foreach $rep (sort{$a <=> $b} keys %{$SAMPLE_IDs{$key}}){
	$this_ID     = $this_SID . "_" . $rep;
	$this_BED_sw = $this_output_dir . $this_ID . "_" . $SMOOTH_WIN . ".bed.gz";

	if($rep == 1){
	    $command .= $COVB . " -counts -a " . $peak_ovl . " -b " . $this_BED_sw;
	}else{
	    $command .= " | sort -V -k1,1 -k2,2 | " . $COVB . " -counts -a stdin -b " . $this_BED_sw;
	}
    }
    # end of each replicate
    $command .= " | sort -k 5gr,5gr | awk 'BEGIN{ OFS= \"\\t\"}{\$4=\"" . $this_ID . "_\"NR; print \$0}' | sort -V -k1,1 -k2,2 | gzip -c > " . $peak_final;
    $command .= "\n\n";

    
    $command .= "rm -f " . $this_output_dir . "*_peaks.xls" . "\n";
    $command .= "rm -f " . $this_output_dir . "*_summits.bed" . "\n";
    $command .= "rm -f " . $this_output_dir . "*.tmp" . "\n";
    $command .= "rm -f " . $this_output_dir . "*.bdg" . "\n";
    
    #age
    $AGE = $this_output_dir . $this_SID . "_macs2.age";
    $MEM = "32G";
    #$MEM = "8G";
    $JID = $this_SID . "_macs2";
    
    &THIS_QSUB();
    print "## run sh " . $AGE . "\n";
    print "## or" . "\n";
    print "## qsub $AGE" . "\n";
    print "\n";

    #system("qsub -v PATH=\$PATH $SGE");
    
}

exit;



sub THIS_QSUB(){

    open(OUT, ">$AGE");
print OUT <<EOF;
#!/bin/bash
#\$ -S /bin/bash
#\$ -o $this_output_dir
#\$ -e $this_output_dir
#\$ -N $JID
#\$ -v LD_LIBRARY_PATH=""
#\$ -cwd -j y -l s_vmem=$MEM -l mem_req=$MEM


if [ -e $peak_final ]; then
   rm $peak_final
fi


$command


echo "DONE"
EOF
    close(OUT);

}

