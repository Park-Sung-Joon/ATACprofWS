#!/usr/bin/env perl

#  all reads are scaled reads
# chr, start, end ID read1 ctrl_read2 fc, fdr
#   common:  "_common.fdr.bed", select 1-4, 7,10,11
# fdr cut <0.01
#   ko_denovo: "_ko_uniq.fdr.bed", select 1-4, 7,10,11
#   nc_denovo: "_nc_uniq.fdr.bed", select 1-4, 7,10,11 

if(@ARGV != 3){
    print "%perl " . $0 . " [scaleFactor.txt] [postMACS2_input_directory] [output_directory]\n";
    exit;
}

require "./ATACprofWS.conf";

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
undef(%SAMPLE_IDs);

open(IN, $scale_file);
while($line = <IN>){
    if($line =~ /^\#/){ next; }
    if($line =~ /_Norm/){ next; }

    chop($line);
    undef(@data);
    @data = split(/\t/, $line);
    
    # sample ID is like "replicate_number"
    $data[0] =~ /(.*)\_(\d+)$/;
    $SID = $1;
    $RID = $2;
    $SAMPLE_IDs{ $SID }{$RID} = 1;
}
close(IN);

$howmany_samples = keys %SAMPLE_IDs;
print "#\t" . $howmany_samples . " samples: ";
print join(" ", keys %SAMPLE_IDs) . "\n";
print "\n";

print "\n";
print "######## Peak Positions\n";
&Peak_Positions();

$WIN =  "10K";
print "######## Windows " . $WIN . "\n";
&Peak_Windows();

exit;

sub Peak_Windows{

    # OutPut directory
    $this_dir = $output_dir. "Peak_Windows_" . $WIN . "/";
    if(!-e $this_dir){
	system(mkdir -p $this_dir);
    }

    $common_dir = $this_dir . "COMMON/";
    if(!-e $common_dir){
	system("mkdir -p $common_dir");
    }

    $denovo_dir = $this_dir . "DENOVO/";
    if(!-e $denovo_dir){
	system("mkdir -p $denovo_dir");
    }
    #------------------------

    $out_numbers = $this_dir . "Numbers.txt";
    open(ON, ">$out_numbers");
    print ON "#Sample\tCommon\tKO_denovo\tNC_denovo\n";

    foreach $SID (sort {$a cmp $b} keys %SAMPLE_IDs){
	if($SID eq $CTRL){ next; }

	$numbers = "";
	# Search files
	$common = $input_dir . $SID . "/" . $SID . "_final.narrowPeak.cut." . $WIN . ".bed.gz.common.fdr.bed";
	$ko_deno = $input_dir . $SID . "/" . $SID .  "_final.narrowPeak.cut." . $WIN . ".bed.gz.ko_denovo.fdr.bed";
	$nc_deno = $input_dir . $SID . "/" . $SID .  "_final.narrowPeak.cut." . $WIN . ".bed.gz.nc_denovo.fdr.bed";
	
	if(!-e $common){
	    print $common . " was not found. Aborted\n"; exit;
	}
	if(!-e $ko_deno){
	    print $ko_deno . " was not found. Aborted\n"; exit;
	}
	if(!-e $nc_deno){
	    print $nc_deno . " was not found. Aborted\n"; exit;
	}

	print $SID . "\n";
	$out_common = $common_dir . $SID . "_" . $WIN . "_common.bed";
	print "\t" . $common . "\n";
	&Get_Cols($common, $out_common, "NULL");

	$line = `grep -v "\#\" $out_common | wc -l`;
	chop($line);
	$numbers .= "\t" . $line;
	
	$out_ko = $denovo_dir . $SID . "_" . $WIN . "_KO_denovo.bed";
	print "\t" . $ko_deno. "\n";
	&Get_Cols($ko_deno, $out_ko, $FDR_THRES);

	$line = `grep -v "\#\" $out_ko | wc -l`;
	chop($line);
	$numbers .= "\t" . $line;
	
	$out_nc = $denovo_dir . $SID . "_" . $WIN . "_NC_denovo.bed";
	print "\t" . $nc_deno . "\n";
	&Get_Cols($nc_deno, $out_nc, $FDR_THRES);

	$line = `grep -v "\#\" $out_nc | wc -l`;
	chop($line);
	$numbers .= "\t" . $line;
	
	print ON $SID . "\t" . $numbers . "\n";
    }
    close(ON);
    print "See " . $out_numbers . "\n";

 }


sub Peak_Positions{
    $HEADER = "#Chr\tStart\tEnd\tID\tKO_Read\tNC_Read\tlog2FC\n";

    # OutPut directory
    $this_dir = $output_dir. "Peak_Positions/";
    if(!-e $this_dir){
	system(mkdir -p $this_dir);
    }

    $common_dir = $this_dir . "COMMON/";
    if(!-e $common_dir){
	system("mkdir -p $common_dir");
    }

    $denovo_dir = $this_dir . "DENOVO/";
    if(!-e $denovo_dir){
	system("mkdir -p $denovo_dir");
    }
    #------------------------

    $out_numbers = $this_dir . "Numbers.txt";
    open(ON, ">$out_numbers");
    print ON "#Sample\tCommon\tKO_denovo\tNC_denovo\n";

    foreach $SID (sort {$a cmp $b} keys %SAMPLE_IDs){
	if($SID eq $CTRL){ next; }

	$numbers = "";

	# Search files
	$common = $input_dir . $SID . "/" . $SID . "_common.fdr.bed";
	$ko_deno = $input_dir . $SID . "/" . $SID . "_ko_uniq.fdr.bed";
	$nc_deno = $input_dir . $SID . "/" . $SID . "_nc_uniq.fdr.bed";
	
	if(!-e $common){
	    print $common . " was not found. Aborted\n"; exit;
	}
	if(!-e $ko_deno){
	    print $ko_deno . " was not found. Aborted\n"; exit;
	}
	if(!-e $nc_deno){
	    print $nc_deno . " was not found. Aborted\n"; exit;
	}
    
	#   common:  "_common.fdr.bed", select 1-4, 7,10,11
	print $SID . "\n";

	# without using FDR threshold
	$out_common = $common_dir . $SID . "_common.bed";
	print "\t" . $common . "\n";
	&Get_Cols($common, $out_common, "NULL");

	$line = `grep -v "\#\" $out_common | wc -l`;
	chop($line);
	$numbers .= "\t" . $line;
	
	$out_ko = $denovo_dir . $SID . "_KO_denovo.bed";
	print "\t" . $ko_deno . "\n";
	&Get_Cols($ko_deno, $out_ko, $FDR_THRES);

	$line = `grep -v "\#\" $out_ko | wc -l`;
	chop($line);
	$numbers .= "\t" . $line;
	
	$out_nc = $denovo_dir . $SID . "_NC_denovo.bed";
	print "\t" . $nc_deno . "\n";
	&Get_Cols($nc_deno, $out_nc, $FDR_THRES);

	$line = `grep -v "\#\" $out_nc | wc -l`;
	chop($line);
	$numbers .= "\t" . $line;
	
	print ON $SID . "\t" . $numbers . "\n";
    }
    close(ON);
    print "See " . $out_numbers . "\n";

}

sub Get_Cols{
    my($infile, $outfile, $fdr_cut) = @_;
    my ($line, $cnt, $fdr, @data, @COLS);

    # this is the case of two replication
    # 0-3: chr start end ID
    # 6: scaled avg read in sample
    # 9: scaled avg read in ctrl
    # 10: log2FC = 6 / 9
    @COLS = (0,1,2,3,6,9,10); 

    #$fdr_col = 12;
    # fdr col is the last one

    open(TOUT, ">$outfile");
    print TOUT "#FDR cut: " . $fdr_cut . "\n";
    print TOUT $HEADER;

    open(TIN, $infile);
    while($line = <TIN>){
	if($line =~ /\#/){ next; }

	chop($line);
	undef(@data);
	@data = split(/\t/, $line);
	
	if($fdr_cut ne "NULL"){
	    if($data[ scalar(@data)-1 ] ne "NA"){
		if($data[ scalar(@data) -1 ] >= $fdr_cut){ next; }
	    }
	}

	$line = ""; #output line
	for($cnt=0; $cnt < scalar(@COLS); $cnt++){
	    $line .= $data[ $COLS[ $cnt ] ] . "\t";
	}
	$line =~ s/\t$//;
	
	print TOUT $line . "\n";
    }
    close(TIN);
    close(TOUT);

}

