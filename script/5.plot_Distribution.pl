#!/usr/bin/env perl

if(@ARGV != 3){
    print "%perl " . $0 . " [input_directory] [abs(FC)|NULL] [output_directory]\n";
    exit;
}

require("./ATACprofWS.conf");

$input_dir  = shift(@ARGV);
$fc_cut   = shift(@ARGV); #not log2
$output_dir = shift(@ARGV);

if( substr($input_dir , -1) ne "/"){
    $input_dir .= "/";
}
if( substr($output_dir , -1) ne "/"){
    $output_dir .= "/";
}

if(!-e $output_dir){
    system("mkdir -p $output_dir");
}

if($fc_cut ne "NULL"){
    # >log2(2) or < -log2(2) 
    $log2_fc_cut = sprintf( "%.1f", log($fc_cut) / log(2) );
}

$tmp1 = $fc_cut;
if($fc_cut eq "NULL"){
    $tmp1 = "all";
}

# Cleaning up
if(!-e $GL_10K){
    print "genome bed file " . $GL_10K . " was not found. Aborted.\n";
    exit;
}

undef(%VALUE_WIN);
open(IN, $GL_10K);
while($line = <IN>){
    chop($line);
    # build ID for this locus
    $ID = $line;
    $ID =~ s/\t/\_/g;

    $VALUE_WIN{ $ID } = $line;
}
close(IN);

$mode = "10K";
$global_profile = $output_dir . "Signals_" . $mode . "_" . $tmp1 . "FC_commonOnly.txt";
#$ps_file = $global_profile . ".ps";
$png_file = $global_profile . ".png";
undef(%VALUE);
%VALUE = %VALUE_WIN;
&CommonOnly();
print "\t" . $global_profile . "\n";
print "\t" . $png_file . "\n";
print "\n";

$global_profile = $output_dir . "Signals_" . $mode . "_" . $tmp1 . "FC_commonAndDenovo.txt";
#$ps_file = $global_profile . ".ps";
$png_file = $global_profile . ".png";
undef(%VALUE);
%VALUE = %VALUE_WIN;
&CommonAndDenovo();
print "\t" . $global_profile . "\n";
print "\t" . $png_file . "\n";
exit;

sub CommonAndDenovo{

    print "#### Common & Denovo\n";

    # from COMMON
    @BEDS = glob($input_dir . "COMMON/*.bed");

    $all_ID = "";
    undef(%COL);
    for($i=0; $i<scalar(@BEDS); $i++){
	$key = $BEDS[ $i];
	$key =~ /(.*\/)(.*)\_$mode\_common.bed/;
	$this_dir = $1;
	$this_ID = $2;
	
	print $this_ID . "\n";
	$ko_deno = $this_dir;
	$ko_deno =~ s/COMMON/DENOVO/;
	$ko_deno .= $this_ID . "_" . $mode. "_KO_denovo.bed";
	
	$nc_deno = $this_dir;
	$nc_deno =~ s/COMMON/DENOVO/;
	$nc_deno .= $this_ID . "_" . $mode. "_NC_denovo.bed";
	
	if(!-e $ko_deno){
	    print $ko_deno . " was not found\n";
	    exit;
	}
	if(!-e $nc_deno){
	    print $nc_deno . " was not found\n";
	    exit;
	}
	
	#check for skip
	#if($SKIP_KO{ $this_ID} eq "yes"){ next; }
	
	$all_ID .= "\t" . $this_ID; #for writing header

	undef(%LOCAL_VALUE);
	open(IN, $key);
	while($line = <IN>){
	    if($line =~ /^\#/){ next; }
	    
	    chop($line);
	    undef(@data);
	    @data = split(/\t/, $line);
	    
	    $ID = $data[0] . "_" . $data[1] . "_" . $data[2] . "_" . $data[3];
	    $fc = sprintf("%.4f", $data[6]); #Fold Change by EdgeR
	    
	    if($fc_cut ne "NULL"){
		if( abs($fc) <= $log2_fc_cut){ next; }
	    }
	    $LOCAL_VALUE { $ID } = $fc;
	}
	close(IN);

	# end of this sample
	foreach $each (keys %VALUE){
	    $fc = $LOCAL_VALUE{ $each };
	    if($fc eq ""){ $fc = 0.0; }
	    $VALUE{ $each } .= "\t" . $fc;
	}
	
    }

    $tmp = $global_profile . ".tmp";
    open(OUT, ">$tmp");
    foreach $key (keys %VALUE){
	print OUT $VALUE{ $key}  . "\n";
    }
    close(OUT);

    open(OUT, ">$global_profile");
    print OUT "Chr\tStart\tEnd\tID" . $all_ID  . "\n";
    close(OUT);

    system("sort -V -k1,1 -k2,2 $tmp >> $global_profile");
    unlink($tmp);

    &R_Plot( $global_profile, $png_file, $log2_fc_cut);
}



sub CommonOnly{
    print "#### Common Only\n";

    # from COMMON
    @BEDS = glob($input_dir . "COMMON/*.bed");

    $all_ID = "";
    undef(%COL);
    for($i=0; $i<scalar(@BEDS); $i++){
	$key = $BEDS[ $i];
	$key =~ /.*\/(.*)\_$mode\_common.bed/;
	$this_ID = $1;

	print $this_ID . "\n";

	#check for skip
	#if($SKIP_KO{ $this_ID} eq "yes"){ next; }
	
	$all_ID .= "\t" . $this_ID; #for writing header

	undef(%LOCAL_VALUE);
	open(IN, $key);
	while($line = <IN>){
	    if($line =~ /^\#/){ next; }
	    
	    chop($line);
	    undef(@data);
	    @data = split(/\t/, $line);
	    
	    $ID = $data[0] . "_" . $data[1] . "_" . $data[2] . "_" . $data[3];
	    $fc = sprintf("%.4f", $data[6]); #Fold Change by EdgeR
	    
	    if($fc_cut ne "NULL"){
		if( abs($fc) <= $log2_fc_cut){ next; }
	    }
	    $LOCAL_VALUE { $ID } = $fc;
	}
	close(IN);
	
	# end of this sample
	foreach $each (keys %VALUE){
	    $fc = $LOCAL_VALUE{ $each };
	    if($fc eq ""){ $fc = 0.0; }
	    $VALUE{ $each } .= "\t" . $fc;
	}
	
    }

    $tmp = $global_profile . ".tmp";
    open(OUT, ">$tmp");
    foreach $key (keys %VALUE){
	print OUT $VALUE{ $key}  . "\n";
    }
    close(OUT);
    
    open(OUT, ">$global_profile");
    print OUT "Chr\tStart\tEnd\tID" . $all_ID  . "\n";
    close(OUT);
    
    system("sort -V -k1,1 -k2,2 $tmp >> $global_profile");
    unlink($tmp);

    &R_Plot( $global_profile, $png_file, $log2_fc_cut);

}


sub R_Plot{
    my ($infile , $png_file, $log2_fc_cut) = @_;
    my $r_cmd  = $png_file . ".r.cmd";
    my $r_out  = $png_file . ".r.out";

#fc_cut = abs(log2 scale) up and down

    unlink($png_file);

open(R, ">$r_cmd");
print R <<EOF;

# config
if( "$fc_cut" == "NULL" ){
     fc_cut <- 0.01    # default resolution. !important
     fc_label <- "all"
}else{
     fc_cut <- $log2_fc_cut
     fc_label <- round($fc_cut, 1)
}
#------------ end of config.


# start
d <- read.table("$infile",header=TRUE, sep="\\t")

# remove all zero cols
idx <- apply( abs( d[, c(5: ncol(d))] ), 1, sum) > 0
UP   <- d[idx, c(5: ncol(d)) ]
DW  <- d[idx, c(5: ncol(d)) ]


UP [ UP <= fc_cut ] <- 0
idx <- apply( abs(UP), 1, sum) > 0
UP <- UP[idx, ]
UP [ UP == 0 ] <- NA #because of calc.ng median
head(UP)

DW [ DW >= -1* fc_cut ] <- 0
idx <- apply( abs(DW), 1, sum) > 0
DW <- DW[idx, ]
DW [ DW == 0 ] <- NA #because of calc.ng median
head(DW)

#postscript("$ps_file", horizontal=TRUE) #FALSE)
png("$png_file", width=800, height=800)
par(mfcol=c(2,2), cex=0.8)

#=============
# keep this order in boxplot
boxMED <- apply(UP, 2, quantile, na.rm=TRUE)["50%",]
mM   <- sort( boxMED, decreasing=FALSE) #TRUE)
bxord   <- names(mM)

boxplot(UP[,bxord], na.rm=TRUE, outline=FALSE,
	las=2,
	col="brown2",
	ylab="log2FCs (Up in KO)",
	main= paste("$mode windows (", fc_label, " FC)", sep="")
)

boxplot(DW[,bxord], na.rm=TRUE, outline=FALSE,
	las=2,
	col="skyblue",
	ylab="log2FCs (Down in KO)"
)
#=============

my.barplot <- function(data, title, lm, fc_label){
    wid <- 0.05
    pos <- barplot(data, las=2,
		 #beside=data,
		   beside=TRUE,
		   space=c(-1,0.5)
		   ,ylim=c(-lm, lm),
		   ,ylab=title,
		   ,col=c("brown2", "skyblue")
		   ,main= paste("$mode windows (", fc_label, " FC)", sep="")
    )

    i<-1
    arrows(pos[i,], data[i,], pos[i,], data[i,]+SD[i,], code=3, angle=90, length=wid)
    i<-2
    arrows(pos[i,], data[i,], pos[i,], data[i,]-SD[i,], code=3, angle=90, length=wid)
}


# mean and sd
# STAT
UPmea <- apply(UP, 2, mean, na.rm=TRUE)
UPsd    <- apply(UP, 2, sd, na.rm=TRUE)
DWmea <- apply(DW, 2, mean, na.rm=TRUE)
DWsd    <- apply(DW, 2, sd, na.rm=TRUE)

ord <- names( sort(UPmea, decreasing=FALSE))

T   <- rbind(UPmea[ord], DWmea[ord])
SD <- rbind(UPsd[ord], DWsd[ord])
lm  <- max( abs(T-SD), T + SD)
head(T)
head(SD)

my.barplot(T, "Mean of log2 FC", lm, fc_label)

#quantile median
UPmea <- apply(UP, 2, quantile, na.rm=TRUE)["50%",]
ord      <- names( sort( UPmea, decreasing=FALSE)  )
UPsd    <- apply(UP, 2, sd, na.rm=TRUE)

DWmea <- apply(DW, 2, quantile, na.rm=TRUE)["50%",]
#ord      <- names( sort( boxMED, decreasing=FALSE)  )
DWsd    <- apply(DW, 2, sd, na.rm=TRUE)

T   <- rbind(UPmea[ord], DWmea[ord])
SD <- rbind(UPsd[ord], DWsd[ord])
lm  <- max( abs(T-SD), T + SD)
my.barplot(T, "Quantile median of log2 FC", lm, fc_label)

dev.off()
q()
EOF

close(R);

    system("R CMD BATCH --no-save $r_cmd $r_out");
    undef(@data);
    @data = `tail -3 $r_out`;
    if( !($data[0] =~ /^\>\sproc\.time.*/) ){
        system("cat $r_out");
    }
    unlink($r_cmd);
    unlink($r_out);
}

