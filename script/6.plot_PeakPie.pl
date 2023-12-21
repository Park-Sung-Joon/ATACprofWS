#!/usr/bin/env perl

if(@ARGV != 3){
    print "%perl " . $0 . " [input_dir: Peak positions] [output_dir] [FC:not log2]\n";
    exit;
}

$input_dir = shift(@ARGV);
$out_dir = shift(@ARGV);
$FC = shift(@ARGV);

if(substr($input_dir, -1) ne "/"){
    $input_dir .= "/";
}

if(substr($out_dir, -1) ne "/"){
    $out_dir .= "/";
}
if(!-e $out_dir){
    system("mkdir $out_dir");
}

#------------------- TARGET files: common
$pref= "_common.bed";

# from COMMON
@BEDS = glob($input_dir . "COMMON/*.bed");

$outfile = $out_dir . "Diff_MACSpeaks_" . $FC . "FC.txt";
#$ps_file = $outfile . ".ps";
$pdf_file = $outfile . ".pdf";

$FDR = 0.01;
open(OUT, ">$outfile");
print OUT "\# FDR < " . $FDR . "\n";
print OUT "\# FC > " . $FC . " Up and FC < -1*" . $FC . " Dw\n";
print OUT "\#Sample\tUP_inKO\tOthers\tDW_inKO\tUniq_inKO\tUniq_inNC\n";

$FCcut1 = log($FC) / log(2);
$FCcut2 = $FCcut1 * -1;

for($i=0; $i<scalar(@BEDS); $i++){
    $key = $BEDS[ $i];
    $key =~ /(.*\/)(.*)\_common.bed/;
    $this_dir = $1;
    $this_ID = $2;

    print $this_ID . "\n";
    $ko_deno = $this_dir;
    $ko_deno =~ s/COMMON/DENOVO/;
    $ko_deno .= $this_ID . "_KO_denovo.bed";

    $nc_deno = $this_dir;
    $nc_deno =~ s/COMMON/DENOVO/;
    $nc_deno .= $this_ID . "_NC_denovo.bed";

    if(!-e $ko_deno){
	print $ko_deno . " was not found\n";
	exit;
    }
    if(!-e $nc_deno){
	print $nc_deno . " was not found\n";
	exit;
    }

   
    undef(@UNC); #unchanged : less than abs(FC)
    undef(@UP);   # up in KO
    undef(@DW);  # down in KO
    undef(@UNIQ_KO);
    undef(@UNIQ_NC);

    $howmany_all = 0;
    if( -e $key){
	$howmany_all = `grep -v "\#" $key | wc -l`;
	chop($howmany_all);

	@UP  = `grep -v "\#" $key | awk '{ if(\$7 > $FCcut1){ print \$0 } }' `;
	@DW = `grep -v "\#" $key | awk '{ if(\$7 < $FCcut2){ print \$0 } }' `;
	chop(@UP);
	chop(@DW);
    }

    if(-e $ko_deno){
	@UNIQ_KO  = `grep -v "\#" $ko_deno | awk '{ if(\$7 > $FCcut1){ print \$0 } }' `;
	chop(@UNIQ_KO);
    }

    if(-e $nc_deno){
	@UNIQ_NC  = `grep -v "\#" $nc_deno | awk '{ if(\$7 < $FCcut2){ print \$0 } }' `;
	chop(@UNIQ_NC);
    }

    $howmany_up = scalar(@UP);
    $howmany_dw = scalar(@DW);
    $others = $howmany_all - ($howmany_up+$howmany_dw);
    $howmany_uniq_ko = scalar(@UNIQ_KO);
    $howmany_uniq_nc = scalar(@UNIQ_NC);

    print OUT $this_ID . "\t";
    print OUT $howmany_up . "\t";
    print OUT $others . "\t";
    print OUT $howmany_dw . "\t";
    print OUT $howmany_uniq_ko . "\t";
    print OUT $howmany_uniq_nc . "\n";
}
close(OUT);
&R_Plot( $outfile, $pdf_file);

exit;


sub R_Plot{
    my ($infile , $pdf_file) = @_;
    my $r_cmd  = $pdf_file . ".r.cmd";
    my $r_out  = $pdf_file . ".r.out";

    unlink($pdf_file);
    
open(R, ">$r_cmd");
print R <<EOF;

#postscript("test.ps", horizontal=TRUE)
#pdf("$pdf_file", width=800, height=800)
pdf("$pdf_file", paper="a4")
par(mfrow=c(3,2), cex=1.0)

# Read the input table
d <- read.table("$infile", row.names=1, sep="\\t")

#-----------------------------------
# UP+Uniq_KO
# Common
# DW+Uniq_NC
T <- cbind( d[,1]+d[,4], d[,2], d[,3]+d[,5])
rownames(T) <- rownames(d)
colnames(T)  <- c("KO_Up+KO_Uniq", "non-diff", "KO_Dw+NC_Uniq")
cl <- c("red", "lightgray", "royalblue")

TOTAL <- apply(T, 1, sum)
T <- T / TOTAL * 100.0
idx <- sort(T[,1], index.return=TRUE)\$ix

barplot( t(T[idx,]), beside=F, ylim=c(0,100), las=2,
	 ylab="Fraction (%)",
	 main = "Differential Peak ($FC fc)",
	 sub="FC: $FC FDR: $FDR",
	 col = cl )

plot(T, type="n", ylab="", xlab="", axes=F)
legend("center", 
       legend=colnames(T), col=cl, pch=15,
       cex=1.0, bty="n", xpd=NA)
#-----------------------------------


#-----------------------------------
#UP, DW, Uniq_KO, Uniq_NC
T <- cbind( d[,1], d[,3], d[,4], d[,5])
rownames(T) <- rownames(d)
colnames(T)  <- c("KO_Up", "KO_Down", "KO_Uniq", "NC_Uniq")
cl <- c("red", "royalblue", "orange", "skyblue")

TOTAL <- apply(T, 1, sum)
T <- T / TOTAL * 100.0
idx <- sort(T[,1], index.return=TRUE)\$ix

barplot( t(T[idx,]), beside=F, ylim=c(0,100), las=2,
	 ylab="Fraction (%)",
	 main = "Differential Peaks ($FC fc)",
	 sub="FC: $FC FDR: $FDR",
	 col = cl )

plot(T, type="n", ylab="", xlab="", axes=F)
legend("center", 
       legend=colnames(T), col=cl, pch=15,
       cex=1.0, bty="n", xpd=NA)
#-----------------------------------

#-----------------------------------
#UP, DW, Uniq_KO, Uniq_NC
T <- cbind( d[,1], d[,3], d[,4], d[,5])
rownames(T) <- rownames(d)
colnames(T)  <- c("KO_Up", "KO_Down", "KO_Uniq", "NC_Uniq")
cl <- c("red", "royalblue", "orange", "skyblue")

TOTAL <- apply(T, 1, sum)
idx <- sort(TOTAL, index.return=TRUE)\$ix

barplot( t(T[idx,]/1000), beside=F, las=2,
	 ylab="#. peaks (x1000)",
	 main = "Differential Peaks ($FC fc)",
	 sub="FC: $FC FDR: $FDR",
	 col = cl )

plot(T, type="n", ylab="", xlab="", axes=F)
legend("center", 
       legend=colnames(T), col=cl, pch=15,
       cex=1.0, bty="n", xpd=NA)
#-----------------------------------



#--------- PIE
par(mfrow=c(2,2), cex=1.0)
#T <- cbind( d[,1], d[,3], d[,4], d[,5])
#colnames(T)  <- c("KO_Up", "KO_Down", "KO_only", "NC_only")
#col <- c("red", "royalblue", "orange", "skyblue")

T <- cbind( d[,4], d[,1], d[,5], d[,3])
rownames(T) <- rownames(d)
colnames(T)  <- c("KO_uniq", "KO_Up", "NC_uniq", "KO_Down")
col <- c("orange", "red", "skyblue", "royalblue")

#### PIE
for (i in 1:nrow(T)){

    prot <- rownames(T)[i]
#remove zero
	idx <- T[prot, ] > 0
	T2 <- T[prot, idx]
	this_CL <- col[idx]

#idx <- sort(T2, index.return=TRUE)\$ix
#DATA <- T2[idx]
#CL     <- this_CL[idx]

DATA <- T2
CL     <- this_CL
LB <- paste( names(DATA), " (", DATA , ")", sep="")

pie( DATA, labels=LB, border="white", col=CL, cex=0.8, clockwise=TRUE)
par(new=TRUE)
pie(1, radius=0.5, col="white", border="white", labels="",clockwise=TRUE, main="$FC fc")
text(0,0,labels=prot, cex=par('cex.main'), col=par('col.main'), font=par('font.main'))

## GET log10 scaled_reads
UP_KO <- DW_KO <- KO_uniq <- NC_uniq <- c()
UP_NC <- DW_NC <- c()

file <- paste("$input_dir", "COMMON/", prot, "_common.bed", sep="")

if(file.exists(file)){
BED <- read.table(file, header=FALSE, sep="\\t")

idx      <- BED[, 7] > $FCcut1   #log2FC 
UP_KO <- log10(BED[idx, 5]) #UP in KO and its reads in KO
UP_NC <- log10(BED[idx, 6]) # counter part in NC

idx <- BED[, 7] < $FCcut2
DW_KO <- log10(BED[idx, 5])          #DW in KO
DW_NC <- log10(BED[idx, 6])  # counter part in NC
}


file <- paste("$input_dir", "DENOVO/", prot, "_KO_denovo.bed", sep="")
if( file.exists(file)){
    res <- try ( BED <- read.table(file, header=FALSE, sep="\\t") )
    if( class(res) != "try-error"){
	idx   <- BED[, 7] > $FCcut1   #log2FC 
	if(sum(idx) > 0){
	    KO_uniq <- log10(BED[idx, 5])
        }
    }
}
    
file <- paste("$input_dir", "DENOVO/", prot, "_NC_denovo.bed", sep="")
if( file.exists(file)){
    res <- try (  BED <- read.table(file, header=FALSE, sep="\\t") )
    if( class(res) != "try-error"){
		idx   <- BED[, 7] < $FCcut2   #log2FC 
		if(sum(idx) > 0){
		    NC_uniq <- log10(BED[idx, 6])
	    }
    }
}
  
SP <- c()
boxplot(UP_KO, UP_NC, KO_uniq,  SP,  DW_KO, DW_NC,  NC_uniq,
	main="$FC fc",
	col=c("red", "gray", "orange", "white", "royalblue", "gray", "skyblue"),
	outline =FALSE,
	ylab="Peak Intensity (log10)",
	las=2, 
	axes=F)
axis(2)

  names=c("KO_UP", "NC_DW", "KO_uniq",  "" , "KO_DW", "NC_UP", "NC_uniq")
  text( c(1:length(names)), 0.5, labels=names, xpd=NA, srt=45, cex=0.8 ) 

}

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
