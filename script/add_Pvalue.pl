#!/usr/local/bin/perl

if(@ARGV != 5){
    print "%perl " . $0 . " [input_bed: *.gz] [output_file] [sample: \"1,2\"] [control:\"3,4\"] [Norm?0|1]\n";
    exit;
}

# gz 
$input_bed = shift( @ARGV );
$output_bed = shift( @ARGV ); #pv and fdr added file 
$sample_col = shift( @ARGV );
$control_col = shift( @ARGV );
$norm_mode = shift(@ARGV ); # 0 is no norm, 1 is norm

$output_fit_bed = $output_bed; #edgeR file
$output_fit_bed =~ s/\.bed$//;
$output_fit_bed .= ".edgeR.txt";

$output_bed_tmp      = $output_bed . ".tmp";

$ADJ = "BH"; #p.adjust.method
unlink($output_bed_tmp); # !!! importantc
&Calc_PV( $input_bed, $output_bed_tmp, $output_fit_bed, $sample_col, $control_col, $ADJ);

if(-e $output_bed_tmp){
    # sorting
    system("grep \"\#\" $output_bed_tmp > $output_bed");
    open(OUT, ">>$output_bed");
    print OUT "# B columns: " . $sample_col . "\t" . "A columns: " . $control_col . "\n";
    print OUT "# log2-fold-change (logFC) = log2[ (B+prior) / (A+prior)]" . "\n";
    print OUT "# p.value adjusted (FDR) by " . $ADJ . "\n";
    close(OUT);

    system("grep -v \"\#\" $output_bed_tmp | sort -V -k1,1 -k2,2 >> $output_bed");
    unlink($output_bed_tmp);

    print $output_bed . " was created\n";
    print $output_fit_bed . " was created\n";
}

exit;


sub Calc_PV{
    my ($infile , $outfile, $outfile_fit, $sample_col, $control_col, $ADJ) = @_;
    my $r_cmd  = $outfile . ".r.cmd";
    my $r_out  = $outfile . ".r.out";
    my @data, $i, $cnt;

    undef(@data);
    @data = split(/\,/, $sample_col);
    $howmany_samples = scalar(@data);

    undef(@data);
    @data = split(/\,/, $control_col);
    $howmany_controls = scalar(@data);

    # the first parameter is for sample B (group B)
    $cname = $gname = "";
    for($i=0; $i<$howmany_samples; $i++){
	$cnt = $i+1;
	$cname .= "\"B" . $cnt . "\",";
	$gname .= "\"B\",";
    }
    # the 2nd parameter is for sample A (group A)
    for($i=0; $i<$howmany_controls; $i++){
	$cnt = $i+1;
	$cname .= "\"A" . $cnt . "\",";
	$gname .= "\"A\",";
    }
    $cname =~ s/\,$//;
    $gname =~ s/\,$//;
        
open(R, ">$r_cmd");
print R <<EOF;

library(edgeR)

# prepare table
d <- read.table("$infile",header=FALSE, sep="\\t")
rownames(d) <- d[, 4] #ID in BED file
T <- as.matrix( d[, c( $sample_col, $control_col) ] )  #pick up the target columns
rownames(T) <- d[, 4] #ID in BED file
idx <- apply( abs(T), 1, sum) > 0 #pick up only >0 at least one col.

count <- T[idx, ] #table is "count"
start_B <- 1; stop_B <- start_B -1  + $howmany_samples                #B columns in count table
start_A <- stop_B + 1; stop_A <- start_A -1 + $howmany_controls   #A columns in count table
#start_B <- 1; stop_B <- 2
#start_A <- 3; stop_A <- 4

#colnames(count) <- c("B1", "B2", "A1", "A2")
#grp <- factor( c("B", "B", "A","A") )
colnames(count) <- c($cname)
grp <- factor( c($gname) )
#-------- End of preparing

# Modeling
#remove(dge)

if($norm_mode == 0){
  # lib.size is replaced to 1
  dge <- DGEList( counts=count, group=grp, lib.size=rep(1, ncol(count)) )
  # skip calc. normalization factors
  #dge <- calcNormFactors(dge)
}else{
  dge <- DGEList( counts=count, group=grp )
  dge <- calcNormFactors(dge)
}


design<- model.matrix( ~ grp ) # Full model
dge <- estimateGLMCommonDisp( dge,design ) # common Dispersion
dge <- estimateGLMTagwiseDisp( dge, design)  # Tage Dispersion

#CPMs are calculated: fitted / (lib.size*normfactor) * 1e+6
#log2-FC uses prior.count
#unshrunk.coeffcients is the version without the usage of prior.count
fit     <- glmFit(dge,design)   # fit counts 
res    <- glmLRT(fit, coef=2) # removal test of grpB
# result table includes logFC, LR, Pvalue, FDR
TAB   <- topTags(res, n=nrow(count), adjust.method="$ADJ", sort.by="none", p.value=1 )\$table #Results

# Results in name(=BED IDs) order
name      <- rownames( TAB ); length(name)

fitted  <- fit\$fitted.values[name,] #fitted counts
prior   <- fit\$prior.count #=0.125
colnames(fitted) <- paste("fitted_", colnames(fitted), sep="")

# log2-fold-change using fitted counts, prior, and scaled
edgeR_logFC   <- TAB[, "logFC" ];    length(edgeR_logFC)
names(edgeR_logFC) <- name

edgeR_PV   <- TAB[, "PValue" ];    length(edgeR_PV)
names(edgeR_PV) <- name

# this FDRs are 'p.adjust( edgeR_PV, method="$ADJ")
edgeR_FDR   <- TAB[, "FDR" ];    length(edgeR_FDR)
names(edgeR_FDR) <- name

# average logCPM
# log2-countsPerMillion based on "fitted / (lib.size*normfactor) * 1e+6"
edgeCPM <- TAB[, "logCPM"]; length(edgeCPM)
names(edgeCPM) <- name

# CPM default prior.count = 2?
# Get back 'average' A and B from logFC and logCPM
edgeR_A <- (stop_A - start_A + 1) * ( ( 2^edgeCPM)  / (1+2^edgeR_logFC) ); length(edgeR_A) #average
edgeR_B <- (stop_B - start_B + 1) *  (2^edgeCPM)  - edgeR_A; length(edgeR_B) #average

if($norm_mode == 0){
  edgeR_A <- ( edgeR_A  * 2 ) / 1e+6 # annd think about disper and offset
  edgeR_B <- ( edgeR_B  * 2 ) / 1e+6 # and think about disper and offset
}
#------------- End of gathering edgeR result

#-------------------------
# Based on input raw counts
# Look PVs and FDRs !!! NO questionable??
count_PV  <- apply(count[name,], 1 , function(i){ return( t.test( i[ start_B:stop_B], i[start_A:stop_A])\$p.value) } )
count_FDR<- p.adjust( count_PV[name], method="$ADJ")

count_logFC <- count[name,]
count_avgB <- apply( count_logFC[, c(start_B: stop_B) ], 1, mean); length(count_avgB) #average
count_avgA <- apply( count_logFC[, c(start_A: stop_A) ], 1, mean); length(count_avgA) #average
# Use prior.count. NOTE that this logFC will be exact equal to edgeR_logFC
count_logFC <- log2( (count_avgB+prior) / (count_avgA+prior) ); length(count_logFC)
#-------------------------

# Final Table (raw B,B,A,A avgB avgA logFC Pvalue FDR, fitted B,B,A,A avgB avgA, logFC Pvalue FDR)
# Look !! edgeR_B and edgeR_A are not equal the average of fitted B and A! Why?? cause GLM!
RESTab <- cbind( count[name,], count_avgB, count_avgA, count_logFC, count_PV,count_FDR,
		 fitted, edgeR_B, edgeR_A, edgeR_logFC, edgeR_PV, edgeR_FDR)
rownames(RESTab) <- name

# Write BED file
BED <- d; name <- rownames(BED)
BED <- cbind( BED, round(edgeR_logFC[name], 6), round(edgeR_PV[name], 6), round(edgeR_FDR[name], 6) )
HEADER <- paste( "# col ", ncol(d)+1, ": log2FC",
		 ", col ", ncol(d)+2, ": PValue", 
		 ", col ", ncol(d)+3, ": FDR", sep="")

# counts + FC, Pv, FDR.   !!!!!!! REMOVE $outfile first!!!!
write.table(HEADER, file="$outfile", append=TRUE, quote=FALSE, sep="\t",
	    row.names=FALSE, col.names=FALSE )
write.table(BED, file="$outfile", append=TRUE, quote=FALSE, sep="\t",
	    row.names=FALSE, col.names=FALSE)


# Write out this table
coln <- t( c( "ID", colnames(RESTab)) )

HEADER <- paste(
    "# p.adjust method(=FDR) = $ADJ", 
    "\n# prior(=psuedo) count = ", prior,
    "\n# raw(B,B,A,A) avgB avgA log2FC PV FDR fitted(b,b,a,a) avgb avga log2FC PV FDR",
    sep="")

write.table(HEADER, file="$outfile_fit", append=FALSE, quote=FALSE, sep="\t",
	    row.names=FALSE, col.names=FALSE )

write.table(coln, file="$outfile_fit", append=TRUE, quote=FALSE, sep="\t",
	    row.names=FALSE, col.names=FALSE )

write.table(RESTab, file="$outfile_fit", append=TRUE, quote=FALSE, sep="\t",
	    row.names=TRUE, col.names=FALSE )

q()
EOF

close(R);

    system("R CMD BATCH --no-save $r_cmd $r_out");
    undef(@data);
    @data = `tail -3 $r_out`;
    if( !($data[0] =~ /^\>\sproc\.time.*/) ){
        system("cat $r_out");
    }
    #system("cat $r_out");
    unlink($r_cmd);
    unlink($r_out);
}
