#!/usr/bin/env perl

if(@ARGV != 2){
    print "%perl " . $0 . " [input_file: mapped read counts] [output_file]\n";
    exit;
}

# TAB separated file.
# SE: single-end, PE: paired-end
# col_1: sample name
# col_2: read counts in SE or fragment counts in PE,  mapped to spikein
# col_3: read counts in SE or fragment counts in PE, mapped to target
$input_file = shift(@ARGV);
if(!-e $input_file){
    print $input_file . " was not found\n";
    print "Aborted\n";
    exit;
}
$output_file = shift(@ARGV);


# Get the read counts
print "##Reading " . $input_file . "\n";
undef(%COUNTS);
undef(@SAMPLE_IDs);
open(IN, $input_file);
while($line = <IN>){
    if($line =~ /\#/){ next; }

    chop($line);
    undef(@data);
    @data = split(/\t/, $line);

    $COUNTS{ $data[0] }{ "SPIKE" }  = $data[1];
    $COUNTS{ $data[0] }{ "TARGET" } = $data[2];
    if($data[1] eq ""){
	print "Error in reading " . $input_file . " for " . $data[0] . "\n";
	print "Aborted\n";
	exit;
    }
    if($data[2] eq ""){
	print "Error in reading " . $input_file . " for " . $data[0] . "\n";
	print "Aborted\n";
	exit;
    }
	
    push(@SAMPLE_IDs, $data[0]);
}
close(IN);
$howmany = keys %COUNTS;

if( $howmany < 1){
    print "Error in reading " . $input_file . "\n";
    print "Aborted\n";
    exit;
}

print "##\t" . $howmany . " samples read\n";

$min_total = ""; #sequencer output in this group
undef(%CASE);
for($i=0; $i<scalar(@SAMPLE_IDs); $i++){
    $this_ID = $SAMPLE_IDs[$i];
		
    $spike  = $COUNTS{$this_ID}{"SPIKE"};
    $target = $COUNTS{$this_ID}{"TARGET"};
    $total_mapped = $spike + $target;
    if($spike eq ""){
	print "Error from " . $this_ID . "\n";
	print "Aborted\n";
	exit;
    }
    if($target eq ""){
	print "Error from " . $this_ID . "\n";
	print "Aborted\n";
	exit;
    }

    # build a output line
    $CASE{ $this_ID } = $this_ID . "\t" . $spike . "\t" . $target . "\t" . $total_mapped;

    #find global min using total_mapped 
    if($min_total eq ""){ 
	$min_total = $total_mapped;
    }else{
	if($min_total > $total_mapped){
	    $min_total = $total_mapped;
	}
    }
}
print "##\t" . "min(total_mapped): " . $min_total . "\n";

$min_spike = $min_spike_norm = "";
for($i=0; $i<scalar(@SAMPLE_IDs); $i++){
    $this_ID = $SAMPLE_IDs[$i];
    undef(@data);
    @data = split(/\t/, $CASE{ $this_ID } );
    
    $total_mapped = $data[ 3 ];
    $target       = $data[ 2 ];
    $spike        = $data[ 1 ];

    # Calc. factors using min output from seq in this group
    #print $CASE{ $this_ID } . "\n";

    $factor0 = sprintf("%.4f", $total_mapped / $min_total);

    $norm1  = sprintf("%.2f", $spike  / $factor0);
    $norm2  = sprintf("%.2f", $target / $factor0);
    $norm3  = $total_mapped / ( $total_mapped  / $min_total);

    if($min_spike eq ""){
	$min_spike = $spike;
    }else{
	if($min_spike > $spike ){
	    $min_spike = $spike;
	}
    }

    # normalized DM6
    if($min_spike_norm eq ""){
	$min_spike_norm = $norm1;
    }else{
	if($min_spike_norm > $norm1){
	    $min_spike_norm = $norm1;
	}
    }
    
    # build a output line
    $CASE{ $this_ID } .= "\t" . $factor0 . "\t" . $norm1 . "\t" . $norm2 . "\t" . $norm3;
}
print "##\t" . "min(spike_mapped): " . $min_spike . "\n";
print "##\t" . "min(spike_norm): " . $min_spike_norm . "\n";

# Scale Factor
for($i=0; $i<scalar(@SAMPLE_IDs); $i++){
    $this_ID = $SAMPLE_IDs[$i];
    
    undef(@data);
    @data = split(/\t/, $CASE{ $this_ID } );
    $ID           = $data[ 0 ];
    $spike        = $data[ 1 ];
    $target       = $data[ 2 ];
    $total_mapped = $data[ 3 ];
    $factor0      = $data[ 4 ];
    $spike_norm   = $data[ 5 ];

    # based on total output from sequencer
    $factor1  = sprintf("%.4f", $min_spike / $spike);
    $factor2  = sprintf("%.4f",  ($min_spike / $spike) / $factor0 );
    # based on normalized DM6 by factor1
    $factor3    = sprintf("%.4f", $min_spike_norm / $spike_norm);
    # using norm data
    $factor4  = sprintf("%.4f",  ($min_spike_norm / $spike_norm) / $factor0 );

    # build a output line
    $CASE{ $this_ID } .= "\t" . $factor1 . "\t". $factor2 . "\t" . $factor3 . "\t" . $factor4;
}


## Write
open(OUT, ">$output_file");

print OUT "# F0 = ALL_Mapped / min(ALL_Mapped)" . "\n";
print OUT "# Spike_Norm  = Spikein_Mapped / F0" . "\n";
print OUT "# Target_Norm = Target_Mapped  / F0" . "\n";
print OUT "# ALL_Norm    = ALL_Mapped     / F0" . "\n";

print OUT "# F1 = min(Spikein_Mapped) / Spikein_Mapped " . "\n";
print OUT "# F2 = 1 / F0 * F1 " . "\n";
print OUT "# F3 = min( Spikein_Norm) / Spikein_Norm" . "\n";
print OUT "# F4 = 1 / F0 * F3" . "\n";

print OUT "Sample\tSpike\tTarget\tALL\tF0\tSpikein_Norm\tTarget_Norm\tALL_Norm";
print OUT "\tF1\tF2\tF3\tF4" . "\n";

for($i=0; $i<scalar(@SAMPLE_IDs); $i++){
    $this_ID = $SAMPLE_IDs[$i];

    print OUT $CASE{ $this_ID } . "\n";
}
close(OUT);

exit;
