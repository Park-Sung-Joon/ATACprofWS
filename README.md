# ATACprofWS (ATAC-seq Profiler with Spike-in)
## In-house scripts for profiling ATAC-seq peaks examined with Spike-in
## These scripts compare the genome-wide accessibility dynamics caused by a TF KO, along with normalizing read counts and by performing statistical tests.  
## These scripts calculate scale factors and detect common/unique ATAC-seq peaks within a KO sample and a control sample.
<img src="https://github.com/Park-Sung-Joon/ATACprofWS/assets/52985953/b5bc30ec-c30d-4f5a-ab61-b42bee4b2b66" width=650>

## Requirements
BEDTools
macs2

perl, BASH, and R libraries
library(edgeR)

## Preparing
1. Edit **"data/mapped_readCount.txt"** with counting your mapped reads either to the Spikein genome or the Target genome. As an example, we prepared three samples with two replications per each (NC, PARP1, and TFDP1).

```
#Sample	Spike	Target
NC_1	4165657	15361208
NC_2	3753537	12121573
PARP1_1	4652827	14347168
PARP1_2	3585646	13066181
TFDP1_1	2909301	10177613
TFDP1_2	2348308	11196191
```
2. Edit **"ATACprofWS.conf"** with your PC environment. This conf file will be imported by PERL script as follows:
```
require("./ATACprofWS.conf");
```
Therefore, this conf file must be located in the current directory.

## Step 1: calculating scale factors
```
%>perl script/1.scaleFactor.pl data/mapped_readCount.txt scaleFactor.txt
##Reading data/mapped_readCount.txt
##	6 samples read
##	min(total_mapped): 13086914
##	min(spike_mapped): 2348308
##	min(spike_norm): 2268896.62
```
You may see **"scaleFactor.txt"**

## Step 2-1: preparing scripts for peak calling
```
%>perl script/2.callpeak_MACS2.pl scaleFactor.txt /your_directory/data/BED/ /your_directory/MACS2
## Reading scale factors from scaleFactor.txt
#	6 SFs
#	3 SFs: TFDP1 PARP1 NC

## run sh /home/park/Miyanari/ATAC/ATACprofWS/MACS2/TFDP1/TFDP1_macs2.age
## or
## qsub /home/park/Miyanari/ATAC/ATACprofWS/MACS2/TFDP1/TFDP1_macs2.age

## run sh /home/park/Miyanari/ATAC/ATACprofWS/MACS2/PARP1/PARP1_macs2.age
## or
## qsub /home/park/Miyanari/ATAC/ATACprofWS/MACS2/PARP1/PARP1_macs2.age

## run sh /home/park/Miyanari/ATAC/ATACprofWS/MACS2/NC/NC_macs2.age
## or
## qsub /home/park/Miyanari/ATAC/ATACprofWS/MACS2/NC/NC_macs2.age
```
Here, we need BED files that include the location of mapped reads. You may convert BAMs to BED as follows:
```
%>bamToBed -i bamfile_of_Bowtie2_mapping | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > output_bed.tmp
%>awk -F $'\t' 'BEGIN {OFS = FS}{ if($6=="+"){$2=$2+4}elseif($6=="-"){$3=$3-5} print $0}' output_bed.tmp > output.bed
```
We have prepared the example BED files ("data/BED/"), which included the top **500,000** lines only in the original BED files.

## Step 2-2: running the peak-calling scripts 
```
%>sh MACS2/NC/NC_macs2.age
%>sh MACS2/PARP1/PARP1_macs2.age 
%>sh MACS2/TFDP1/TFDP1_macs2.age 
```
This process merges the two replicates after running MACS2. 

## Step 3-1: preparing scripts for capturing final peak sets
```
%>perl script/3.post_callpeak_MACS2_v4-1.pl scaleFactor.txt /your_directory/MACS2 /your_directory/postMACS2_v4-1
```
This process adjusts the read counts at each MACS2 peak using scale Factors.
This process tests statistical significance. The thresholds were defined in **ATACprofWS.conf**. 
This process also generates files for de novo peaks (unique in a KO sample) and common peaks (found in both KO and Control samples).
This process also addresses peaks within 10k-bin genomic windows.

## Step 3-2: running the scripts
```
%>sh postMACS2_v4-1/PARP1/PARP1_pstMACS2.age 
%>sh postMACS2_v4-1/TFDP1/TFDP1_pstMACS2.age
```

## Step 4: gathering all the results to be used for further analysis
```
%>perl script/4.get_Profile.pl scaleFactor.txt /home/park/Miyanari/ATAC/ATACprofWS/postMACS2_v4-1/ /home/park/Miyanari/ATAC/ATACprofWS/Profile
```
Now, the directory **"Profile"** includes the common and de-novo "Peak_Positions" and "Peak_Windows_10K". The file "Numbers.txt" shows the number of ATAC peaks in each category.
For example, "TFDP1_KO_denovo.bed" is for the unique peaks that appeared after TFDP1 KO."TFDP1_NC_denovo.bed" is for the disappeared peaks in the TFDP1-KO sample (= unique peaks in the Control sample).

## Step 5: downstream analysis
```
%>perl script/5.plot_Distribution.pl Profile/Peak_Windows_10K/ NULL Profile/Boxplot
%>perl script/5.plot_Distribution.pl Profile/Peak_Windows_10K/ 2 Profile/Boxplot
%>perl script/6.plot_PeakPie.pl Profile/Peak_Positions/ Profile/Diff_Peaks/ 2
```
By using the files in the directory "Profile/", these scripts look at the genome-wide fold-change of accessibility of KO samples, differential peaks, distribution of unique/common peaks, etc.

## Example: 
```
%>tar xfvz example.tar.gz
```
For guidance purposes, the output results are archived in this example directory.

