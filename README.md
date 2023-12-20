# ATACprofWS (ATAC-seq Profiler with Spike-in)
<img src="https://github.com/Park-Sung-Joon/ATACprofWS/assets/52985953/b5bc30ec-c30d-4f5a-ab61-b42bee4b2b66" width=650>


## under construction.....

## Step 1
```
%>perl script/1.scaleFactor.pl data/mapped_readCount.txt scaleFactor.txt
```

## Step 2
```
%>perl script/2.callpeak_MACS2.pl scaleFactor.txt /your_directory/data/BED/ /your_directory/MACS2
%>sh MACS2/NC/NC_macs2.age 
%>sh MACS2/PARP1/PARP1_macs2.age 
%>sh MACS2/TFDP1/TFDP1_macs2.age 
```

## Step 3
```
%>perl script/3.post_callpeak_MACS2_v4-1.pl scaleFactor.txt /your_directory/MACS2 /your_directory/postMACS2_v4-1
%>sh postMACS2_v4-1/PARP1/PARP1_pstMACS2.age 
%>sh postMACS2_v4-1/TFDP1/TFDP1_pstMACS2.age
```

## Step 4
