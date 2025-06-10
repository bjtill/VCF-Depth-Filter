# VCF-Depth-Filter (VDF)
A tool to filter a multi-sample VCF by minimum and maximum depth. 
____________________________________________________________________________________________________________

Use at your own risk. I cannot provide support. All information obtained/inferred with this script is without any implied warranty of fitness for any purpose or use whatsoever.

ABOUT: 

This program provides a quick way to filter a multi-sample VCF to retain variants where all samples meet a minimum and maximum depth. Samples with a no-call (./.) are not considered. This is a simpler tool than MVFF (https://github.com/bjtill/MVFF-GUI, https://github.com/bjtill/Multi-sample-Filter-for-Depth-in-VCFs-CLI), but faster and designed to be used as a module in other programs.  Unlike MVFF, VDF provides a maximum depth filter.  

PREREQUISITES:

1. C++ compiler (for example GCC/G++) 

INSTALLATION:

This program should work on all systems. Download the .cpp file and compile according to your system. For Linux Ubuntu, compile with g++ (g++ -o vcf_dpfilter vcf_depth_filter_V2.cpp -lz). 

TO RUN:

The program works with compressed or uncompressed VCFs. Below are command line examples of usage. 

Read compressed VCF, write compressed output

./vcf_dpfilter --input variants.vcf.gz --output filtered.vcf.gz --min-depth 10 --max-depth 100

Read compressed VCF, write uncompressed output

./vcf_dpfilter --input variants.vcf.gz --output filtered.vcf --min-depth 10 --max-depth 100

Read uncompressed VCF, write compressed output

./vcf_dpfilter --input variants.vcf --output filtered.vcf.gz --min-depth 10 --max-depth 100

