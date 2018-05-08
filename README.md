# vcf2fasta
This python script takes a vcf file and creates fasta and paml input files for regions specified in a bed file for all samples

When encountering a heteroygote, a random allele is chosen. Sequence not noted in the vcf file is inferred from the reference genome

Pysam and Pyvcfmust be installed in the version of python used.

vcf files must be bgzipped and tabix indexed. Reference genome must also be indexed

Usage is ./vcf2fasta.py <target_region_file.bed> <out_stem> <reference_genome>

