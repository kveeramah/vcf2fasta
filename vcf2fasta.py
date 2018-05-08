#!/usr/bin/env python
# -*- coding: ASCII -*-

#####This python script takes a vcf file and creates fasta and paml input files for regions specified in a bed file for all samples
#####When encountering a heteroygote, a random allele is chosen. Sequence not noted in the vcf file is inferred from the reference genome
#####Pysam and Pyvcfmust be installed in the version of python used.
#####vcf files must be bgzipped and tabix indexed. Reference genome must also be indexed
#####The region bed file must have the tab seperated fields in the following order: chromosome, start position, end position 
#####Written (poorly) by Krishna Veeramah (krishna.veeramah@stonybrook.edu)

#####usage is ./vcf2fasta.py <target_region_file.bed> <out_stem> <reference_genome>


###import libraries
import string
import vcf
import pysam
from sys import argv
from random import randint

###Input arguments
VCFin=argv[1]# tabindex and bgzipped vcf file
filenameout=argv[2] #stem of output file
region_file=argv[3] #must have the tab seperated fields in the following order chromosome, start position, end position (zero_based)
ref_file=argv[4] #reference genome


###Read region file
file=open(region_file,'r')
regions=file.read()
regions=string.split(regions,'\n')
if regions[-1]=='':
    del(regions[-1])


###Load up reference genome
ref=pysam.FastaFile(ref_file)

###Load up vcf file
vcf_reader = vcf.Reader(open(VCFin, 'r'))

###Get sample names
samples=vcf_reader.samples
paml_names=[]
for g in range(len(samples)):
    if len(samples[g])<10:
        paml_names.append(samples[g]+' '*(10-len(samples[g])))
    else:
        paml_names.append(samples[g][:10])
    
    
nb_samps=len(samples)


###Iterate through regions
for g in range(len(regions)):
    print 'running '+regions[g]
    k=string.split(regions[g])
    chromo=k[0]
    start=int(k[1])
    end=int(k[2])

    base=ref.fetch(chromo,start,end)
    base=base.upper()

    seqs=[]
    for gg in range(nb_samps):
        seqs.append(list(base))


    for record in vcf_reader.fetch(chromo, start, end):
        alleles=[record.REF]
        for gg in range(len(record.ALT)):
            alleles.append(str(record.ALT[gg]))

        pos=record.POS
        index=pos-start-1
        for gg in range(len(record.samples)):
            geno=record.samples[gg]['GT']
            if (geno<>'0/0') and (geno<>None):
                all12=[alleles[int(geno[0])],alleles[int(geno[-1])]]
                change_all=all12[randint(0,1)]
                seqs[gg][index]=change_all

    region_string=string.join(string.split(regions[g],'\t'),"_")

    fileout1=open(filenameout+'_'+region_string+'.fasta','w')
    fileout2=open(filenameout+'_'+region_string+'.paml','w')
    out2='  '+str(len(paml_names))+'\t'+str(len(base))+'\n'
    fileout2.write(out2)
    for gg in range(len(seqs)):
        out1='>'+samples[gg]+'\n'
        for ggg in range(0,len(seqs[gg]),60):
            out1=out1+string.join(seqs[gg][ggg:ggg+60],'')+'\n'
        fileout1.write(out1)
        out2=paml_names[gg]+'  '+string.join(seqs[gg],'')+'\n'
        fileout2.write(out2)
    fileout1.close()


print 'Finished. Laterz'
