# Population Genetics Analysis Pipeline
![Update](https://img.shields.io/badge/Update-31/07/2024-green?logo=github)
![Author](https://img.shields.io/badge/Author-Xu.Wang-orange)
![Email](https://img.shields.io/badge/Email-571720850@qq.com-blue?)

## Data preparation
**Input Data**
* WGS data: `.R1.fastq.gz` and `.R2.fastq.gz`
* VCF File

**Dependency**

The detials of all tools can be available in their **offical website** as followed and most of them can quickly install using [Anaconda](https://anaconda.org/):
* [Fastp](https://github.com/OpenGene/fastp), [BWA](https://anaconda.org/bioconda/bwa), [GATK](https://anaconda.org/bioconda/gatk4), [GTX](www.gtxlab.com/product/cat), [Delly](https://anaconda.org/bioconda/delly)
* [VCFtools](https://anaconda.org/bioconda/vcftools), [PLINK2.0](https://anaconda.org/bioconda/plink2)
* [Admixture](https://anaconda.org/bioconda/admixture), [Cytoscape](https://cytoscape.org/)
 
## Variant Calling
### [GTX](www.gtxlab.com/product/cat) (version 4.2.3.0) for SNPs and Indels
This step is to conduct data accusation and remove duplicate reads.
```
i="SRRXXX"
mkdir clean index sort finish gvcf

fastp -i ${i}_1.fastq.gz -o ./clean/${i}_1.clean.fastq.gz \
      -I ${i}_2.fastq.gz -O ./clean/${i}_2.clean.fastq.gz \
      -w 6 --html ./clean/${i}.html

bwa mem -t 8 -R '@RG\tID:"${i}"\tLB:"${i}"\tSM:"${i}"\tPL:illumina' \
        ./index/CSwithmtpt.fa \
        ./clean/${i}_1.clean.fastq.gz ./clean/${i}_2.clean.fastq.gz | samtools sort -O BAM -@ 6 -o ./sort/${i}.sort.bam

gatk MarkDuplicates --REMOVE_DUPLICATES true \
                     -I ./sort/${i}.sort.bam \
                     -M ./finish/${i}.sort.markdup_metrics.txt \
                     -O ./finish/${i}.sort.removedup.bam

samtools index ./finish/${i}.sort.removedup.bam
```

The reference genome for 'Cabernet Sauvignon' has been uploaded to [Zenodo](https://zenodo.org/doi/10.5281/zenodo.8278185), and is named `CSwithmtpt.sort.fa`.
```
gtx index ./index/CSwithmtpt.sort.fa
gtx vc -r ./index/CSwithmtpt.sort.fa \
       -i ./finish/${i}.sort.removedup.bam \
       -o ./gvcf/SRRXXX.g.vcf.gz -g
```

If you have finished variant calling of multiple samples, you can combine them together.
```
gtx joint --reference ./index/CSwithmtpt.fa -t 10 -s ./joint.txt -o all_samples.vcf.gz

-----------
./joint.txt:
B10    /your/path/to/data/gvcf/B10.g.vcf.gz
B11    /your/path/to/data/gvcf/B11.g.vcf.gz
...
```

### [Delly](https://anaconda.org/bioconda/delly) for SVs Calling

The details of SV calling can be found on their official GitHub page: [https://github.com/dellytools/delly](https://github.com/dellytools/delly)

## VCF filtering
Firstly, we should filter the low quality variations using [VCFtools](https://anaconda.org/bioconda/plink2).
```
vcftools --gzvcf ./all_samples.vcf.gz \
         --max-missing 0.8 \
         --minGQ 20 \
         --min-alleles 2 \
         --max-alleles 2 \
         --minDP 4 \
         --maxDP 1000 \
         --maf 0.0001 \
         --recode \
         --out all_miss0.8GQ20maf0.0001
```

Further, we need to add the rs numbers (ID: chr_position) to your VCF files. This will help us filter variations.
* You could find this script `VCF_addID.perl` in this repository.
```
perl VCF_addID.perl all_miss0.8GQ20maf0.0001.vcf all_miss0.8GQ20maf0.0001.id.vcf
```
This **VCF file** `all_miss0.8GQ20maf0.0001.id.vcf` will be the core data for subsequent anlyses, such as polygenetic tree, PCA, admixture, IBD, pi, Fst statistics, introgression (fd statistics), PBScan, GWAS, and more. In this repository, we also provide the example files for plotting.
The filter tool we used in this study is [PLINK2.0](https://anaconda.org/bioconda/plink2) and [VCFtools](https://anaconda.org/bioconda/vcftools).
## Polygenetic tree
```
plink --vcf ./all_miss0.8GQ20maf0.0001.id.vcf \
--vcf-half-call m \
--indep-pairwise 50 5 0.2 \
--geno 0.2 \
--maf 0.0005 \
--recode vcf-iid \
--allow-extra-chr \
--const-fid \
--out 548o_geno0.2gq20maf0.0005LD20.5.0.2

plink --vcf 548o_geno0.2gq20maf0.0005LD20.5.0.2.vcf \
--vcf-half-call m \
--extract 548o_geno0.2gq20maf0.0005LD20.5.0.2.prune.in \
--recode vcf-iid \
--allow-extra-chr \
--const-fid \
--out 548o_geno0.1mq20maf0.01LD20prunein

vcftools --vcf 548o_geno0.1mq20maf0.01LD20prunein.vcf \
--remove-indels \
--recode \
--out 548o_geno0.1mq20maf0.01LD20prunein.rmindel
```
250,821 high quality SNPs and 548 people pass filters and QC. 
* You could find this script `vcf2other.py` in this repository.
```
# vcf2phy
python ./vcf2other.py -f \
-i ./548o_geno0.1mq20maf0.01LD20prunein.rmindel.recode.vcf \
-o 548tree

iqtree -s 548tree.min4.fasta -T 30 -m GTR+I+G -bb 1000 -bnni -alrt 1000 --prefix 548tree.nwk
```
Next, we used iTols [(https://itol.embl.de/)](https://itol.embl.de/) to polish the polygenetic tree.

## PCA analysis
For PCA, we need to remove outgroup to reduce dispersion degree.

```
plink --vcf ./all_miss0.8GQ20maf0.0001.id.vcf \
--allow-extra-chr \
--const-fid \
--threads 1 \
--keep ./keep541.txt \
--pca 20 \
--out 548o_geno0.2gq20maf0.0001

---------------------
./keep541.txt:
0  B1
0  B2
...
```
You can check and download the R script for ploting above in this repository.

## Admixture
Similarly, the admixture analysis also need to remove outgroup.
```
plink --vcf ./all_miss0.8GQ20maf0.0001.id.vcf \
--vcf-half-call m \
--recode 12 \
--keep ./keep541.txt \
--allow-extra-chr \
--const-fid \
--out ./finish

for i in {2..10};do
admixture --cv ./finish.ped $i -j16 | tee log${i}.out ; done

grep -h 'CV' log*.out > CV.out
```
In CV.out, the lowest value of cross-validation error is the best K value. 

```
cat 548tree.nwk.treefile | sed 's/(//g' | sed 's/)//g' | sed 's/:[^,]*,/\n/g' | sed 's/:[^,]*;//g' > treeID.order.txt
```
This code generates the input data according to the tree order for admixture plotting (See example file).

## IBD analysis
Here, we select all seedless grapes for IBD analysis from ***Vitis vinifera*** and ***Vitis vinifera x Vitis labrusca***.
```
plink --vcf ./all_miss0.8GQ20maf0.0001.id.vcf \
--allow-extra-chr \
--const-fid \
--threads 1 \
--keep ./keep46.txt \
--genome \
--out IBD_result

#  ibd.genome  Output
...
FID1           IID1 FID2           IID2 RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO
   0             B7   0            GB1 OT     0  0.1263  0.5415  0.3322  0.6030  -1  0.863903  1.0000 21.2571
   0             B7   0            GB6 OT     0  0.0000  0.0430  0.9570  0.9785  -1  0.992082  1.0000      NA
   0             B7   0     SRR7662541 OT     0  0.8044  0.1503  0.0453  0.1204  -1  0.740131  0.6952  2.0803
...
```
PI_HAT = 1.0 indicates that the two samples are either the same variety or somatic variants.
PI_HAT = 0.5 indicates that the two samples have a parent-offspring relationship.
PI_HAT < 0.5 indicates a distant relative relationship.

The visualization of IBD-PPI network using [Cytoscape](https://cytoscape.org/). GitHub: [https://github.com/cytoscape/cytoscape](https://github.com/cytoscape/cytoscape).

## Fst and pi statistics
We should prepare the ID of two groups (VV, VVxVL) for calculation.

```
# Fst
vcftools --vcf ./all_miss0.8GQ20maf0.0001.id.vcf --weir-fst-pop keep_VV35.txt --weir-fst-pop keep_VVVL11.txt --fst-window-size 20000 --fst-window-step 20000 --out fst_VV35vsVVVL11.txt

# Pi
vcftools --vcf ./all_miss0.8GQ20maf0.0001.id.vcf --keep keep_VV35.txt --window-pi 20000 --window-pi-step 20000 --out pi.keep_VV35.txt
vcftools --vcf ./all_miss0.8GQ20maf0.0001.id.vcf --keep keep_VVVL11.txt --window-pi 20000 --window-pi-step 20000 --out pi.keep_VVVL11.txt
```














