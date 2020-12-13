## the pipeline includes two major sections: 1) genotype calling from RNAseq samples 2) Allele specific expression analysis with reference bias correction

## 1) Calling genotypes from RNAseq samples 
#perform alignment with Hisat2, convert, and sort

## mark duplicate reads
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load picard/2.20.6

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${line}_merged.bam O=${line}_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${line}_output.metrics 

## use sortsam to sort dedupped bam by coordinate (required)
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load picard/2.20.6

java -jar $EBROOTPICARD/picard.jar SortSam I=${line}_dedupped.bam O=${line}_dedupped_sorted.bam SORT_ORDER=coordinate
## split N cigar reads (reads that span intronic regions)
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module laod gatk/4.1.6.0

gatk SplitNCigarReads --input ${line}_dedupped_sorted.bam --output ${line}_split.bam --reference /scratch/aa9gj/wrk_backup/genomes/gencode_v34_grch38/GRCh38.primary_assembly.genome.fa

##add read groups to split file
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load picard/2.20.6

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${line}_split.bam O=${line}_split_RG.bam RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${line}

##base recalibration to correct errors. Make sure to download the --known-sites file from ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/ 
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load gatk/4.1.6.0

gatk BaseRecalibrator --input ${line}_split_RG.bam --output ${line}_recal_data.table --reference /scratch/aa9gj/wrk_backup/genomes/gencode_v34_grch38/GRCh38.primary_assembly.genome.fa --known-sites /scratch/aa9gj/wrk_backup/genomes/gencode_v34_grch38/common_all_20180418.vcf 

##apply recalibration results from above
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module laod gatk/4.1.6.0

gatk ApplyBQSR --input ${line}_split_RG.bam --output ${line}_recalibrated.bam --reference /scratch/aa9gj/wrk_backup/genomes/gencode_v34_grch38/GRCh38.primary_assembly.genome.fa --bqsr-recal-file ${line}_recal_data.table

##calling snps
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load gatk/4.1.6.0

gatk HaplotypeCaller --input ${line}_recalibrated.bam --output ${line}.vcf --reference /scratch/aa9gj/wrk_backup/genomes/gencode_v34_grch38/GRCh38.primary_assembly.genome.fa --dbsnp /scratch/aa9gj/wrk_backup/genomes/gencode_v34_grch38/common_all_20180418.vcf

##selecting only snps and removing indels because indels cause gatk to throw and error 
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load gatk/4.1.6.0

gatk SelectVariants --variant ${line}.vcf --reference /scratch/aa9gj/wrk_backup/genomes/gencode_v34_grch38/GRCh38.primary_assembly.genome.fa --output ${line}_select.vcf --select-type-to-include SNP

##we ended up with a vcf file

##2) Allele specific expression analysis with reference bias correction

##add readgroups to the original bam files 

#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load picard/2.20.6

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${line}_merged.bam O=${line}_RG.bam RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${line}

##create a reference vcf file for each chromosome in order for WASP to work. use this code to end up with a partitioned file. make sure your input vcf file doesn't have .gz. In my case I put the files in chrom_vcf directory
#!/bin/bash

#SBATCH -J separate
#SBATCH -n 1
#SBATCH -t 20:00:00
#SBATCH -o out
#SBATCH -e err
#SBATCH -p standard
#SBATCH -A cphg-farber
#SBATCH --mem=10000

module load samtools
module load bcftools

VCF="common_all_20180418.vcf"
VCFGZ="${VCF##*/}.gz"  # get basename and add .gz compression extension

bgzip -c $VCF > $VCFGZ  #compress vcf
tabix -p vcf $VCFGZ  # index compressed vcf
tabix --list-chroms $VCFGZ > chromosomes.txt  # save all the chromosome names into a file

while IFS= read -r line; do
  tabix $VCFGZ $line > $line.vcf;
done < chromosomes.txt  # make an individual vcf for each chromosome

##extract snps from resultant vcf files in the chrom_vcf directory (needed as input for WASP)
/scratch/aa9gj/wrk_backup/human_ASE/alignment_RF/merged_all_bam/WASP_corrected_ASE/WASP/mapping/extract_vcf_snps.sh chrom_vcf chrom_vcf

####we did this make sure you have all these options, the output is what's needed for the realignment and downstream analysis. notice we're putting directories for outputs and snp_dir
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load gcc/9.2.0  mvapich2/2.3.3
module load python/3.7.7
import tables

python /scratch/aa9gj/wrk_backup/human_ASE/alignment_RF/merged_all_bam/WASP_corrected_ASE/WASP/mapping/find_intersecting_snps.py --is_paired_end --is_sorted --output ./WASP_analysis --snp_dir ./chrom_vcf /scratch/aa9gj/wrk_backup/human_ASE/alignment_RF/merged_all_bam/${line}_RG.bam

##realign with hisat2 and make sure to convert sort and index. note that the input files here have fq1 fq2 because that's the output from the other command
#!/bin/bash

module load gcc/7.1.0
module load hisat2/2.1.0

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

hisat2 2>./${line}_alignment_summary.txt --rna-strandness RF -x /scratch/aa9gj/wrk_backup/genomes/gencode_v34_grch38/latest/latest -1 $(ls | grep -e ${line} | grep -e "remap.fq1") -2 $(ls | grep -e ${line} | grep -e "remap.fq2") --dta -S ${line}_hisat2.sam


##filter remapped step. Make sure to keep keep.bam from both steps and name them differently. You'll see that we need both in the last merge step below 
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load gcc/9.2.0  mvapich2/2.3.3
module load python/3.7.7
import tables

python /scratch/aa9gj/wrk_backup/human_ASE/alignment_RF/merged_all_bam/WASP_corrected_ASE/WASP/mapping/filter_remapped_reads.py ${line}_RG.to.remap.bam ${line}_wasp.sort.bam ${line}_RG_filter_remap.keep.bam

##merge both files to make a new bam, then sort it and index it. The file is needed as input for gatk ASEreadcounter
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load samtools

samtools merge ${line}.keep.merge.bam ${line}_RG_filter_remap.keep.bam ${line}_RG.keep.bam

##perform ASE read counter
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk_backup/human_ASE/run23_uniq)

module load gatk/4.1.6.0

gatk ASEReadCounter --input ${line}.keep.merge.sort.bam --variant /scratch/aa9gj/wrk_backup/human_ASE/alignment_RF/merged_all_bam/${line}_select.vcf --output ${line}_wasp_counter -R /scratch/aa9gj/wrk_backup/genomes/gencode_v34_grch38/GRCh38.primary_assembly.genome.fa



