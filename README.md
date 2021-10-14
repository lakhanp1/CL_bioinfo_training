# ChIPseq data processing training

## Install required softwares using conda package manager

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n omics_py37 python=3.7
conda deactivate
conda activate omics_py37
conda install -y perl=5.26.2 perl-lwp-simple perl-encode-locale 
conda install -c bioconda -y htslib=1.10.2  samtools=1.10 bcftools=1.10.1
conda install -c bioconda -y bowtie2=2.3.5.1 bwa=0.7.17  hisat2=2.2.0 stringtie=2.1.2
conda install -c bioconda -y bedtools=2.29.2 deeptools=3.4.3
conda install -c bioconda -y ucsc-bedsort ucsc-bedtobigbed ucsc-bedtopsl ucsc-bedgraphtobigwig
conda install -c bioconda -y macs2=2.2.7.1

```

### Process TF ChIPseq data

```bash

mkdir TF1_OE_16h_HA_ChIPMix64_3
cd TF1_OE_16h_HA_ChIPMix64_3

## mapping raw data using bowtie2
bowtie2 -p 2 --trim5 8 --local  -x ~/20211013_CL_bioinfo_training/database/A_nidulans_FGSC_A4/bowtie2_index/A_nidulans_FGSC_A4_version_s10-m04-r03_chromosomes.fasta -U ~/20211013_CL_bioinfo_training/raw_data/TF1_OE_16h_HA_ChIPMix64_3_R1.fastq.gz | samtools view -bS - | samtools sort  -O bam -o TF1_OE_16h_HA_ChIPMix64_3_bt2.bam

## extract mapping stats
samtools index TF1_OE_16h_HA_ChIPMix64_3_bt2.bam
samtools flagstat TF1_OE_16h_HA_ChIPMix64_3_bt2.bam > alignment.stats

mappedReads=`grep -P ' 0 mapped \(' alignment.stats | grep -P -o '^\d+'`
scale=`perl -e "printf('%.3f', 1000000/$mappedReads)"`

##macs2 pileup with 200bp extension
macs2 pileup --extsize 200 -i TF1_OE_16h_HA_ChIPMix64_3_bt2.bam -o TF1_OE_16h_HA_ChIPMix64_3_pileup.bdg

##normalize
printf "Normalizing TF1_OE_16h_HA_ChIPMix64_3_pileup.bdg with factor %s\n" $scale
macs2 bdgopt -i TF1_OE_16h_HA_ChIPMix64_3_pileup.bdg -m multiply -p $scale -o temp_normalized.bdg

##Remove the first line
sed -n '2,$p' temp_normalized.bdg > TF1_OE_16h_HA_ChIPMix64_3_normalized.bdg
rm temp_normalized.bdg

##bedSort and bedGraph to bigWig conversion
bedSort TF1_OE_16h_HA_ChIPMix64_3_normalized.bdg TF1_OE_16h_HA_ChIPMix64_3_normalized.bdg
bedGraphToBigWig TF1_OE_16h_HA_ChIPMix64_3_normalized.bdg ~/20211013_CL_bioinfo_training/database/A_nidulans_FGSC_A4/reference/genome.size TF1_OE_16h_HA_ChIPMix64_3_normalized.bw


##macs2 narrowPeak calling, without using input:
macs2 callpeak -t TF1_OE_16h_HA_ChIPMix64_3_bt2.bam --name TF1_OE_16h_HA_ChIPMix64_3.macs2 --outdir macs2_withoutCtrl_narrow -g 30e6 --nomodel --extsize 200 -B --SPMR

```


### Process TF ChIPseq data
```bash

bowtie2 -p 2  --trim5 8 --local  -x ~/20211013_CL_bioinfo_training/database/A_nidulans_FGSC_A4/bowtie2_index/A_nidulans_FGSC_A4_version_s10-m04-r03_chromosomes.fasta -U ~/20211013_CL_bioinfo_training/raw_data/WT_16h_polII_ChIPMix66_3_R1.fastq.gz | samtools view -bS - | samtools sort  -O bam -o WT_16h_polII_ChIPMix66_3_bt2.bam

##reference config
conf_polIIFeatures=""
##

## mapping stats
samtools index WT_16h_polII_ChIPMix66_3_bt2.bam
samtools flagstat WT_16h_polII_ChIPMix66_3_bt2.bam > alignment.stats

mappedReads=`grep -P ' 0 mapped \(' alignment.stats | grep -P -o '^\d+'`
scale=`perl -e "printf('%.3f', 1000000/$mappedReads)"`

##macs2 pileup with 200bp extension
macs2 pileup --extsize 200 -i WT_16h_polII_ChIPMix66_3_bt2.bam -o WT_16h_polII_ChIPMix66_3_pileup.bdg -f BAM

##normalize
printf "Normalizing WT_16h_polII_ChIPMix66_3_pileup.bdg with factor %s\n" $scale
macs2 bdgopt -i WT_16h_polII_ChIPMix66_3_pileup.bdg -m multiply -p $scale -o temp_normalized.bdg

##Remove the first line
sed -n '2,$p' temp_normalized.bdg > WT_16h_polII_ChIPMix66_3_normalized.bdg
rm temp_normalized.bdg

##bedSort and bedGraph to bigWig conversion
bedSort WT_16h_polII_ChIPMix66_3_normalized.bdg WT_16h_polII_ChIPMix66_3_normalized.bdg
bedGraphToBigWig WT_16h_polII_ChIPMix66_3_normalized.bdg ~/20211013_CL_bioinfo_training/database/A_nidulans_FGSC_A4/reference/genome.size WT_16h_polII_ChIPMix66_3_normalized.bw

perl ~/20211013_CL_bioinfo_training/scripts/zqWinSGR-v2.pl -feature_file ~/20211013_CL_bioinfo_training/database/A_nidulans_FGSC_A4/annotation/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_Unique.bed -socre_file WT_16h_polII_ChIPMix66_3_normalized.bdg -chrom_column 1 -start_column 2 -end_column 3  -direction_column 6 -bin_count 1 -output_folder $PWD -outout_name WT_16h_polII_ChIPMix66_3_polii_expr.tab


```
