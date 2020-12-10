# Study of the transcriptome during intense exercises in highlands 

## Project goals and objectives

The aim of this project is to study differential genes expression of 19 sportsmans during physical and psychological stress before and after running in extreme highlands conditions (2450-3450 m, Elbrus m.) and also in "start" point before arrival at competition (
St. Petersburg).

### Tasks

1. Processing and evaluating the quality of raw reads
2. Alignment to the human reference genome ()
3. Count gene and isoform expression level
4. Analysis of differential gene expression
5. Functional analysis of differentially expressing genes
6. Сluster analysis

## Materials and methods

1. For some samples we had several pairs of reads, so this files were merged with merger.sh script.
2. Quality of raw reads was checked using FastQC tool (version ....).
3. For processing alignment we used STAR (version) and [GENCODE reference genome Release 36 (GRCh38.p13)](https://www.gencodegenes.org/human/) primary assembly. 

We run following command to generate genome indexes:
```bash
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /path/to/genomeDir --genomeFastaFiles /path/to/genome/fasta --sjdbGTFfile /path/to/annotations.gtf --sjdbOverhang 99
```
Then we run following command to process alignment:
```bash
STAR --genomeDir /path/to/genomeDir --sjdbGTFfile /path/to/annotations.gtf --readFilesCommand zcat --readFilesIn /path/to/read_R1.fastq.gz /path/to/read_R2.fastq.gz  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN 8 --outFileNamePrefix "/path/to/out/files"
```
