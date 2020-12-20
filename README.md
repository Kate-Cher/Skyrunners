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

1. For some samples we had several pairs of reads, so this files were merged with [merger.sh](https://github.com/Kate-Cher/Skyrunners/blob/main/merger.sh) script.
2. Quality of raw reads was checked using FastQC (v0.11.9).
3. For processing alignment we used [STAR](https://github.com/alexdobin/STAR) (v2.7) and [GENCODE reference genome Release 36 (GRCh38.p13)](https://www.gencodegenes.org/human/) primary assembly. 

We run following command to generate genome indexes:
```bash
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /path/to/genomeDir \
     --genomeFastaFiles /path/to/genome/fasta \
     --sjdbGTFfile /path/to/annotations.gtf \
     --sjdbOverhang 99
```
Then we run following command to process alignment:
```bash
STAR --genomeDir /path/to/genomeDir \
     --sjdbGTFfile /path/to/annotations.gtf \
     --readFilesCommand zcat \
     --readFilesIn /path/to/read_R1.fastq.gz /path/to/read_R2.fastq.gz  \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 16000000000 \
     --outSAMunmapped Within \
     --outFilterMultimapNmax 1 \
     --quantMode TranscriptomeSAM \
     --runThreadN 8 \
     --outFileNamePrefix "/path/to/out/file."
```
4. The most interesting file from previous step was file.Aligned.toTranscriptome.out.bam because we could count gene and isoform expression using it. We also used [RSEM](https://deweylab.github.io/RSEM/) (v1.3.3) to perform this analysis:
```bash
rsem-calculate-expression --paired-end \
                          --bam \
                          --no-bam-output \
                          -p 8 \
                          /path/to/file.Aligned.toTranscriptome.out.bam \
                          /path/to/genome/index out_file_prefix
```
5. Analysis of differential expressing genes was performed in [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) (v1.30.0) R package. We had two R-scripts: for [gene analysis](https://github.com/Kate-Cher/Skyrunners/blob/main/genes_difexp_analysis/script.r) and for [isoform analysis](https://github.com/Kate-Cher/Skyrunners/blob/main/isoform_difexp_analysis/dif_exp.R).
6. Lists of differential expressing genes and isoforms were analysed using [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/annotate.jsp) and [GeneQuery](https://ctlab.itmo.ru/genequery/searcher/).
7. Cluster analysis was processed in [Phantasus](https://ctlab.itmo.ru/phantasus/).

## Results

We compared three types of conditions in pairs and lists of differentially expressing genes were obtained. More information about functional analysis of genes lists may be found [here](https://docs.google.com/presentation/d/1EzHlLa3YpSoFgdo5eapeyzFIxL91gq_KU4ZYvp_kF3k/edit#slide=id.ga9517ca1e1_0_43). Next is planned a time series analysis and donor effect correction and also we are going to find out more about gene sets responsible for neurodegenerative diseases.
