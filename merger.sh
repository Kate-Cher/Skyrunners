#! /bin/bash
cd ../../rnaseq/RNAskyrunners/RNA_full/
i=87
while [ $i -gt 68 ]
files=()
do
for var in $(ls | grep -E "^[condition_number]_CH027${i}_S[[:alnum:]]*_[[:alnum:]]*_R1/2")
do
files+=( $var )
done
zcat ${files[0]} ${files[1]} | gzip -c > ~/../../rnaseq/RNAskyrunners/alignment_ready/[condition_number]_CH027${i}_R1/2_merged.fastq.gz
i=$[ $i - 1 ]
done

