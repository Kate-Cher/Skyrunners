#! /bin/bash
cd ../../rnaseq/RNAskyrunners/RNA_full/
i=87
while [ $i -gt 68 ]
files=()
do
for var in $(ls | grep -E "^III_CH027${i}_S[[:alnum:]]*_[[:alnum:]]*_R2")
do
files+=( $var )
done
zcat ${files[0]} ${files[1]} | gzip -c > ~/../../rnaseq/RNAskyrunners/alignment_ready/III_CH027${i}_R2_merged.fastq.gz
i=$[ $i - 1 ]
#if [ $i -eq 68 ]
#then
#break
#fi
done

