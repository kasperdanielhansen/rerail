#!/bin/bash
set -e -u

# for CHR in $(seq 1 22); do
#     echo "${CHR}"
#     gunzip -c extdata/sra/sraintrons/all_SRA_introns.tsv.gz | grep -P "^chr${CHR}\t" | gzip - > extdata/sra/sraintrons/chr${CHR}_SRA_introns.tsv.gz
# done

for CHR in "X" "Y" "M"; do
    echo "${CHR}"
    gunzip -c extdata/sra/sraintrons/all_SRA_introns.tsv.gz | grep -P "^chr${CHR}\t" | gzip - > extdata/sra/sraintrons/chr${CHR}_SRA_introns.tsv.gz
done
    
