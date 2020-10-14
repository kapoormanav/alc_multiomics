#!/bin/bash

for i in 11:47758531:A:G 11:47397353:G:A

do

#rm -r snp.list.txt
touch snp.list.$i.txt

echo $i >> snp.list.$i.txt


cat <<EOF > run.gcta.${i}.lsf

#!/bin/bash
#BSUB -J coga
#BSUB -P acc_COGA
#BSUB -q premium
#BSUB -n 2
#BSUB -R "span[ptile=2]"
#BSUB -R "rusage[mem=6000]"
#BSUB -W 00:35
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

cd /sc/arion/projects/COGA/Manav/Manav_6/from_scratch_SMR/alcohol_munged/output/conditional_analysis

module load plink2
module load gcta


gcta64  --bfile /sc/arion/projects/LOAD/Dado/projects/2018-11-22.ldref.adgc/output/ADGC_2014.chr11.CPRA_b37  --maf 0.01 --diff-freq 0.5 --cojo-file chr11_spi1_METAANALYSIS_MVP_AUD_PGC_COGA.ma --cojo-cond snp.list.$i.txt --out aud.results.conditioned.on.$i

EOF
        
        bsub <run.gcta.${i}.lsf
        rm run.gcta.${i}.lsf
        #rm snp.list.txt
        
        done
