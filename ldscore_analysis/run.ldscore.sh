#!/bin/bash    
for j in  Multi_tissue_chromatin

do

cat <<EOF > $j.ldscore.lsf

#!/bin/bash
#BSUB -J coga
#BSUB -P acc_LOAD
#BSUB -q premium
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -R "rusage[mem=4000]"
#BSUB -W 08:55
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

cd /sc/orga/projects/COGA/Manav/Manav_6/ldscore_regression/alcohol

module load python/2.7.16
module load py_packages

/sc/orga/projects/COGA/Manav/Manav_6/ldscore_regression/ldsc-master/ldsc.py \
    --h2-cts METAANALYSIS_MVP_AUD_PGC_COGA1.munged.ma.sumstats.gz \
    --ref-ld-chr ../1000G_EUR_Phase3_baseline/baseline. \
    --out AUD_meta_$j \
    --ref-ld-chr-cts $j.ldcts \
    --w-ld-chr ../weights_hm3_no_hla/weights. 

EOF

        bsub <$j.ldscore.lsf
        rm $j.ldscore.lsf

        done
