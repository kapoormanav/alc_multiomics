#!/bin/bash

for file in {1..22}

do

#for j in IFN.47231.367.b.PEER_20.expression.sub.ids.txt LPS2.47231.261.b.PEER_20.expression.sub.ids.txt LPS24.47231.322.b.PEER_20.expression.sub.ids.txt
for j in CD14.47231.414.b.PEER_20.expression.sub.ids.txt

do

out=`echo $j | sed 's/.PEER_20.expression.sub.ids.txt//g'`

cat <<EOF> run.mat.eQTL.$file.$out.lsf

#!/bin/bash
#BSUB -J FF_eqtl
#BSUB -P acc_LOAD
#BSUB -q premium
#BSUB -n 1
#BSUB -R span[ptile=1]
#BSUB -R "rusage[mem=15000]"
#BSUB -W 02:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

cd /hpc/users/kapoom02/6/from_scratch_SMR/Fairfax_hrc/EQTL_analysis/peer_expression

module load R

R --vanilla <run.matrix.eqtl.chr$file.$out.R

EOF

        
        
        bsub <run.mat.eQTL.$file.$out.lsf
        rm run.mat.eQTL.$file.$out.lsf
        
        done


done