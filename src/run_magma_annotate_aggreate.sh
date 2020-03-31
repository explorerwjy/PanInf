#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N MAGMA
#$ -l h_rt=4:00:00
#$ -l h_vmem=5G
#$ -cwd

while getopts i:r:a:l:t:PFH opt; do
    case "$opt" in
        i) ListFil="$OPTARG";;
        a) ArrNum="$OPTARG";;
        l) LogFil="$OPTARG";;
    esac
done

if [[ -z "${ArrNum}" ]]
then
    ArrNum=$SGE_TASK_ID
fi

SNPLOC="/ifs/scratch/c2b2/dv_lab/jw3514/NBT/GWAS/src/makeSNPLOC.py"
GENELOC="/ifs/scratch/c2b2/dv_lab/jw3514/NBT/GWAS/dat/MAGMA_DATA/NCBI37.3.gene.name.loc"
ONEKGEUR="/ifs/scratch/c2b2/dv_lab/jw3514/NBT/GWAS/dat/MAGMA_DATA/g1000_eur"

InpFil=`readlink -f $(tail -n+$ArrNum $ListFil | head -n 1 | cut -f1)`
NAME=$(basename $InpFil|sed s/.assoc.tsv.gz//g)

#gunzip
gunzip -c $InpFil > $NAME.tsv

#SNPLOC
python $SNPLOC --id rsid --var variant -i $InpFil -o $NAME.snploc

#Annotation
magma --annotate window=5,5 --snp-loc $NAME.snploc --gene-loc $GENELOC --out magma.$NAME

#Gene Test
magma --bfile $ONEKGEUR --gene-annot magma.$NAME.genes.annot --pval $NAME.tsv use=2,9 ncol=nCompleteSamples --out $NAME.genetest

# Remove Tmp files
