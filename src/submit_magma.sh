#InpFIl="test.list"
InpFil="/ifs/scratch/c2b2/dv_lab/jw3514/NBT/GWAS/dat/ukbb.assoc.list"
qsub -t 1-2419 -tc 100 run_magma_annotate_aggreate.sh -i $InpFil
#qsub -t 1-2  run_magma_annotate_aggreate.sh -i test.list
