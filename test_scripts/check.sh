check="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/build/check"

ts="$(date +%m%d%H%M%S)"
cp ${check} "${check}_${ts}"
check="${check}_${ts}"

# SRR1562082_1/SPRING
#zip_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/SRR1562082_1_spring.seeds"
#mem_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/SRR1562082_1_spring.seeds"
# SRR1562082_1/PgRC
zip_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/SRR1562082_1_pgrc.seeds"
mem_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/SRR1562082_1_pgrc.seeds"

# ERP001775_s5_1/SPRING
zip_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/ERP001775_s5_1_spring.seeds"
mem_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/ERP001775_s5_1_spring.seeds"

pstat n "${check}" seed ${zip_seeds} ${mem_seeds}

rm "${check}"