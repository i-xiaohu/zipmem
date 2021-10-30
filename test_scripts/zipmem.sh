zipmem="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/build/zipmem"

hg19_index="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/index_HG19/ref"
ecoli_index="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/index_ecoli/ref"

# SPRING/ERP001775_s5_1
#reads="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring/ERP001775_s5_1.fq.gz"
#seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/spring_ERP001775_s5_1.seeds"
#sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/spring_ERP001775_s5_1.sam"
#log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/zipmem_spring_ERP001775_s5_1"

# PgRC/ERP001775_s5_1
reads="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/pgrc/ERP001775_s5_1.reads"
seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/pgrc_ERP001775_s5_1.seeds"
sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/pgrc_ERP001775_s5_1.sam"
log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/zipmem_pgrc_ERP001775_s5_1"

# SPRING/SRR1562082_1
#reads="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring/SRR1562082_1.fq.gz"
#seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/spring_SRR1562082_1.seeds"
#sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/spring_SRR1562082_1.sam"

# PgRC/SRR1562082_1
#reads="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/pgrc/SRR1562082_1.reads"
#seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/pgrc_SRR1562082_1.seeds"
#sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/pgrc_SRR1562082_1.sam"

ts="$(date +%m%d%H%M%S)"
cp ${zipmem} "${zipmem}_${ts}"
zipmem=${zipmem}_${ts}

pstat n "${zipmem}" -t 16 ${hg19_index} ${reads} 1> ${sam} 2> ${log}.log
pstat n "${zipmem}" seeding -t 16 ${hg19_index} ${reads} 1> ${seeds} 2> ${log}_seeding.log
pstat n "${zipmem}" extend -t 16 ${hg19_index} ${seeds} ${reads} 1> ${sam} 2> ${log}_extend.log

#pstat n "${zipmem}" seeding -t 16 ${ecoli_index} ${reads} > ${seeds}
#pstat n "${zipmem}" extend -t 16 ${ecoli_index} ${seeds} ${reads} > ${sam}.mate
#pstat n "${zipmem}" -t 16 ${ecoli_index} ${reads} > ${sam}

rm "${zipmem}"