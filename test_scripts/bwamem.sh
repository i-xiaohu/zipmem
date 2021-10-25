bwa="/vol1/agis/ruanjue_group/jifahu/software/bwa/bwa"
bwamem="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/build/bwamem"

hg19_index="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/index_HG19/ref"
human_reads1="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring/ERP001775_s5_1.fq.gz"
human_reads2="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring/ERP001775_s5_2.fq.gz"
human_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/spring_ERP001775_s5_1.seeds"
human_sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/spring_ERP001775_s5_1.sam"

ecoli_index="/vol1/agis/ruanjue_group/jifahu/ecoli/ref/ecoli.fa"
ecoli_reads1="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring/SRR1562082_1.fq.gz"
ecoli_reads2="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring/SRR1562082_2.fq.gz"
ecoli_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/spring_SRR1562082_1.seeds"
ecoli_sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/spring_SRR1562082_1.sam"

ts="$(date +%m%d%H%M%S)"

cp ${bwamem} "${bwamem}_${ts}"
bwamem="${bwamem}_${ts}"
pstat "${bwamem}" seeding ${ecoli_index} "/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/4875816.fq" \
	> "/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/4875816.seeds"
pstat "${bwamem}" extend ${ecoli_index} \
	"/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/4875816.seeds" \
	"/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/4875816.fq" \
	> "/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/4875816.sam"
#pstat n "${bwamem}" seeding -t 24 ${hg19_index} ${human_reads1} > ${human_seeds}
#pstat n "${bwamem}" extend -t 16 ${hg19_index} ${human_seeds} ${human_reads1} > ${human_sam}
#pstat n "${bwamem}" extend -t 16 ${ecoli_index} ${ecoli_seeds} ${ecoli_reads1} > ${ecoli_sam}
rm "${bwamem}"