# After seeding pipeline finish, run this script to extend seeds and check alignments

bwamem="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/bin/bwamem"
check="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/bin/check-clip"
vdir="m23"

index="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/index_ecoli/ref"
com_node=$(hostname)
read_id="SRR1562082_1"
comp="spring"
reads="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/${comp}/${read_id}.reads"

zip_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/${read_id}_${comp}.seeds"
mem_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/${read_id}_${comp}.seeds"
zip_sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/${read_id}_${comp}.sam.gz"
mem_sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/${read_id}_${comp}.sam"
zip_log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/${vdir}/zipmem_${read_id}_${comp}"

echo -e "$(date +%Y-%m-%d\ %H:%M:%S): Run on ${com_node}" > ${zip_log}_extend.log
pstat "${bwamem}" extend -t 24 ${index} ${zip_seeds} ${reads} 2>> ${zip_log}_extend.log | gzip > ${zip_sam}

check_log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/${vdir}/check_${read_id}_${comp}"

pstat "${check}" seed ${zip_seeds} ${mem_seeds} 2> ${check_log}_seed.log
pstat "${check}" sam ${zip_sam} ${mem_sam} 2> ${check_log}_sam.log
pstat n "Extend+check_done"