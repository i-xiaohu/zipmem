diff="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/build/study_diff"

read_id="ERP001775_s4_1"
comp="spring"
reads="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/${comp}/${read_id}.reads"
zip_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/${read_id}_${comp}.seeds"
mem_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/${read_id}_${comp}.seeds"
zip_sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/${read_id}_${comp}.sam"
mem_sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/${read_id}_${comp}.sam"
log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/diff_${read_id}_${comp}"

pstat n ${diff} ${reads} ${mem_seeds} ${zip_seeds} ${mem_sam} ${zip_sam} 1> ${log}.out