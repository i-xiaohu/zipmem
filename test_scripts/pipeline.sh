# Zipmem version 5.2-ram
# BWA-MEM version 0.7.17-r1188
zipmem="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/bin/zipmem-5.2-ram"
bwamem="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/bin/bwamem"

index="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/index_HG19/ref"
read_id="ERP001775_s4_1"
comp="pgrc"
reads="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/${comp}/${read_id}.reads"
threads_n=12

zip_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/${read_id}_${comp}.seeds"
zip_sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/${read_id}_${comp}.sam"
zip_log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/5_2_ram/zipmem_${read_id}_${comp}"
pstat n "${zipmem}" -t ${threads_n} ${index} ${reads} 1> ${zip_seeds} 2> ${zip_log}_seeding.log
pstat n "${bwamem}" extend -t 20 ${index} ${zip_seeds} ${reads} 1> ${zip_sam} 2> ${zip_log}_extend.log

mem_seeds="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/${read_id}_${comp}.seeds"
mem_sam="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/bwamem/${read_id}_${comp}.sam"
mem_log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/bwamem_${read_id}_${comp}"
pstat n "${bwamem}" seeding -t ${threads_n} ${index} ${reads} 1> ${mem_seeds} 2> ${mem_log}_seeding.log
pstat n "${bwamem}" extend -t 20 ${index} ${mem_seeds} ${reads} 1> ${mem_sam} 2> ${mem_log}_extend.log

# Please check the consistency of seeds and final alignments.