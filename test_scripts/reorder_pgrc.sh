#read_id="ERP001775_s5_1"
#read_file="/vol1/agis/ruanjue_group/jifahu/human/${read_id}.fq"
#com_file="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/pgrc/${read_id}.pgrc"
#dec_file="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/pgrc/${read_id}.reads"

read_id="SRR1562082_1"
read_file="/vol1/agis/ruanjue_group/jifahu/ecoli/${read_id}.fq"
com_file="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/pgrc/${read_id}.pgrc"
dec_file="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/pgrc/${read_id}.reads"

pgrc="/vol1/agis/ruanjue_group/jifahu/software/PgRC/build/PgRC"
threads_n=16
com_node=$(hostname)
com_log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/pgrc_${read_id}_com.log"
dec_log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/pgrc_${read_id}_dec.log"

echo "PgRC compresses ${read_id}, stderr written to ${com_log}"
echo -e "$(date +%Y-%m-%d\ %H:%M:%S): enable ${threads_n} threads on ${com_node}" > "${com_log}"
pstat n ${pgrc} -t ${threads_n} -i ${read_file} ${com_file} 2>> "${com_log}"

echo "PgRC decompresses ${read_id}.pgrc, stderr written to ${dec_log}"
echo -e "$(date +%Y-%m-%d\ %H:%M:%S): enable ${threads_n} threads on ${com_node}" > "${dec_log}"
pstat n ${pgrc} -t ${threads_n} -d ${com_file} 2>> "${dec_log}"
mv ${com_file}_out ${dec_file}

