read_id="ERP001775_s5"
in_file1="/vol1/agis/ruanjue_group/jifahu/human/${read_id}_1.fq.gz"
in_file2="/vol1/agis/ruanjue_group/jifahu/human/${read_id}_2.fq.gz"
com_file="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring/${read_id}.spring"
dec_file1="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring/${read_id}_1.fq.gz"
dec_file2="/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring/${read_id}_2.fq.gz"
spring="/vol1/agis/ruanjue_group/jifahu/software/SPRING/build/spring"
threads_n=16
com_node=$(hostname)
com_log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/spring_${read_id}_com.log"
dec_log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/spring_${read_id}_dec.log"


echo "Spring compress ${read_id}, stderr written to ${com_log}"
echo -e "$(date +%Y-%m-%d\ %H:%M:%S): enable ${threads_n} threads on ${com_node}" > "${com_log}"
echo "" >> "${com_log}"
pstat n ${spring} \
	-c \
	-t ${threads_n} \
	-r \
	-g \
	-i ${in_file1} ${in_file2}\
	-o ${com_file} \
	1>> "${com_log}" \
	2>> "${com_log}"


echo "Spring decompress ${read_id}, stderr written to ${dec_log}"
echo -e "$(date +%Y-%m-%d\ %H:%M:%S): enable ${threads_n} threads on ${com_node}" > "${dec_log}"
echo "" 2>> "${dec_log}"
pstat n ${spring} \
	-d \
	-t ${threads_n} \
	-r \
	-g \
	-i ${com_file} \
	-o ${dec_file1} ${dec_file2}\
	1>> "${dec_log}" \
	2>> "${dec_log}"


#pstat n "/vol1/agis/ruanjue_group/jifahu/software/Spring-reorder-only/build/spring-reorder" \
#	-t 12 \
#	--gzipped-input \
#	--gzipped-output \
#	-i "/vol1/agis/ruanjue_group/jifahu/ecoli/SRR1562082_1.fq.gz" \
#	   "/vol1/agis/ruanjue_group/jifahu/ecoli/SRR1562082_2.fq.gz" \
#	-o "/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring_reorder/${read_id}_a1.fq.gz" \
#	   "/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring_reorder/${read_id}_a2.fq.gz"