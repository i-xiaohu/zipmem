# Minicom is a script that only works in this directory.
cd /vol1/agis/ruanjue_group/jifahu/software/minicom || exit

read_id="ERP001775_s4_1"
read_file="/vol1/agis/ruanjue_group/jifahu/human/${read_id}.fq"
ln -s ${read_file} ${read_id}.fastq
read_file="${read_id}.fastq"
com_node=$(hostname)
threads_n=16
com_log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/minicom_${read_id}_com.log"
dec_log="/vol1/agis/ruanjue_group/jifahu/project/zip_seeding/test_results/minicom_${read_id}_dec.log"

echo "Minicom compresses ${read_id}, stderr written to ${com_log}"
echo -e "$(date +%Y-%m-%d\ %H:%M:%S): enable ${threads_n} threads on ${com_node}" > "${com_log}"
pstat n ./minicom -r ${read_file} -t ${threads_n} 2>> "${com_log}"
if [[ ! -f ${read_id}_comp.minicom ]]; then
	echo "Compression failed, exiting..."
	rm -rf ${read_id}_comp output_* ${read_file}
	exit
fi

echo "Minicom decompresses ${read_id}.minicom, stderr written to ${dec_log}"
echo -e "$(date +%Y-%m-%d\ %H:%M:%S): enable ${threads_n} threads on ${com_node}" > "${dec_log}"
mv ${read_id}_comp.minicom ${read_id}.minicom
pstat n ./minicom -d ${read_id}.minicom -t ${threads_n} 2>> "${dec_log}"
mv ${read_id}_dec.reads ${read_id}.reads
mv ${read_id}.minicom "/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/minicom"
mv ${read_id}.reads   "/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/minicom"

# Removing temporary files if any part crashed down
rm -rf ${read_id}_comp output_* ${read_file}
