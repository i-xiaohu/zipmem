# Shell wrapper to stage scripts to be run, protecting them from modifications before or in running.
# If you modify this file, don't forget to run command `dos2unix`

if [[ $# -eq 0 ]]; then
	echo "Usage: ./run_wrapper.sh <script1.sh> [script2.sh] ..."
	exit 0
fi

RED='\033[1;36m'
GREEN='\033[1;32m'
END='\033[0m'

# Copy the scripts to be run
jobs=()
for t in "$@"; do
	if [[ ! -f ${t} ]]; then
		echo -e "${t}" does not exist, skip it...
		continue
	fi
	if [[ ! "${t##*.}" = 'sh' ]]; then
		echo -e "${t}" is not a shell file, skip it...
		continue
	fi
	cpt="NULL"
	for ((i = 1; i <= 20; i++)); do
        if [[ ! -f ${t%.*}_running_$i.sh ]]; then
            cpt="${t%.*}_running_$i.sh"
            break
        fi
    done
    if [[ "${cpt}" = "NULL" ]]; then
    	echo -e There are 20 running jobs of "$t", no more allowed.
    	continue
    fi
	cp "$t" "$cpt"
	dos2unix "$cpt" 2> /dev/null
	chmod -w "$cpt"
	echo -e Script "${t}" is copied to job file "${cpt}"
	jobs+=("$cpt")
done
echo ""

n=${#jobs[*]}
if [[ $n -eq 0 ]]; then
	echo -e "${RED}No job collected${END}"
	exit 1
fi

job_id=1
for job_name in ${jobs[*]}; do
	echo -e "${GREEN}$(date +%Y-%m-%d\ %H:%M:%S) [$job_id/$n] $job_name is running${END}"
	chmod +x "$job_name"
	"./$job_name"
	rm -f "$job_name"
	echo ""
	job_id=$((job_id+1))
done
echo -e "${GREEN}$(date +%Y-%m-%d\ %H:%M:%S) All jobs finished${END}"