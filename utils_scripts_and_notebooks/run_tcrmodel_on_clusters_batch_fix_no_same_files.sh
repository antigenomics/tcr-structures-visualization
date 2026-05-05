#!/bin/bash

#SBATCH --job-name=tcrmodel_clusters        # Job name
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb                 # Job memory request
#SBATCH --time=20-10:30:00           # Time limit hrs:min:sec
#SBATCH --output=JobName.%j.log   # Standard output and error log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luppov.dv@phystech.edu
#SBATCH --partition=infinite
#SBATCH --constraint=gpu
#SBATCH --gres=gpu:1

input_file=$1

if [ ! -f "$input_file" ]; then
    echo "File not found: $input_file"
    exit 1
fi

# Read lines after header into an array
mapfile -t lines < <(tail -n +2 "$input_file")

for line in "${lines[@]}"; do
    job_id=$(echo "$line" | cut -f16)
    tcra_seq=$(echo "$line" | cut -f11)
    tcrb_seq=$(echo "$line" | cut -f12)
    peptide_seq=$(echo "$line" | cut -f8)
    mhc_seq=$(echo "$line" | cut -f14)
    mhc_gen=$(echo "$line" | cut -f17)

    # Полный путь к выходной папке конкретного job_id
    full_output_dir="/projects/structures/clusters/${mhc_gen}/${job_id}"

    # Если папка уже существует — пропускаем
    if [ -d "$full_output_dir" ]; then
        echo "=== SKIPPED: $job_id (output directory already exists) ==="
        continue
    fi

    # Родительская папка для всех job_id с одинаковым mhc_gen
    parent_dir="/projects/structures/clusters/${mhc_gen}"
    mkdir -p "$parent_dir"

    echo "$job_id"
    echo "$peptide_seq"
    echo '-----'
    echo "$mhc_seq"

    python ./run_tcrmodel2.py \
        --job_id="$job_id" \
        --output_dir="$parent_dir" \
        --tcra_seq="$tcra_seq" \
        --tcrb_seq="$tcrb_seq" \
        --pep_seq="$peptide_seq" \
        --mhca_seq="$mhc_seq" \
        --ori_db=/projects/tcrmodel2/
done
