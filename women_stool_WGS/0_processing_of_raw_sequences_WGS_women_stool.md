# Processing of Raw WGS Stool Samples (Mothers) — Preeclampsia Project
# University of Massachusetts — HPC
# ~100 samples (mothers only)
# Pipeline: KneadData → MetaPhlAn4 (absolute & relative abundance) → GTDB taxonomy conversion → Diversity analysis

# -----------------------------
# 1. Connect to the HPC
ssh user@hpc-server-address

# -----------------------------
# 2. Download data from Illumina BaseSpace to the HPC
module load basespace-cli/
date=$(date +"%Y%m%d")
bsub -q long -W 16:00 -o bs_download_${date}.log \
  "bs list dataset --input-run <RUN_ID> --terse | xargs -n1 -I{} bs dataset download --id {} -o {}"

# -----------------------------
# 3. (Optional) Count FASTQ files
find . -type f -name "*.fastq.gz" | wc -l

# -----------------------------
# 4. Gather R1 FASTQ file paths
find /path/to/WGS_mom_raw -type f -name "*_R1_*.fastq.gz" > file_moms_R1.txt
# R2 will be paired automatically by the scripts.

# -----------------------------
# 5. Run KneadData (human read removal + QC)
# Script is a wrapper for KneadData; adjust memory if needed.
kneaddata_PARTA_human_fromLIST.sh -i file_moms_R1.txt -o output_kneaddata_moms

# -----------------------------
# 6. Collect KneadData output paths
find /path/to/output_kneaddata_moms/*_knd/*kneaddata.fastq.gz \
  -type f > file_moms_knead_output.txt

# -----------------------------
# 7. Run MetaPhlAn4 for absolute counts
metaphlan4_fromList_counts.sh -i file_moms_knead_output.txt \
  -o output_metaphlan_abs_counts_moms

# -----------------------------
# 8. Run MetaPhlAn4 for relative abundance
metaphlan4_fromList_rel_ab.sh -i file_moms_knead_output.txt \
  -o output_metaphlan_rel_ab_moms

# -----------------------------
# 9. Transfer results to local machine
scp -r user@hpc-server-address:/path/to/output_metaphlan_rel_ab_moms/*.txt .
scp -r user@hpc-server-address:/path/to/output_metaphlan_abs_counts_moms/*.txt .

# -----------------------------
# 10. Convert MetaPhlAn4 SGB profiles to GTDB taxonomy
for file in /path/to/rel_ab/*.mph.txt; do
    base_name=$(basename "$file" .mph.txt)
    sgb_to_gtdb_profile.py -i "$file" -o "${base_name}_gtdb.txt"
done

# -----------------------------
# 11. Merge GTDB profiles
merge_metaphlan_tables.py --gtdb_profiles "/path/to/gtdb/"*_gtdb.txt \
  > "/path/to/gtdb/gtdb_merged_rel_ab_table.txt"

# -----------------------------
# 12. Diversity analysis (R script)
# Using relative abundance GTDB tables
Rscript calculate_diversity.R \
  -f /path/to/gtdb/gtdb_merged_rel_ab_table.txt \
  --taxon_separator s__ -d beta -m aitchison

Rscript calculate_diversity.R \
  -f /path/to/gtdb/gtdb_merged_rel_ab_table.txt \
  --taxon_separator s__ -d alpha -m shannon
