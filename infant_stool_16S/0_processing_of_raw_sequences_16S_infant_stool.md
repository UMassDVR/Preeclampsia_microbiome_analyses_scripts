# Processing of Raw 16S V3–V4 Sequences (Infant Stool / Meconium) — QIIME 2 + DADA2
# University of Massachusetts — Preeclampsia Project
# 79 infant meconium samples + 1 negative control
# Target: 16S rRNA V3–V4 region
# Pipeline: QIIME 2, DADA2, Greengenes2 Bayesian classifier

# -----------------------------
# 1. Connect to the HPC
ssh user@hpc-server-address

# -----------------------------
# 2. Create working directory and move into sequencing folder
mkdir 2025-04-14_DCH_16S_NextSeq
cd /path/to/sequencing_amplicon/2025-04-14_DCH_16S_NextSeq

# -----------------------------
# 3. Transfer raw FASTQ files from Illumina BaseSpace to HPC
# (Use a BaseSpace CLI transfer script, e.g., transfer_files_illumina_to_UMass.sh)

# -----------------------------
# 4. (Optional) Run quality check on HPC
# fastqc and multiqc scripts adapted for UMass HPC

# -----------------------------
# 5. Count number of sample folders
ls ds.* | wc -l

# -----------------------------
# 6. Activate QIIME 2 environment on HPC
export PATH="/path/to/anaconda3/bin:$PATH"
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate qiime2-2023.5

# -----------------------------
# 7. Copy FASTQ files into working space
mkdir -p /path/to/16S_meconium_2025/seqs
find /path/to/sequencing_amplicon/2025-04-14_DCH_16S_NextSeq/ \
  -type f -name "*.fastq.gz" -path "*/ds.*/*" \
  -exec cp {} /path/to/16S_meconium_2025/seqs/ \;

# -----------------------------
# 8. Remove undetermined reads
rm /path/to/16S_meconium_2025/seqs/Undetermined_S0_L001_R*

# -----------------------------
# 9. Import data into QIIME 2
bsub -q long -n 4 -R "rusage[mem=16000]" -W 4:00 -o qiime_import.log \
"qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /path/to/16S_meconium_2025/seqs/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza"

# -----------------------------
# 10. Summarize demultiplexed data
export OPENBLAS_NUM_THREADS=1
qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux-paired-end.qzv

# -----------------------------
# 11. DADA2 denoising
# Parameters: trim-left-f=17, trunc-len-f=290, trim-left-r=21, trunc-len-r=250
bsub -q long -n 4 -R "rusage[mem=16000]" -W 60:00 -o qiime_dada.log \
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trim-left-f 17 \
  --p-trunc-len-f 290 \
  --p-trim-left-r 21 \
  --p-trunc-len-r 250 \
  --o-representative-sequences rep_seqs_demux.qza \
  --o-table table_demux.qza \
  --o-denoising-stats denoising_stats_demux.qza

# -----------------------------
# 12. Transfer results to local machine
scp user@hpc-server-address:/path/to/16S_meconium_2025/rep_seqs_demux.qza .
scp user@hpc-server-address:/path/to/16S_meconium_2025/table_demux.qza .
scp user@hpc-server-address:/path/to/16S_meconium_2025/denoising_stats_demux.qza .
scp user@hpc-server-address:/path/to/16S_meconium_2025/demux-paired-end.qza .

# -----------------------------
# 13. Activate QIIME 2 on local machine
conda activate qiime2-amplicon-2024.10

# -----------------------------
# 14. Taxonomic classification (Greengenes2, full-length)
qiime feature-classifier classify-sklearn \
  --i-classifier gg2_full16S_trained_classifier/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza \
  --i-reads rep_seqs_demux.qza \
  --o-classification taxonomy/taxonomy.qza

# -----------------------------
# 15. Export results
qiime tools export \
  --input-path taxonomy/taxonomy.qza \
  --output-path taxonomy

qiime tools export \
  --input-path rep_seqs_demux.qza \
  --output-path rep_seqs_demux
