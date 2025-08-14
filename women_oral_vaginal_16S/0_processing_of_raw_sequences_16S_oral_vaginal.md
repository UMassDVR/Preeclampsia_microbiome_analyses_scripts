# Processing of Raw 16S V1–V3 Sequences (Oral & Vaginal) — QIIME 2 + DADA2

# University of Massachusetts — Preeclampsia Microbiome Project
# 202 samples (200 biological + 2 controls) — 404 paired-end FASTQ files
# Target region: 16S rRNA V1–V3
# Pipeline: QIIME 2, DADA2, Greengenes2 Bayesian classifier

# 1. Connect to the HPC cluster
ssh user@hpc-server-address

# 2. (Optional) Run FastQC on one sample to check quality
mkdir fastqc_results
for file in *.fastq.gz; do
    fastqc "$file" -o fastqc_results/
done

# 3. Activate Conda environments
# Local
conda activate qiime2-amplicon-2024.5
# HPC
export PATH="/path/to/anaconda3/bin:$PATH"
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate qiime2-2023.5

# 4. (Optional) Move/copy all FASTQ files into a single folder
# for dir in /path/to/raw_fastq_directory/ds.*/
# do
#    mv "$dir"/*.fastq.gz /path/to/fastq_run_folder/ 2>/dev/null
# done

# 5. Import data into QIIME 2
bsub -q long -n 4 -R "rusage[mem=16000]" -W 4:00 -o qiime_import.log \
"qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /path/to/fastq_run_folder/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza"

# 6. Subsample (20%) due to large dataset size
bsub -q long -n 4 -R "rusage[mem=16000]" -W 16:00 -o qiime_subsample02.log \
qiime demux subsample-paired \
  --i-sequences demux-paired-end.qza \
  --p-fraction 0.2 \
  --o-subsampled-sequences demux-paired-end_subsample02.qza

# 7. DADA2 denoising
# Parameters: trim-left-f=10, trim-left-r=10, trunc-len-f=300, trunc-len-r=300
bsub -w "done(JOB_ID)" -q long -n 4 -R "rusage[mem=16000]" -W 60:00 -o qiime_dada_subsample02.log \
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end_subsample02.qza \
  --p-trim-left-f 10 \
  --p-trunc-len-f 300 \
  --p-trim-left-r 10 \
  --p-trunc-len-r 300 \
  --o-representative-sequences rep-seqs_subsample02.qza \
  --o-table table_subsample02.qza \
  --o-denoising-stats denoising-stats_subsample02.qza

# 8. Transfer results to local machine
scp user@hpc-server-address:/path/to/rep-seqs_subsample02.qza .
scp user@hpc-server-address:/path/to/table_subsample02.qza .
scp user@hpc-server-address:/path/to/denoising-stats_subsample02.qza .

# 9. Taxonomic classification (local)
qiime feature-classifier classify-sklearn \
  --i-classifier classifier_full_length_gg_ext.qza \
  --i-reads rep-seqs_subsample02.qza \
  --o-classification results/taxonomy_full_length_subsample02_gg_ext.qza

# 10. Export results
qiime tools export \
  --input-path results/taxonomy_full_length_subsample02_gg_ext.qza \
  --output-path taxonomy_full_length_subsample02_gg_ext
qiime tools export \
  --input-path rep-seqs_subsample02.qza \
  --output-path rep-seqs_subsample02

# 11. (Optional) Export FASTQ from a .qza file
mkdir -p exported_fastq
bsub -q short -n 1 -R "rusage[mem=8000]" -W 0:30 -o qiime_export.log \
"source /path/to/anaconda3/etc/profile.d/conda.sh; \
 conda activate qiime2-2023.5; \
 qiime tools export --input-path demux-paired-end.qza --output-path exported_fastq"
