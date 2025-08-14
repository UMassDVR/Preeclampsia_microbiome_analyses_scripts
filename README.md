# Preeclampsia Microbiome Scripts

This repository contains R and RMarkdown files used for microbiome analyses in a study of preeclampsia.  
The code is organized by body site and sequencing strategy (e.g., 16S rRNA and WGS).  

---

## Additional Contents

In addition to individual analysis scripts, this repository provides:

- **Ready-to-use phyloseq objects** containing processed microbiome data for direct analysis.
- **Full scripts** for generating phyloseq objects from raw or intermediate data, allowing complete reproducibility if regeneration is required.
- **Pipelines for sequence processing** from raw FASTQ files to fully processed phyloseq objects, ensuring transparency and reproducibility.

---

## Requirements

- **R** (≥ 4.0 recommended)

---

## Usage

1. Clone or download this repository.
2. Navigate to the desired directory and open the corresponding `.R` script or `.Rmd` notebook.
3. Review the header comments for input requirements and set any local file paths as needed.
4. Choose your starting point:
   - **Direct analysis**: Use the provided phyloseq objects.
   - **Full regeneration**: Run the object-generation scripts or sequence-processing pipelines starting from raw data.
5. When numeric prefixes are present in file names, run the scripts in that order.

---

## Data Availability

- **Raw sequencing data**: Deposited in a public repository such as the NCBI Sequence Read Archive (SRA) cited in the manuscript.
- **Processed phyloseq objects**: Provided in this repository for direct analysis.
- **Object-generation scripts and pipelines**: Included for full reproducibility if objects need to be regenerated from raw or intermediate data.

