# Pre-xQTL scan RNA-seq sample genotype processing

This is for pulling the genotypes of sample subsets, identifying unrelated samples, and performing genotype PCA. These results are needed when selecting samples for xQTL scans, and may be used as scan covariates.

## Input

This code requires the following:
* Paths to full TOPMed BCF files
* A list of NWD IDs to analyze
* Ancestry estimates for those NWD IDs
* A list of SNPs to use for genotype PCA

## Output

* Subsetted BCF and plink files
* Unrelated sample set
* Genotype PCs

## Running

You must only have NextFlow (>= v. 21.04) and Singularity (v. 3) installed. NextFlow should be configured as appropriate for your computing platform.

### 1. Generate subsetted VCF and plink files

```bash
nextflow run -resume --results /path/to/pipeline_results \
                        --samples samples.txt \
                        --bcf_glob '/path/to/topmed/chr*.bcf' \
                        fetch-genotypes.nf
```

where `samples.txt` is a list of sample IDs to be pulled from the TOPMed BCF files, and `--bcf_glob` is a shell glob capturing the TOPMed BCF files (each named: {chrom}.bcf). This should be in quotes as the glob itself should be passed to the pipeline.

### 2. Generate a set of unrelated samples

```bash
nextflow run -resume --results /path/to/pipeline_results \
                        --vcf_glob '/path/to/subsetted/*.vcf' \
                        get-unrelated-samples.nf
```

where `--vcf_glob` is a shell glob capturing the subsetted VCF files from step one (each named: {chrom}.vcf). This should be in quotes as the glob itself should be passed to the pipeline.

### 3. Calculate genotype PCs

```bash
nextflow run -resume --results /path/to/pipeline_results \
                        --samples samples.txt \
                        --variants variants.txt \
                        --ancestries ancestries.txt \
                        --vcf_glob '/path/to/subsetted/*.vcf' \
                        --vcf_index_glob '/path/to/subsetted/*.vcf.tbi' \
                        genotype-pca-using-predefined-snps.nf
```

where `samples.txt` is the unrelated sample set output from step 2; `variants.txt` is a list of variants to use for the genotype PCA (not all will necessarily be used); `ancestries.txt` is a file denoting the per-ancestry fraction for each sample (used for plotting only; file has a header; first column is NWD ID, all remaining columns correspond to an ancestry, and each value is between 0 and 1, representing the estimated fraction of that ancestry in that sample); and `--vcf_glob` and `vcf_index_glob` are shell globs capturing the subsetted VCF files from step one (each named: {chrom}.vcf and {chrom}.vcf.tbi)