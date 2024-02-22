# Single_cell_Allignments
Here I describe how to reach the same output cluster results, and the same number of cells with Alevinfry and CellRanger in Pbmc_1k dataset 

## Alevinfry

### Directory and folders content:
> ~/Downloads/single_cell_alevin$ ls
#### pbmc_10k_v3_fastqs.tar  pbmc_1k_v3_fastqs.tar  refdata-gex-GRCh38-2020-A  simpleaf

> ~/Downloads/single_cell_alevin/simpleaf$ ls
#### analyze_alignment.R  genes.gtf  genome.fa  index_dir  linux_commands_single_cell_alevinfry  output_dir  output_dir_  pbmc_1k_v3_fastqs  plist  r_script_simpleaf.R  simpleaf_info.json

### Install alevin-fry with conda:
> conda install -c bioconda alevin-fry

### Install specific version of salmon (=> 1.4.0):
> conda install salmon=1.4.0

### conda create: This is the command to create a new conda environment.
### -n af: This specifies the name of the new environment, in this case, "af".
### -y: This flag indicates that you want to automatically confirm the installation without asking for user confirmation.
### -c bioconda: This specifies a channel (a repository of packages) from which to install packages. In this case, it's the "bioconda" channel, which contains bioinformatics-related packages.
### -c conda-forge: This specifies another channel, "conda-forge," which is a community-driven collection of conda packages.
### simpleaf: This is the name of a package to be installed in the new environment. It appears to be a package related to bioinformatics.
### piscem: This is another package to be installed in the new environment, likely also related to bioinformatics.
> conda create -n af -y -c bioconda -c conda-forge simpleaf piscem

### activate the conda environment previously named as "af"
> conda activate af

### --strip-components=1: This option tells tar to strip the leading directory component from file names. For example, if the tar archive contains files like directory/file1, directory/file2, etc., using --strip-components=1 would extract them as file1, file2, etc., removing the leading directory/ part.
### -C $FASTQ_DIR: This option specifies the directory where the extracted files should be placed. $FASTQ_DIR is likely a variable that holds the directory path where you want to extract the files.
> tar xf - --strip-components=1 -C $FASTQ_DIR

### Sets an environment variable named ALEVIN_FRY_HOME to the value of the $Downloads variable.
> export ALEVIN_FRY_HOME="$Downloads"

### Sets the maximum number of open file descriptors (file handles) for the current shell session to 4096.
> ulimit -n 4096

### Extracts the second line from the decompressed FASTQ file pbmc_10k_v3_S1_L002_R2_001.fastq.gz and counts the number of bytes in that line. It's essentially getting the length of the second line of the uncompressed FASTQ file. Just to determine the rlen parameter in the next function as: this output minus 1.
> gunzip -c pbmc_10k_v3_S1_L002_R2_001.fastq.gz | head | sed -n '2p' | wc -c

### simpleaf index: This is the main command for indexing using the simpleaf tool.
### --output index_dir: This specifies the directory where the index files will be saved.
### --fasta genome.fa: This specifies the input genome FASTA file (genome.fa) to be used for indexing.
### --gtf genes.gtf: This specifies the gene annotation file (genes.gtf) to be used for indexing.
### --rlen 91: This specifies the read length (rlen) as 91 base pairs. It's likely used for indexing calculations.
### --threads 28: This specifies the number of threads or CPU cores to use for the indexing process.
### --use-piscem: This flag indicates to use Piscem, which is likely another tool or method integrated with simpleaf for certain indexing tasks.
> simpleaf index --output index_dir --fasta genome.fa --gtf genes.gtf --rlen 91 --threads 28 --use-piscem

### Sets or configures paths related to the simpleaf tool. 
> simpleaf set-paths

### Quantification on paired-end RNA sequencing data using simpleaf, with specific settings for read files, index files, library chemistry, resolution, and output directory. 
### --threads 20: This specifies the number of threads or CPU cores to use for the quantification process.
### --resolution cr-like: This specifies the resolution of the quantification. The term "cr-like" might refer to a specific resolution or normalization method used in the analysis.
### --unfiltered-pl: This flag indicates to output unfiltered pseudo-likelihood (PL) values.
### --expected-ori fw: This specifies the expected orientation of the reads. In this case, it's "fw", which likely stands for "forward".
### --t2g-map index_dir/index/t2g_3col.tsv: This specifies the transcript-to-gene mapping file used for quantification.
> simpleaf quant --reads1 pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L001_R1_001.fastq.gz,pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L002_R1_001.fastq.gz --reads2 pbmc_10k_v3_fastqs/  pbmc_10k_v3_S1_L001_R2_001.fastq.gz,pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L002_R2_001.fastq.gz --threads 20 --index index_dir/index --chemistry 10xv3 --resolution cr-like --unfiltered-pl --expected-ori fw --t2g-map index_dir/index/t2g_3col.tsv --output output_dir

### In R, we have to do the same filtering as CellRanger, that is more or less "cells_with_800_umi <- cell_sums >= 800": 
> a <- fishpond::loadFry("/home/jrivas/Downloads/single_cell_alevin/simpleaf/output_dir/af_quant/",
                         outputFormat = "scrna", nonzero = T)

> cell_sums <- colSums(counts(a))

> cells_with_800_umi <- cell_sums >= 800

> a_filtered <- a[,cells_with_800_umi]

> seurat_obj <- Seurat::CreateSeuratObject(counts(a_filtered), )

### To continue the scRNA-seq pipeline in R, see my own scRNA-seq workflow: https://github.com/Enriquedlrm16/scRNA-seq_Workflow

## CellRanger

### Directory and folders content:
> ~/Downloads/cellranger_probe_sc$ ls
#### cellranger-7.2.0  cellranger-7.2.0.tar.gz  cellranger-7.2.0.tar.tar  pbmc_1k_v3_fastqs  pbmc_1k_v3_fastqs.tar  refdata-gex-GRCh38-2020-A  refdata-gex-GRCh38-2020-A.tar.gz

### Fetches the Cell Ranger software version 7.2.0 tarball from the specified URL and saves it as cellranger-7.2.0.tar.gz in the current directory.
> curl -o cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1708126137&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=kjBAMOEuKNyNXKXOgMuGF2n4-3XAaa7-hJvEnVB2P9jbXj14f6bG46FbuhUN06C14JmIJAicozUpEbYOst1VsdbmStedVU1b5IOKuaaKrVxWSLgzk4W~hvmBLEAyCuta~4iskl1Tehg6EVHbFzZ64I9q0gQpFeY8QI89p1nztjRoWL8MTuQjufVu1VSs2NEmlMIfHGdxL20Z18L8Fc1DIlRooM0CFu5URCwch3XZDvPtJh6ZyMWYYQAKEKXcfH9udWACdYNPH7-iiarliRQAU~w4ZSO-Yf0Fz4m68pqccN8JPylG9AaTRPegbqYFctaapk~0bvofYtKmNr9D5W8arA__"

### This command downloads the specified tarball and saves it to the current directory with the same name as the remote file (refdata-gex-GRCh38-2020-A.tar.gz)
> curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"

### Extract the contents of the cellranger-7.2.0.tar.gz file, printing out the names of the files as they are extracted
> tar -xzvf cellranger-7.2.0.tar.gz

### Extract the contents of the refdata-gex-GRCh38-2020-A.tar.gz file, printing out the names of the files as they are extracted
> tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

### Modifies the PATH environment variable to include the directory where the cellranger-7.2.0 executable is located.
> export PATH=/home/jrivas/Downloads/cellranger_probe_sc/cellranger-7.2.0:$PATH

### Sets the maximum number of open file descriptors (file handles) for the current shell session to 1024.
> ulimit -n 1024

### Run the cellranger count pipeline on the specified input FASTQ files, using the specified transcriptome reference data. It will produce an output directory named run_count_1kpbmcs containing the results of the gene expression analysis.
> /home/jrivas/Downloads/cellranger_probe_sc/cellranger-7.2.0/bin/cellranger count --id=run_count_1kpbmcs --fastqs=/home/jrivas/Downloads/cellranger_probe_sc/pbmc_1k_v3_fastqs --sample=pbmc_1k_v3 --transcriptome=/home/jrivas/Downloads/cellranger_probe_sc/refdata-gex-GRCh38-2020-A

### In R
>  b <- Read10X("/home/jrivas/Downloads/cellranger_probe_sc/run_count_1kpbmcs/outs/filtered_feature_bc_matrix/")

### To continue the scRNA-seq pipeline in R, see my own scRNA-seq workflow: https://github.com/Enriquedlrm16/scRNA-seq_Workflow
