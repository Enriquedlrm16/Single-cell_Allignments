# Single_cell_Allignments
Here I describe how to reach the same output cluster results, and the same number of cells with Alevinfry and CellRanger in Pbmc_1k dataset 

## Alevinfry

### Directory and folders content:
> ~/Downloads/single_cell_alevin$ ls
### pbmc_10k_v3_fastqs.tar  pbmc_1k_v3_fastqs.tar  refdata-gex-GRCh38-2020-A  simpleaf

> ~/Downloads/single_cell_alevin/simpleaf$ ls
### analyze_alignment.R  genes.gtf  genome.fa  index_dir  linux_commands_single_cell_alevinfry  output_dir  output_dir_  pbmc_1k_v3_fastqs  plist  r_script_simpleaf.R  simpleaf_info.json

### Install alevin-fry with conda:
> conda install -c bioconda alevin-fry

### Install specific version of salmon (=> 1.4.0):
> conda install salmon=1.4.0

### 
> conda create -n af -y -c bioconda -c conda-forge simpleaf piscem

### 
> conda activate af

### 
> tar xf - --strip-components=1 -C $FASTQ_DIR

### 
> export ALEVIN_FRY_HOME="$Downloads"

### 
> ulimit -n 4096

###
> gunzip -c pbmc_10k_v3_S1_L002_R2_001.fastq.gz | head | sed -n '2p' | wc -c

### 
> simpleaf index --output index_dir --fasta genome.fa --gtf genes.gtf --rlen 91 --threads 28 --use-piscem

### 
> simpleaf set-paths

### 
> simpleaf quant --reads1 pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L001_R1_001.fastq.gz,pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L002_R1_001.fastq.gz --reads2 pbmc_10k_v3_fastqs/  pbmc_10k_v3_S1_L001_R2_001.fastq.gz,pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L002_R2_001.fastq.gz --threads 20 --index index_dir/index --chemistry 10xv3 --resolution cr-like --unfiltered-pl --expected-ori fw --t2g-map index_dir/index/t2g_3col.tsv --output output_dir

## CellRanger

### 
> curl -o cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1708126137&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=kjBAMOEuKNyNXKXOgMuGF2n4-3XAaa7-hJvEnVB2P9jbXj14f6bG46FbuhUN06C14JmIJAicozUpEbYOst1VsdbmStedVU1b5IOKuaaKrVxWSLgzk4W~hvmBLEAyCuta~4iskl1Tehg6EVHbFzZ64I9q0gQpFeY8QI89p1nztjRoWL8MTuQjufVu1VSs2NEmlMIfHGdxL20Z18L8Fc1DIlRooM0CFu5URCwch3XZDvPtJh6ZyMWYYQAKEKXcfH9udWACdYNPH7-iiarliRQAU~w4ZSO-Yf0Fz4m68pqccN8JPylG9AaTRPegbqYFctaapk~0bvofYtKmNr9D5W8arA__"

###
> 
