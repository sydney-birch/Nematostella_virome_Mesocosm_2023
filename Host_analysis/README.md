# Host analysis workflow 

we are picking up here on Step 8 of the overall analysis to analyze the host transcriptome

## Step 8) Host (Nematostella) Analysis  

### 8.1) Use Salmon to align reads and get counts

  - First, align the mapped (Nematostella) reads to the UK transcriptome to get counts and sequence IDs for the reads 

  - 8.1.A) Copy over transcriptome: GCF_932526225.1/rna.fna (I renamed to uk_transcriptome.fa)

  - 8.1.B) submit Salmon slurm: sbatch 8.1.B_salmon-mapped.slurm
    - `#create the index (only run 1 time)
salmon index -t uk_transcriptome.fa -i mapped_index -k 31`

    - #run salmon script to align each replicate
      - `./8.1.B_run_salmon.py -a ../4_RNA_filt/mapped_fastq_files -b ../../8_salmon -c aligned `

    - example line of code of salmon alignment line: 
      - `salmon quant -p 12 --seqBias --gcBias -i mapped_index -l A -1 {0}/{1}_paired1.fastq.gz -2 {0}/{1}_paired2.fastq.gz -o {1}".format(fastq_dir,sample)`


### 8.2) Differential Gene Expression     
  - 8.2.A) copy over Salmon output files to computer for R analyses - place in a folder called: mapping
    - `scp -r sbirch1@hpc.charlotte.edu:./../../scratch/sbirch1/Nematostella_transcriptomics/8_salmon/aligned_NH_T* ./`

  - 8.2.B) Run edgeR script in R to get Differentially Expressed Genes between Timepoints (T0 vs T96) for each location/strain - write out into text files **IN PROGRESS**
    - script: edgeR_meso_2022_MAPPED_14_groups.R 

  - 8.2.C) Input DEG lists back to the terminal and get header information for each DEG to be used with selectSeqs **IN PROGRESS**

    - First run this script to turn the edgeR input into a list of headers to retrieve: 
      - `./run_get_headers_from_edgeR.sh --> This runs run_get_headers_from_edgeR.py `

    - Next, move all header lists into a new dir called mapped_selectseqs_header_lists: 
      - `mv *-header_list ../mapped_selectseqs_header_lists/`

    - Now use script to get full headers from initial accession id to run with selectSeqs: 
      - `./3_get_full_headers.sh uk_transcriptome.fa --> this will run: script 3.B_get_full_headers.py`

    - Run selectSeqs to extract the sequences for each DEG from the transcriptome using the full headers: 

      - Run select seq bash script: 
        - `./run_selectseqs.sh --> this will run ./selectSeqs.pl`

  - 8.2.D) Run DEGs on Metacerberus to get potential functional annotations **IN PROGRESS**
      - copy output to computer - turn counts into spreadsheets and make figures in R using script: Metacerberus_stats_output.R
    
  - 8.2.E) Examine Host Expression of Nematostella Immune gene set
      -  We interogated the host gene expression of 56 previously identified Nematostella immune genes. This required: 
         * Text file of libraries_to_stages.txt (two column file annotating which reps belong to which groups)
         * Dataframe of Immune gene accid annotations: Immune_genes_56.csv
         * Rscript: Nematostella_Immune_heatmap_NEE.R

  - 8.2.F) Make DEG heatmaps **IN PROGRESS**
    - run script to turn list of DEGs into a string to import to R: 
      - `./prep_for_R_string_deg.py`

    - Make heatmaps in R using: Nematostella_mapped_DEG_heatmaps.R

### 8.3) WGCNA    **IN PROGRESS**
 


