# Host analysis workflow 

we are picking up here on Step 8 of the overall analysis to analyze the host transcriptome

## 8) Host (Nematostella) Analysis  
### Prep) Use Salmon to align reads and get counts

First, align the mapped (Nematostella) reads to the UK transcriptome to get counts and sequence IDs for the reads 

1.A) Copy over transcriptome: GCF_932526225.1/rna.fna (I renamed to uk_transcriptome.fa)

1.B) submit Salmon slurm: sbatch 8_salmon-mapped.slurm
`#create the index (only run 1 time)
salmon index -t uk_transcriptome.fa -i mapped_index -k 31`

`#run salmon script to align each replicate
./8_run_salmon.py -a ../4_RNA_filt/mapped_fastq_files -b ../../8_salmon -c aligned `

`#example line of code of salmon alignment line: 
salmon quant -p 12 --seqBias --gcBias -i mapped_index -l A -1 {0}/{1}_paired1.fastq.gz -2 {0}/{1}_paired2.fastq.gz -o {1}".format(fastq_dir,sample)`


### A) Differential Gene Expression     
1.A) copy over Salmon output files to computer for R analyses - place in a folder called: mapping

`scp -r sbirch1@hpc.charlotte.edu:./../../scratch/sbirch1/Nematostella_transcriptomics/8_salmon/aligned_NH_T* ./`

As a side quest - We first interogated the host gene expression of 56 previously identified Nematostella immune genes. This required: 
   * Text file of libraries_to_stages.tx (two column file annotating which reps belong to which groups)
   * Dataframe of Immune gene accid annotations: Immune_genes_56.csv
   * Rscript: Nematostella_Immune_heatmap_NEE.R

1.B) Run edgeR script in R to get Differentially Expressed Genes between Timepoints (T0 vs T14) for each location - write out into text files 
script: edgeR_meso_2022_MAPPED_14_groups.R

1.C) Input DEG lists back to the terminal and get header information for each DEG to be used with selectSeqs 

`#first run this script to turn the edgeR input into a list of headers to retrieve: 
./run_get_headers_from_edgeR.sh --> This runs run_get_headers_from_edgeR.py `

`#Next, move all header lists into a new dir called mapped_selectseqs_header_lists: 
mv *-header_list ../mapped_selectseqs_header_lists/`

`#Now use script to get full headers from initial accession id to run with selectSeqs: 
./3_get_full_headers.sh uk_transcriptome.fa --> this will run: script 3.B_get_full_headers.py`

1.D) Run selectSeqs to extract the sequences for each DEG from the transcriptome using the full headers:

`#Run select seq bash script: 
./run_selectseqs.sh --> this will run ./selectSeqs.pl`

1.E) Run DEGs on Metacerberus to get potential functional annotations 
`#Field: 
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/FIELD_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_FIELD_T0vT14_DEGs
#FL
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/FL_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_FL_T0_v_T14_DEGs        
#MA
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/MA_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_MA_T0_v_T14_DEGs  
#ME
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/ME_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_ME_T0_v_T14_DEGs        
#NH
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/NH_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_NH_T0_v_T14_DEGs
#NS
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/NS_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_NS_T0_v_T14_DEGs
#SC
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/SC_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_SC_T0_v_T14_DEGs
`

1.F) copy output to computer - turn counts into spreadsheets and make figures in R using script: Metacerberus_stats_output.R

1.G) Make DEG heatmaps 
`#run script to turn list of DEGs into a string to import to R: 
./prep_for_R_string_deg.py`

Make heatmaps in R using: Nematostella_mapped_DEG_heatmaps.R



### B) WGCNA    
 
*In Progress*
