# 16S analysis workflow

We are picking up here on Step 9 of the overall analysis to analyze the microbiome using QIIME2. 

We have a total of 117 samples (12 strains, 2 time points [T0, T96]): ME.F-T0 has 3 reps, ME.M-T0 has 4 reps, all other samples have 5 reps.  
  - Clonal lines: 
    - FL.Male  ; FL.Female  ;  NC.Male  ;  NC.Female  ;  ME.Male  ;  ME.Female
        
  - KO Lines:
    - WT  ;  RLRa  ;  RLRb  ;  MAVs  ;  OAS  ;  RNaseL

## 9) 16s rRNA Analysis (QIIME2)   
I'm using the following tutorials in this analysis: 
  - Moving Pictures tutorial: https://docs.qiime2.org/2020.6/tutorials/moving-pictures/
  - Atacama soil tutorial (paired end): https://docs.qiime2.org/2024.10/tutorials/atacama-soils/
  - Gut to Soil axis tutorial: https://amplicon-docs.qiime2.org/en/stable/tutorials/gut-to-soil.html

### Prep work: 
   * Current data is in a shared projects folder - copy over the files into your working dir and change names   
`./9.0.A_change_fq_name_and_copy_over.py -a name_change.txt -b /projects/areitze2_research/03_10_25_16S_QK_SB/`   

   * Remove underscores in name
`./9.0.B_remove_underscores_in_name.py -b ../raw_reads/`   
   * Make a metadata file and copy to terminal
     `scp Meso_2023_ALL_16s_metadata.tsv sbirch1@hpc.charlotte.edu:../../projects/meso_2023/9_qiime2`
     
### Step 1: Import your data by creating a qiime artifact - follow Casava 1.8 paired-end demultiplexed fastq   

```
module load qiime2/amplicon-2024.10

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path raw_reads \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza
```
Visualize output (go to https://view.qiime2.org): 
```
qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux.qzv
```


### Step 2: Sequence quality control and feature table construction - I'm using option 1: DADA2

```
#part a - quality filter and truncate(denoise) 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 300 \
  --p-trunc-len-r 300 \
  --o-table table-dada2.qza \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-denoising-stats denoising-stats-dada2.qza

#part b - Make visualizations and summaries 
qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2.qzv \
  --m-sample-metadata-file metadata-file-2022_v2.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs-dada2.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats-dada2.qza \
  --o-visualization denoising-stats-dada2.qzv

 	#send to computer to visualize

	#rename files for next section  
		mv rep-seqs-dada2.qza rep-seqs.qza
		mv table-dada2.qza table.qza
```

### Step 3 Generate a tree for phylogenetic diversity analysis 
```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

Next chose a sampling depth based on table.qzv:
  - 61428 --> Retained 7,125,648 (36.12%) features in 116 (99.15%) samples at the specifed sampling depth. (4 reps OAS T96 and keeps all 3 reps of ME.F.T0 if go above at 2 reps)


*In Progress*
