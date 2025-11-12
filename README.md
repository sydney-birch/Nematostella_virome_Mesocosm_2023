# Nematostella_virome_Mesocosm_2023
This repository contains the virome, host transcriptome, and 16s rRNA work for our 2023 Mesocosm study examining how 1) different male and female Nematostella clonal lines and 2) knockout lines of immunity genes, interact with microbes and viruses from a natural, foreign environment (SC estuary). This mesocosom study was conducted in SC at Belle W. Baruch Marine Field Lab in Aug 2023. 

### Study Info and Sample Structure: 
- We have 94 total samples (12 strains, 2 time points [T0, T96], and 4 Reps for each sample [MAVs T96 and OAS T96 have 3 reps]):
-   Clonal lines: 
      - FL.Male  ; FL.Female  ;  NC.Male  ;  NC.Female  ;  ME.Male  ;  ME.Female
        
- KO Lines:
  - WT  ;  RLRa  ;  RLRb  ;  MAVs  ;  OAS  ;  RNaseL

Nematostella were sampled at our intial time point (T0) before going into the mesocosm and 96 hours (T96) after being exposed to natural estuary water from SC. The mesocosms were bins filled with water from the estuary that was bag filtered to remove larger debris. Water changes occured every other day, and each mesocosm bin had a pump to produce water flow for the animals. 

### Pipeline overview:  
#### Prep:
1) We extracted RNA and DNA using the Allprep RNA/DNA mini kit
2) We created libraries for total RNAseq using the Tecan Universal Plus total RNA-seq with NuQuant kit with an N. vectensis rRNA depletion step
3) We sequenced total RNA-seq samples using NovaseqX
4) Hypervariable regions V1-V2 of bacterial 16S rRNA genes were amplified with the extracted DNA from Step 1
5) We sequenced 16S samples using NextSeq
 
### Virome:
1) Prep work: Run Fastqc, trim reads, and run fastqc again
2) Map reads to the genome to get two pots of data --> mapped (Nematostella reads) and unmapped (virome/microbe reads)    
3) Map the unmapped reads back to the genome until 0 reads align
4) Align reads to rRNA database comprised of SILVA LSU & SSU, GTDB SSU, and Rfam --> unmapped reads presumably viral
5) Assemble viral transcriptome using Spades rnaviral - created:
    - A) 24 collapsed replicate assemblies,
    - B) 93 individual assemblies, and
    - C) Total T0 and Total T96 assemblies & calculated normalized read count **(USED in Sharoni et. al., 2025 - *in review*)**
6) Analyze Function    
   - A) Run MetaCerberus on assemblies    
   - B) Get functional terms: KEGG, COG, VOG, PHROG    
   - C) Visualize in R     
7) Analyze Taxonomy    
   - A) BLAST transcriptomes to RefSeq Viruses    
   - B) Filter BLAST output (>= 70% identity)    
   - C) NCBI Taxon Kit to get taxon IDs and lineage info    
   - D) visualize in R &  run diversity metrics
  
### Host: 
1) Use the mapped reads from Prep step 5 to align to the Nematostella genome using Salmon
2) Aanalyze Differentially expressed Genes (DEGs)
    - A) Use EdgeR to conduct pariwise comparisons (T0 vs T14)
    - B) Run MetaCerberus on DEGs
    - C) Get functional terms: KEG, FOAM
    - D) Visualize in R 
3) Analyze WGCNA
    - A) Run WGCNA in R
    - B) Run MetaCerberus on Modules from WGCNA
    - C) Get functional terms: KEGG, FOAM
    - D) Visualize in R
  
### Microbe (16s QIMME2 Analysis): 
1) Import reads into QIIME2
2) Sequence quality control and feature table construction using DADA2 plugin
3) Generate a tree for phylogenetic diversity analysis
4) Choose a sampling depth and run alpha and beta diversity analysis
5) Run taxonomic analysis
6) Import QIMME2 artifacts to R for visualization

  
