#  Workflow for Nematostella mesocosom 2023 analysis

## Prep work: 
A) Transfer over all data and adjust names by running the change_raw_fq_file_names.py with name_change.txt --> has the conversion of names from admera    

The admera fastq name Format: 21081FL-07-01-13_S62_L003_R2_001.fastq.gz    
Changed to: SC_T0-WT_B2_S62_L003_R2_001.fastq.gz  

`./0_change_raw_fq_file_names.py -a 0_name_change.txt -b raw_reads_10-7-24/`

   
B) Download genome 
 `wget https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_932526225.1/`

 
	
## 1) Run fastqc
Use 1_fastqc.py to loop thru and run fastqc for each sample --> run in raw reads dir (1_fastqc.slurm)  
  `./1_fastqc.py -a ../1_fastqc/before_trim`    
  
      The actual line of code for fastqc: 
      fastqc {fastq_file} -o ../1_fastqc/before_trim
*output in fastqc dir --> before_trim dir*

 
	
## 2) Trimm reads and re-run fastqc 
  A) Run trimmomatic to trim the adaptors (run in 2_trimmomatic dir) - the script loops thru and runs for each sample (2.A_trimmommatic.slurm)

`./2.A_trimmomatic.py -a ../raw_reads_10-7-24 -b ../2_trimmomatic`   

	The actual line of code for trimmomatic used: 
 
 	java -jar /apps/pkg/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 -threads 4 {path_to_R1} {path_to_R2} -baseout {Sample_name}_filtered.fq.gz ILLUMINACLIP:/apps/pkg/trimmomatic/0.39/adapters/TruSeq3-PE-2.fa:1:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36    
  
	Output - you get 4 files in this format in the 2_trimmomatic dir for each sample:
	    SC_T0-WT_B2_filtered_1P.fq.gz
	    SC_T0-WT_B2_filtered_1U.fq.gz
	    SC_T0-WT_B2_filtered_2P.fq.gz
	    SC_T0-WT_B2_filtered_2U.fq.gz
		#We are interested in the _1P and _2P files (paired files)


   B) Re-Run Fastqc in the fastqc dir (2.B_fastqc.slurm) 
`./2.B_fastqc.py -a ../1_fastqc/after_trim -b ../2_trimmomatic`
*output goes to fastqc dir - after_trim dir*  

 
	
## 3) Map reads to genome to get two pots of data --> mapped and unmapped    

   ### A) Get unmapped reads (non Nematostella reads aka viral/microbial reads)

	Map reads (HISAT2 and Bowtie) and use sam tools to get mapped (1st alignment) and unmapped reads
		To get unmapped reads (All non-nematostella reads), Re-align the unmapped multiple times (4x total) to the genome - run 2 alignments with HISAT2 then 2 alignments with bowtie using high sensitivity settings

  
	#1st alignment - HISAT2 
		#1st Ali - (run slurm in 1st_ali_2-20-24): 
			./3.A_initial_hisat2_v2.py -a ../../2_trimmomatic -b 1 -c 2nd_ali -d ../3_HISAT2/1st_ali_2-20-24/
                         
			 The actual hisat2 code run in script: 
			 hisat2 -q -p 12 -x ../Nematostella_genome -1 {path_to_R1}_filtered_1P.fq.gz -2 {path_to_R2}_filtered_2P.fq.gz -S Nematostella_genome_ali_{alignment_#}_{Sample_Name}.sam
				## Output: sam file with _ali_1 in a dir in next ali dir (2nd_ali dir) - it makes the second ali dir and all sub dirs which contain sam files (Nematostella_genome_ali_1_{1}.sam)		

  
		# run sam tools - make fastq files to be used in second ali of all reads that did NOT map to the genome (Run slurm in 2nd_ali dir)
			sbatch 3.B_samtools_processing.slurm
			./3.B_samtools_processing.py -b 1

                         The actual samtools code run in script:  
			 	samtools view Nematostella_genome_ali_{0}_*.sam -f 0x4 -h -b -o unalign_reads.bam
				samtools collate -u -O unalign_reads.bam | samtools fastq -1 unaligned_{0}_{1}_ali_paired1.fq -2 unaligned_{0}_{1}_ali_paired2.fq -0 /dev/null -s /dev/null -n
	
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_1_ali_paired1.fq)
						


	#2nd alignment - HISAT2 
		#2nd Ali - Run slurm in 2nd_ali dir (Have to hard code where index is - In 1st ali dir - adjust name if needed)
			sbatch 3.C_hisat2_2_ali.slurm
			./3.C_hisat2_next_ali_v2.py -b 2 -c 3rd_ali_BT 
		        	## Output: sam file with _ali_2 in a subdir in next ali dir (3rd_ali dir) - it makes the third ali dir and all sub dirs which contain sam files (Nematostella_genome_ali_2_{1}.sam)

	    
		# run sam tools - make fastq files to be used in third ali (Run slurm in 3rd_ali dir) 
			sbatch 3.B_samtools_processing.slurm
			./3.B_samtools_processing.py -b 2
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_2_ali_paired1.fq)     	
		        	
		        	
		        	
	#3rd alignment - bowtie
		#First - copy over genome file ending in .fna (in 1st_ali dir) - need to make a bowtie index
			
		#3rd Ali - Run slurm in 3rd_ali dir (Have to hard code where index is - In 1st ali dir - adjust name if needed)
			sbatch 3.D_bowtie_ali_3.slurm
			./3.D_bowtie_ali_initial.py -b 3 -c 4th_ali_BT

      			The actual code run in script: 
	 			bowtie2 -q -p 12 --very-sensitive-local -x ../Nematostella_genome -1 {sample_name}_paired1.fq -2 {sample_name}_paired2.fq -S Nematostella_genome_BT_ali_{0}_{1}.sam".format
    
				# Output: sam file with BT_ali_3 in a dir in next ali dir (4th_ali dir) - it makes the fourth ali dir and all sub dirs which contain sam files (Nematostella_genome_BT_ali_3_{1}.sam)

   
		# run sam tools - make fastq files to be used in fourth ali (Run slurm in 4th_ali dir)	
			sbatch 3.E_samtools_processing_BT.slurm
			./3.E_samtools_processing_BT.py -b 3
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_3_ali_paired1.fq)
		
	
	#4th alignment - bowtie (hopefully final alignment - if it is rename 5th_ali to final_sam_files)
		#4th Ali - Run slurm in 4th_ali dir (Have to hard code where index is - In 1st ali dir - adjust name if needed)
			sbatch 3.F_bowtie_4_ali.slurm
			./3.F_bowtie_next_ali.py -b 4 -c 5th_ali
				## Output: sam file with BT_ali_4 in a dir in next ali dir (5th_ali dir) - it makes the fifth ali dir and all sub dirs which contain sam files (Nematostella_genome_BT_ali_4_{1}.sam)
    
			#### rename 5th_ali to final_unmapped_sam_files #### 
							
					
		# run sam tools - make fastq files  (Run slurm in 5th_ali dir/ final_sam_files dir)	
			sbatch 3.E_samtools_processing_BT.slurm
			./3.E_samtools_processing_BT.py -b 4
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_4_ali_paired1.fq)   


	#After each alignment run the count reads script and add counts to excel sheet
		#the script will iterate through each sample directory and count the number of reads in the fastq files
			sbatch 3.B.2_count_reads_after_samtools.slurm
			./3.B.2_count_reads_after_samtools.py
				


 ### B) Get mapped reads (Nematostella reads) and move unmapped reads into next step dir (RNA filtration)    
First, make 3 dirs:    
     `mkdir /Nematostella_transcriptomics/4_RNA_filt`    
     `mkdir /Nematostella_transcriptomics/4_RNA_filt/mapped_fastq_files`     
     `mkdir /Nematostella_transcriptomics/4_RNA_filt/unmapped_fastq_files`


Next, Run samtools to make fastq files for each sample of the MAPPED reads to the genome from the first alignment (in 2nd_ali dir)     
        `sbatch 3.H_samtools_mapped_reads.slurm`    
	`./3.H_samtools_mapped_reads.py -b 1 `    
	    `#Output: 2 fastq files (aligned_{sample_name}_1_ali_paired1.fastq) placed in 4_RNA_filt/mapped_fastq_files`


Move Unmapped reads into RNA filtration main dir        	  
        `sbatch 3.G_mv_fastqs.slurm`      
	`./3.G_mv_fastqs_v2.py -b 4 -c 4_RNA_filt `     
		`#Output: all unmapped fastqs will be moved to the RNA_filt unmapped_fastq_files dir`

  ## 4) rRNA Filtration of Viral/Microbial Reads
  Ultimately, in this step we want to have presumably all viral reads, so we need to filter out any prokaryotic and eukaryotic reads.    

  
  A) Download files to create your rRNA database     
  `./4.A_RNAdb_download_files.sh`    

  
  B) Concatenate all the downloads and run cd-hit to remove duplicate sequences.     
  `sbatch 4.B_cat_cd-hit.slurm` this runs: `./4.B_cat_cd-hit.sh RNA_db rRNAdb`     

  
  C) Map reads to rRNA db using hisat2       
  `sbatch 4.C_filtered_hisat2.slurm` this runs: `./4.C_filtered_hisat2.py -a filtered_unmapped_fastq_files -b rRNAdb_ed1.fa -c unmapped_fastq_files`     

  
  D) Use samtools to make fastq files of the reads that did *not* map to rRNA db - these are the primarily viral reads      
  `sbatch 4.D_samtools_processing.slurm` this runs: `./4.D_samtools_processing.py -a filtered_unmapped_fastq_files`      

  
  E) As an additional filtration step, we also mapped the reads to the NCBI nr database with viruses removed     
  `sbatch 4.E_nr_filtered_hisat2.slurm` this runs: `./4.E_nr_filtered_hisat2.py -a ncbi_nr_filtered_unmapped_fastq_files -b ncbi_nr_DB.fa -c filtered_unmapped_fastq_files`      

  
  F) Use samtools to make fastq files of the reads that did *not* map to nr database       
  `sbatch 4.F_samtools_processing.slurm ` this runs: `./4.F_samtools_processing.py -a ncbi_nr_filtered_unmapped_fastq_files`     
  
  
 ## 5) Make Viral Assemblies 

 Here, we are going to generate three types of assemblies: A) Collapsed assemblies - where we collapse the replicates to make an assembly for each sample ; B) Individual assemblies - where we will assemble each individual replicate  ; and C) Total assemblies - where we will combine all T0 samples into a T0 assembly and same for T96. The Total assemblies will be used to get normalized read counts. 

 ### A) Collapsed replicate Assemblies 
assemble replicates for each sample (24 total)
 `sbatch 5.A.2_run_spades_collapse_reps.slurm` this runs: `./5.A.2_run_spades_collapse_reps.py -i ../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/` 
 
 Use spades --rnaviral to generate the viral assemblies, an example line of code:     
 ```spades.py --rnaviral --pe1-1 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B1_4_ali_paired1.fq.gz --pe2-1 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B2_4_ali_paired1.fq.gz --pe3-1 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B4_4_ali_paired1.fq.gz --pe1-2 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B1_4_ali_paired2.fq.gz --pe2-2 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B2_4_ali_paired2.fq.gz --pe3-2 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B4_4_ali_paired2.fq.gz -o NS_T0_output```     
    

### B) Individual replicate Assemblies
assemble each individual replicate (93 total)
 `sbatch 5.A_run_spades.slurm` this runs: `./5.A_run_spades.py -c ../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/` 

 Use spades --rnaviral to generate the viral assemblies, an example line of code:     
 ```spades.py --rnaviral -1 {0}_paired1.fq.gz -2 {0}_paired2.fq.gz -o ../../5_assemblies/{1}_spades_output```     
    
### C) Total Assemblies 
assemble all T0 samples into a T0 assembly and all T96 samples into a T96 assembly (2 total)

##### C.1) Make two dirs SC_T0_fastqs and SC_T96_fastqs
Copy over all fastq files to respective directory

##### C.2) Concatenate all T0 and all T96 fastq files (run on both)
`sbatch A.3_make_totals_files.slurm` runs: `./5.A.3_make_totals_files.py -a ../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/SC_T0_fastqs -b filtered -c ../../../5_assemblies/ -d SC_T0`

##### C.3) Run spades to create two total assemblies
`sbatch ./5.A.4_run_spades_totals.slurm` runs: 

	spades.py --rnaviral -1 SC_T0-total_R1.fq.gz -2 SC_T0-total_R2.fq.gz -o SC_T0_Total_spades_output
	spades.py --rnaviral -1 SC_T96-total_R1.fq.gz -2 SC_T96-total_R2.fq.gz -o SC_T96_Total_spades_output`

##### C.4) Using the two Total Assemblies - Align fastq files to get normalized read counts
Next, align the rRNA filtered unmapped reads (presumably viral reads) to the T0 and T96 total assemblies respectively using HISAT2

    #Align T0 filtered reads to T0 total assembly: 
	 ./5.A.5_hisat2_total_assemblies.py -a T0_total_mapping_fastqs -b total-SC-T0_assembly.fa_FIX.fa -c ../../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/ -d T0
	 
  	#Align the T96 filtered reads to the T96 total assembly:
    ./5.A.5_hisat2_total_assemblies.py -a T96_total_mapping_fastqs -b total-SC-T96_assembly.fa_FIX.fa -c ../../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/ -d T96

	#use samtools to create fastq files of the mapped reads:
	sbatch 5.A.6_samtools_mapped.slurm
    ./5.A.6_samtools_mapped.py -a T0_total_mapping_fastqs 
    ./5.A.6_samtools_mapped.py -a T96_total_mapping_fastqs

	#get counts of mapped reads
     ./4.G_count_reads_in_fastqs.py -a T0_total_mapping_fastqs
     ./4.G_count_reads_in_fastqs.py -a T96_total_mapping_fastqs

Next, calculate the Normalized Read Count in excel:     

*Number of mapped viral reads to the total assemblies / Number of total raw reads (number of reads after sequencing) * 1000*


 
 ## 6) Viral Taxonomy Analysis     



 ## 7) Viral Functional Analysis       


 
## 8) Host (Nematostella) Analysis    

### A) Differential Gene Expression     


### B) WGCNA    
 
  
## 9) 16s rRNA Analysis (QIIME2)   
	
