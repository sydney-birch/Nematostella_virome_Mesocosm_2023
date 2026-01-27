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
To investigate the viral taxa found in our analysis, we ultimately used NCBI taxon kit. We first had to identify viral protiens and find their taxon IDs so we used BLAST and retrieved taxids using the Refseq Catalog (much more efficent than using bioentrez like before).   
We conducted this analysis on both the individual assemblies and the collapsed assemblies      


### 1) Run BLAST on Individual Assemblies 

1.A) Download/copy over the refseq viral database.    
`wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz `     
*This has 683,238 viral protiens*      


1.B) Make blast databases for each of the 93 viral assemblies     
`./6.0_blastdb.sh`     


1.C) Run a BLASTp using the RefSeq viral database against each of the 93 viral asssemblies.
`sbatch 6.1_blast.slurm`  runs --> `./1_blast.sh viral.1.protein.faa`   

This runs a blast search --> `blastp -query $1 -db $fasta -out ${fasta%.}_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 0.00001 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1`        

Get a count of the number of hits:   
`1.A_get_blast_hit_counts.py -a blastout`   


1.D) Decided to trim the blastout to 50 percent identity - more managable and more stringent (originally ran at 70, 60, 50, 40, 30, and 10 percent id - chose 50 percent id):    
```
			# run 50%: 
				./1.B_run_trim_blastout_50pid.sh
					#runs: 1.B_trim_blastout_50_pident.py
				
				# Get counts: 
					./1.A_get_blast_hit_counts.py -a 50_pi_blastout
```


1.E) Get viral accession IDs for the 50_pi_blastout table from ncbi file (check this later)
`./2_get_accessions.sh`


1.F) Get the full headers from the accession IDs to run with selectSeqs   
```
./3_get_full_headers.sh viral.1.1.genomic.fna
		#this will run script 3.B_get_full_headers.py

#Submit slurm 
sbatch 6.B_get_full_headers.slurm

#copy blastout 50% accids and full headers from ternimal to computer - input into spreadsheet
```   

1.G) Get Taxids using the Refseq Catalog
   * Access Refseq catalog: https://ftp.ncbi.nlm.nih.gov/refseq/release/README
      * download: `wget https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release230.catalog.gz`
      * unzip: `gunzip RefSeq-release230.catalog.gz`
      * will search accid and return taxonomy id

   * First split the T0s and T96s, and KO and Clonal into 4 dirs
     ```
     mkdir T0_Clonal_50_pi_hit1_ALL_accessions-taxids
	 mkdir T96_Clonal_50_pi_hit1_ALL_accessions-taxids
			
	 mkdir T0_KOs_50_pi_hit1_ALL_accessions-taxids
	 mkdir T96_KOs_50_pi_hit1_ALL_accessions-taxids
     ```
     * Copy over accession IDs
       ```
        cp 50_pi_hit1_accessionss/SC_T0*_hits T0_KOs_50_pi_hit1_ALL_accessions-taxids
		cp 50_pi_hit1_accessions/SC_T96*_hits T96_KOs_50_pi_hit1_ALL_accessions-taxids
			
		mv T0_KOs_50_pi_hit1_ALL_accessions-taxids/SC_T0*.F T0_Clonal_50_pi_hit1_ALL_accessions-taxids
		mv T0_KOs_50_pi_hit1_ALL_accessions-taxids/SC_T0*.M T0_Clonal_50_pi_hit1_ALL_accessions-taxids
			
		mv T96_KOs_50_pi_hit1_ALL_accessions-taxids/SC_T96*.F T96_Clonal_50_pi_hit1_ALL_accessions-taxids
		mv T96_KOs_50_pi_hit1_ALL_accessions-taxids/SC_T96*.M T96_Clonal_50_pi_hit1_ALL_accessions-taxids
	   ```

   * Submit slurm script to run get_taxids.sh from the RefSeq catalog
     ```
			#T0 Clonal
			#Adjust the bash script for T0
			sbatch 4_get_taxids.slurm
			./4_get_taxids.sh RefSeq-release230.catalog	
			
			#T96 Clonal
			#Adjust the bash script for T96
			sbatch 4_get_taxids.slurm
			./4_get_taxids.sh RefSeq-release230.catalog	

			#T0 KO
			#Adjust the bash script for T0
			sbatch 4_get_taxids.slurm
			./4_get_taxids.sh RefSeq-release230.catalog	
			
			#T96 KO
			#Adjust the bash script for T96
			sbatch 4_get_taxids.slurm
			./4_get_taxids.sh RefSeq-release230.catalog	

     #copy over taxids_accids to computer

     #Next In chatgpt make excel files of the output files: two seperate ones KOs and Clonals for T0 and T96 seperately (combine after):
	
			Ask: can you put these files into an excel sheet where each file is a tab titled the name of the file?
				 can you put these files into an excel sheet where each file is a tab titled the name of the file? Can you not add header formatting
			Attach the 8 files (10 max) - combine in to a totals spreadsheet
				Ask: can you combine these three workbooks with no alterations to the tab names
			check counts after
     ```

1.H) Run Taxon kit to get lineage information    
   * Copy over taxon_link database info from Meso_22 dir   
   * Copy over taxid files into final_taxids dir   

Run taxon kit: 
`sbatch 6.5_run_taxonkit.slurm` this runs --> `./5_run_taxonkit.sh final_taxids`   

example of code:   
`cat final_taxids/FIELD-T0_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> FIELD-T0_linage.txt`   

output will be a dir lineage_files that has all info in it --> copy to computer


1.I) Analyze Taxonomy data in R 
   * First you'll need to process the data in excel - make a total_Lineage files (two separate KO and Clonal) - import each lineage file and adjust data to columns
   	    * Each tab in these file is a location/strain
        * Have chatGPT combine the spreadsheets into KO and Clonal spreadsheets
   * Additionally, make a spreadsheets called (clonal or KO)_genus.csv and copy each location into a column in this sheet - this will be used as the input in R script
   * Make a third spreadsheet of the genetic composition of the viruses present using the VMR spreadsheet using xlookup
        * Run two xlookups - one using the species column, the other using virus name column to conduct your search
        * Then run an if statement to merge the data: `=IF(V2=0,W2,V2)` where V2 = species col and W2 = virus name col
        * remove s__ by find and replace
        * Then do a find and replace of 0 (entire cell) with Unknown (entire workbook)


Now run R scripts to look at taxanomic overlaps and run diversity statistics: Individual_genus_work.R


### 2) Run BLAST on Collapsed Assemblies 

2.A) Download/copy over the refseq viral database.    
`wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz `      

2.B) Make blast databases for each of the 24 collapsed viral assemblies     
`./6.0_blastdb.sh`     


2.C) Run a BLASTp using the RefSeq viral database against 24 collapsed viral asssemblies.
`sbatch 6.1_blast.slurm`  runs --> `./1_blast.sh viral.1.protein.faa`   

This runs a blast search --> `blastp -query $1 -db $fasta -out ${fasta%.}_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 0.00001 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1`        

Get a count of the number of hits:   
`1.A_get_blast_hit_counts.py -a blastout`   


2.D) Decided to trim the blastout to 50 percent identity - more managable and more stringent (originally ran at 70, 60, 50, 40, 30, and 10 percent id - chose 50 percent id):    
```
			# run 50%: 
				./1.B_run_trim_blastout_50pid.sh
					#runs: 1.B_trim_blastout_50_pident.py
				
				# Get counts: 
					./1.A_get_blast_hit_counts.py -a 50_pi_blastout
```


2.E) Get viral accession IDs for the 50_pi_blastout table from ncbi file (check this later)
`./2_get_accessions.sh`


2.F) Get the full headers from the accession IDs to run with selectSeqs   
```
./3_get_full_headers.sh viral.1.1.genomic.fna
		#this will run script 3.B_get_full_headers.py

#Submit slurm 
sbatch 6.B_get_full_headers.slurm

#copy blastout 50% accids and full headers from ternimal to computer - input into spreadsheet
```   

2.G) Get Taxids using the Refseq Catalog

   * First split the T0s and T96s, into 2 dirs
     ```
     mkdir T0_ALL_accessions_50_pi
	 mkdir T96_ALL_accessions_50_pi
     ```
     * Copy over accession IDs
       ```
        cp 50_pi_hit1_accessions/SC_T0* T0_ALL_accessions_50_pi/
		cp 50_pi_hit1_accessions/SC_T96* T96_ALL_accessions_50_pi/
	   ```

   * Submit slurm script to run get_taxids.sh from the RefSeq catalog
     ```
			#T0 
			#Adjust the bash script for T0
			sbatch 6.4_get_taxids.slurm 
			./4_get_taxids.sh RefSeq-release230.catalog	
			
			#T96 
			#Adjust the bash script for T96
			sbatch 6.4_get_taxids.slurm 
			./4_get_taxids.sh RefSeq-release230.catalog	

     #copy over taxids_accids to computer

     #Next, In chatgpt make excel files of the output files: two seperate ones KOs and Clonals for T0 and T96 seperately (combine after):
	
			Ask: can you put these files into an excel sheet where each file is a tab titled the name of the file?
			Attach the 6 files
			check counts after
     ```

2.H) Run Taxon kit to get lineage information    
   * Copy over taxon_link database info from Meso_22 dir   
   * Copy over taxid files into final_taxids dir   

Run taxon kit: 
`sbatch 6.5_run_taxonkit.slurm` this runs --> `./5_run_taxonkit.sh final_taxids`   

example of code:   
`cat final_taxids/FIELD-T0_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> FIELD-T0_linage.txt`   

output will be a dir lineage_files that has all info in it --> copy to computer


2.I) Analyze Taxonomy data in R 
   * First you'll need to process the data in excel - make a total_Lineage files (two separate KO and Clonal) - import each lineage file and adjust data to columns
   	    * Each tab in these file is a location/strain
        * Have chatGPT combine the spreadsheets into KO and Clonal spreadsheets
   * Additionally, make a spreadsheets called (clonal or KO)_genus.csv and copy each location into a column in this sheet - this will be used as the input in R script
   * Make a third spreadsheet of the genetic composition of the viruses present using the VMR spreadsheet using xlookup
        * Run two xlookups - one using the species column, the other using virus name column to conduct your search
        * Then run an if statement to merge the data: `=IF(V2=0,W2,V2)` where V2 = species col and W2 = virus name col
        * remove s__ by find and replace
        * Then do a find and replace of 0 (entire cell) with Unknown (entire workbook)


Now run R scripts to look at taxanomic overlaps and run diversity statistics: Collapsed_genus_work.R


 ## 7) Viral Functional Analysis       


 
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
 
  
## 9) 16s rRNA Analysis (QIIME2)   
	
