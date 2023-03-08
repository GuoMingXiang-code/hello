version 1.0

# import task  form Mytasks.wdl
import "Task.wdl" as ALL_Task

# import sub workflow from  VaraintsCalling_for_knownsnps

import  "VaraintsCalling_for_knownsnps.wdl" as sub_VFK


## call varaints based on known snp and then imputation

workflow VaraintsCalling_basedon_knownsnps {
	meta {
		description: " Preprocessing the paired fq files through trimatic, bwa mem and markduplicate, finally generate knowsites vcf record for bqsr"
	}

	
	input {
			File fastq1_files
			File fastq2_files
			File sample_file
			File RefFasta
			File RefIndex
			File RefDict
			File Refamb
			File Refann
			File Refbwt
			File Refpac
			File Refsa
			String species
			File RG_file
			
			File scattered_calling_intervals_list
			
			# Runtime parameters
			String docker
			Int NUM_THREAD
			Int MEMORY
			Int DISK
	}
	
	Array[String] SampleNam = read_lines(sample_file)
	Array[String] RG_array = read_lines(RG_file)
	Array[File] scattered_calling_intervals=read_lines(scattered_calling_intervals_list)

	# bult known snps
	# output:  bam file for per sample; 
	call sub_VFK.VaraintsCalling_for_knownsnps {
		input :
			fastq1_files=fastq1_files,
			fastq2_files=fastq2_files,
			sample_file=sample_file,
			RefFasta=RefFasta,
			RefIndex=RefIndex,
			RefDict=RefDict,
			Refamb=Refamb,
			Refann=Refann,
			Refbwt=Refbwt,
			Refpac=Refpac,
			Refsa=Refsa,
			species=species,

			RG_file=RG_file,
			
			scattered_calling_intervals_list=scattered_calling_intervals_list,
			
			# Runtime parameters
			docker=docker,
			NUM_THREAD=NUM_THREAD,
			MEMORY=MEMORY,
			DISK=DISK
	}
	# BQSR; input,sorted markDuplicated bam perl sample; based on knowSites
	# bam Array file from subworkflow
	scatter (j in range(length(VaraintsCalling_for_knownsnps.sorted_bam))) {
		call ALL_Task.run_BQSR {
		input:
		docker = docker,
		NUM_THREAD = 5,  
		MEMORY = 50,
		DISK = 500,
		
		RefFasta = RefFasta,
		RefIndex = RefIndex,
		RefDict = RefDict,
		Refamb = Refamb,
		Refann = Refann,
		Refbwt = Refbwt,
		Refpac = Refpac,
		Refsa = Refsa,
		
		sample_name = SampleNam[j],
		input_bam=VaraintsCalling_for_knownsnps.sorted_bam[j],
		input_bam_index=VaraintsCalling_for_knownsnps.sorted_bam_index[j],
		
		knownSites_snps=VaraintsCalling_for_knownsnps.for_bqsr_snps,
		knownSites_snps_index=VaraintsCalling_for_knownsnps.for_bqsr_snps_index,
		knownSites_indels=VaraintsCalling_for_knownsnps.for_bqsr_indels,
		knownSites_indels_index=VaraintsCalling_for_knownsnps.for_bqsr_indels_index
		}
		
		
		# HaplotypeCaller, base on 40 genome region; 
		# bam to vcf per file, split genome into 40 regions
		
		scatter(interval_file in scattered_calling_intervals) {
			call ALL_Task.Haplotypecaller{
			input :
			input_bam=run_BQSR.basrBAM,
			input_bam_index=run_BQSR.basrBAI,
			interval_list=interval_file,
			ref_dict=RefDict,
			ref_fasta=RefFasta,
			ref_fasta_index=RefIndex,
			
			docker=docker,
			NUM_THREAD=4,
			MEMORY=50,
			DISK=100
			}
		
		}

		# merge 40 gvcf per sample to one vcf
		call ALL_Task.MergeGVCFs {
		input:
		input_vcfs=Haplotypecaller.output_gvcf,
		input_vcfs_indexes=Haplotypecaller.output_gvcf_index,
		output_filename=SampleNam[j],
		docker=docker,
		NUM_THREAD=4,
		MEMORY=15,
		DISK=30
		}
	}	
	## 第二部分
	## joint calling based split genome into 40 regions
	scatter (interval_file in scattered_calling_intervals) {
		call ALL_Task.ImportGVCFs { 
		input:
		gvcf_file_list=MergeGVCFs.output_vcf,
		gvcf_file_list_index=MergeGVCFs.output_vcf_index,
		interval=interval_file,
		ref_fasta=RefFasta,
		ref_fasta_index=RefIndex,
		ref_dict=RefDict,
		docker=docker,
		DISK=DISK,
		NUM_THREAD=4,
		MEMORY=80
		}                   
		call ALL_Task.GenotypeGVCFs {
		input:
		conbineGVCF_file=ImportGVCFs.output_combineGVCF,
		conbineGVCF_file_index=ImportGVCFs.output_combineGVCF_index,
		interval=interval_file,
		ref_fasta=RefFasta,
		ref_fasta_index=RefIndex,
		ref_dict=RefDict,
		docker=docker
		}
	}
	
	call ALL_Task.MergeVCFs {
	input:
	input_vcfs=GenotypeGVCFs.output_vcf,
	input_vcfs_indexes=GenotypeGVCFs.output_vcf_index,
	output_filename=species + ".vcf.gz",
	docker=docker
	}

	# varaints filter and VQSR for per sample vcf
	# "run VQSR for vcf files, include VariantsFilteration & SelectVariants, XY chromosome repeatRegion"
	call ALL_Task.run_VQSR {
	input:
		docker =docker,
		NUM_THREAD = 4,
		MEMORY = 30,
		DISK = 500,
		
		RefFasta=RefFasta,
		RefIndex=RefIndex,
		RefDict=RefDict,
		Refamb=Refamb,
		Refann=Refann,
		Refbwt=Refbwt,
		Refpac=Refpac,
		Refsa=Refsa,
		
		species = species,
		merged_vcf = MergeVCFs.output_vcf,
		merged_vcf_index = MergeVCFs.output_vcf_index
	}                   
	
	
	# annotation and imputation
	call ALL_Task.annotation_and_imputation {
	input:
		input_vcf = run_VQSR.final_vcf,
		species = species,
		docker = docker,
		
		NUM_THREAD=NUM_THREAD,
		MEMORY=MEMORY,
		DISK=DISK
	}           
	
	# workflow output	
	output {
		File result_vcf = run_VQSR.final_vcf
		File imputation_vcf = annotation_and_imputation.imputation_vcf
		#File annotation_gzvcf = annotation_and_imputation.annotation_gzvcf
	}           
}
