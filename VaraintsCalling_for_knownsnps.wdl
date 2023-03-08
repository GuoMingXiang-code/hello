version 1.0

import "Task.wdl" as All_Task

## variants calling for generating known snps

workflow VaraintsCalling_for_knownsnps {
	meta {
		description: " Preprocessing the paired fq files through trimatic, bwa mem and markduplicate, finally generate knowsites vcf record for bqsr"
	}
	input {
		File fastq1_files
		File fastq2_files
		File sample_file
		File RG_file
		File RefFasta
		File RefIndex
		File RefDict
		File Refamb
		File Refann
		File Refbwt
		File Refpac
		File Refsa
		String species
		
		File scattered_calling_intervals_list
		
		# Runtime parameters
		String docker
		Int NUM_THREAD
		Int MEMORY
		Int DISK
	}
	
	Array[File] Fastq1 = read_lines(fastq1_files)
	Array[File] Fastq2 = read_lines(fastq2_files)
	Array[String] SampleNam = read_lines(sample_file)
	Array[String] RG_array = read_lines(RG_file)
	Array[File] scattered_calling_intervals=read_lines(scattered_calling_intervals_list)
	
	scatter (i in range(length(Fastq1))) {
		call All_Task.preprocessing {
			input :
				docker=docker,
				NUM_THREAD=NUM_THREAD,
				MEMORY=MEMORY,
				DISK=DISK,
				RefFasta=RefFasta,
				RefIndex=RefIndex,
				RefDict=RefDict,
				Refamb=Refamb,
				Refann=Refann,
				Refbwt=Refbwt,
				Refpac=Refpac,
				Refsa=Refsa,

				sample_name=SampleNam[i],
				fastq_1=Fastq1[i],
				fastq_2=Fastq2[i],
				RG=RG_array[i]
		
		}
		
		# bam to vcf per file, split genome into 40 regions
		scatter(interval_file in scattered_calling_intervals) {	
			call All_Task.Haplotypecaller {
				input :
					input_bam=preprocessing.output_bam,
					input_bam_index=preprocessing.output_bam_index,
					interval_list=interval_file,
					ref_dict=RefDict,
					ref_fasta=RefFasta,
					ref_fasta_index=RefIndex,
					
					docker=docker,
					NUM_THREAD=4,
					MEMORY=15,
					DISK=80			

			}

		}
		
		# merge 40 gvcf per sample
		call All_Task.MergeGVCFs {
			input:
				input_vcfs=Haplotypecaller.output_gvcf,
				input_vcfs_indexes=Haplotypecaller.output_gvcf_index,
				output_filename=SampleNam[i],
				docker=docker,
				NUM_THREAD=4,
				MEMORY=100,
				DISK=500
		}

	}

		## 第二部分
		## joint calling based split genome into 40 regions
		scatter (interval_file in scattered_calling_intervals) {
		call All_Task.ImportGVCFs {
		input:
			gvcf_file_list=MergeGVCFs.output_vcf,
			gvcf_file_list_index=MergeGVCFs.output_vcf_index,
			interval=interval_file,
			ref_fasta=RefFasta,
			ref_fasta_index=RefIndex,
			ref_dict=RefDict,
			docker=docker,
			DISK=1000,
			MEMORY=200,
			NUM_THREAD=4
		}
		call All_Task.GenotypeGVCFs {
		input:
			conbineGVCF_file=ImportGVCFs.output_combineGVCF,
			conbineGVCF_file_index=ImportGVCFs.output_combineGVCF_index,
			interval=interval_file,
			ref_fasta=RefFasta,
			ref_fasta_index=RefIndex,
			ref_dict=RefDict,
			docker=docker,
			DISK=500,
			MEMORY=100,
			NUM_THREAD=4
			}
		}

		call All_Task.MergeVCFs {
		input:
			input_vcfs=GenotypeGVCFs.output_vcf,
			input_vcfs_indexes=GenotypeGVCFs.output_vcf_index,
			output_filename=species + ".vcf.gz",
			docker=docker,
			DISK=1000,
			MEMORY=30,
			NUM_THREAD=4
		}


		## 第三部分
		# filter_for_bqsr_knowsites
		call All_Task.filter_for_bqsr_knowsites {
		input:
			docker=docker,
			DISK=1000,
			MEMORY=100,
			NUM_THREAD=6,
			RefFasta=RefFasta,
			RefIndex=RefIndex,
			RefDict=RefDict,
			Refamb=Refamb,
			Refann=Refann,
			Refbwt=Refbwt,
			Refpac=Refpac,
			Refsa=Refsa,
			species=species,
			
			vcf=MergeVCFs.output_vcf,
			vcf_index=MergeVCFs.output_vcf_index
	  
		}	

		## workflow output

		output {
			File for_bqsr_snps=filter_for_bqsr_knowsites.bqsr_snps
			File for_bqsr_snps_index=filter_for_bqsr_knowsites.bqsr_snps_index
			File for_bqsr_indels=filter_for_bqsr_knowsites.bqsr_indels
			File for_bqsr_indels_index=filter_for_bqsr_knowsites.bqsr_indels_index
			
			Array[File] sorted_bam=preprocessing.output_bam
			Array[File] sorted_bam_index=preprocessing.output_bam_index
			Array[File] duplicate_metrics=preprocessing.duplicate_metrics
		}	

}
