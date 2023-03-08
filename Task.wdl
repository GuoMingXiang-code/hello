version 1.0

task preprocessing {
	meta {
	description: " Preprocessing the paired fq files through trimatic, bwa mem and markduplicate"
	}
	input{
	String docker
	Int NUM_THREAD 
	Int MEMORY
	Int DISK

	File RefFasta
	File RefIndex
	File RefDict
	File Refamb
	File Refann
	File Refbwt
	File Refpac
	File Refsa

	String sample_name
	File fastq_1
	File fastq_2
	String RG
	}

	command <<<
	set -e

	sample=~{sample_name}
	RG=~{RG}
	temp="split"
	mkdir -p $temp
	ref=~{RefFasta}
	time seqkit split2 -1 ~{fastq_1} -2 ~{fastq_2} -s 1000000 -O $temp -f -j ~{NUM_THREAD}  && echo "**~{sample_name} split done**"
	counts=`ls split|cut -f3 -d '_' |cut -f1 -d '.' |sort -n|tail -n1|awk '{print int($0)}'`
	echo $counts;

	for((k=1;k<=${counts};k++));
	do
	var=$(printf "%03d" "$k")
	read1_split=$temp"/"$sample"_1.part_"$var".fastq.gz"
	read2_split=$temp"/"$sample"_2.part_"$var".fastq.gz"
	trim1=$temp"/"$sample"_1_part_"$var".trim_paried.fastq.gz"
	trim2=$temp"/"$sample"_2_part_"$var".trim_paried.fastq.gz"
	trim1_unparied=$temp"/"$sample"_1_part_"$var".trim_unparied.fastq.gz"
	trim2_unparied=$temp"/"$sample"_2_part_"$var".trim_unparied.fastq.gz"
	untreated=$temp"/untreated_part_"$var".logfile"

	bamFilePart=$temp"/"$sample"_"$var"_part.bam" 
	
	# trimmomatic, 只对paired结果进行比对，输出成bam: trimmomatic instand of trimmomatic`
	time java -jar /software/bin/trimmomatic-0.39.jar \
	PE -phred33 -trimlog $untreated \
	$read1_split $read2_split \
	$trim1 $trim1_unparied $trim2 $trim2_unparied \
	SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50 -threads ~{NUM_THREAD} \
	&& echo "**$sample_$var trimmomatic done**" \
	&& rm $read1_split $read2_split \
	&& time bwa mem -t ~{NUM_THREAD} -M -Y -R $RG $ref $trim1 $trim2 | samtools view -b -S -@ ~{NUM_THREAD} ->$bamFilePart \
	&& echo "**$sample_$var bwa mapping done **"
	done

	# path declaration
	bamFile=$temp"/"$sample".bam"
	bamSortFile=$temp"/"$sample".sorted.bam"
	bamDir=".";
	bamSortMarkupFile=$bamDir"/"$sample".sorted.markup.bam"
	metricFile=$bamDir"/"$sample".sorted.markdup_metrics.txt"
	#pushd $temp
	echo "test $$$$$$$$"
	bamlist=`ls $temp"/"*part.bam | tr '\n' ' '`
	echo $bamlist
	echo "test end $$$$$$$$"
	JavaMem=`echo ~{MEMORY-4}|bc`
	echo "****************JavaMem   $JavaMem *****************"

	# merge，sort and markDuplicates
	time samtools merge $bamFile $bamlist -@ ~{NUM_THREAD} && rm $bamlist  && echo "**$sample merge done **" 
	time samtools sort -@ ~{NUM_THREAD} -O bam -o $bamSortFile $bamFile && rm $bamFile && echo "**$sample sort is done**" 
	ulimit -n 100000
	time gatk  --java-options "-Xmx${JavaMem}G -XX:+UseParallelGC -XX:ParallelGCThreads=~{NUM_THREAD}" MarkDuplicates -I $bamSortFile -O $bamSortMarkupFile -M $metricFile
	echo "**$sample mark duplicates is done**" 
	time samtools index $bamSortMarkupFile -@ ~{NUM_THREAD}

	>>>
	
	runtime {
	docker: docker 
	cpu: "${NUM_THREAD}" 
	memory: "${MEMORY} GB"
	disk: "${DISK} GB"
	}

	output {
	File output_bam="${sample_name}.sorted.markup.bam"
	File output_bam_index="${sample_name}.sorted.markup.bam.bai"
	File duplicate_metrics="${sample_name}.sorted.markdup_metrics.txt"
	}

}


task run_BQSR {
	meta {
	description: "BQSR for the input bam"
	}
	input{
	String docker
	Int NUM_THREAD
	Int MEMORY
	Int DISK

	File RefFasta
	File RefIndex
	File RefDict
	File Refamb
	File Refann
	File Refbwt
	File Refpac
	File Refsa

	String sample_name
	File input_bam
	File input_bam_index
	################# changed  ###########
	File knownSites_snps
	File knownSites_snps_index
	File knownSites_indels
	File knownSites_indels_index	
	}

	command {
	set -e

	time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=~{NUM_THREAD}" BaseRecalibrator \
	-R ${RefFasta} \
	-I ${input_bam} \
	-O ${sample_name}.table \
	--known-sites ${knownSites_snps} \
	--known-sites ${knownSites_indels} \
	--bqsr-baq-gap-open-penalty 30 && echo "** ${sample_name}.recal_data.table done **" \
	&& time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=~{NUM_THREAD}" ApplyBQSR \
	-R ${RefFasta} \
	-I ${input_bam} \
	-bqsr ${sample_name}.table \
	--static-quantized-quals 10 --static-quantized-quals 20 \
	--static-quantized-quals 30 -O ${sample_name}.bqsr.bam &&echo "**${sample_name} ApplyBQSR done **"
	
	samtools index -@ ${NUM_THREAD} ${sample_name}.bqsr.bam	
	}
	
	runtime {
	docker: docker 
	cpu: "${NUM_THREAD}" 
	memory: "${MEMORY} GB"
	disk: "${DISK} GB"
	}

	output {
	File basrBAM="${sample_name}.bqsr.bam"
	File basrBAI="${sample_name}.bqsr.bam.bai"
	}

}

# 以合并后的vcf做为起始
task filter_for_bqsr_knowsites {
	meta {
	description: "generate knowsites vcf record for bqsr"
	}
	input{
	String docker
	Int NUM_THREAD 
	Int MEMORY
	Int DISK

	File RefFasta
	File RefIndex
	File RefDict
	File Refamb
	File Refann
	File Refbwt
	File Refpac
	File Refsa

	String species
	File vcf 
	File vcf_index 

	}

	command {

	organism=~{species}
	ref=~{RefFasta}

	filterDir="filterDir"
	mergeVcf=~{vcf}
	mkdir -p $filterDir	
	
	########## changed #####################
	raw_snps="raw_snps.vcf.gz"
	raw_indels="raw_indels.vcf.gz"
	filtered_snps="filtered_snps.vcf.gz"
	filtered_indels="filtered_indels.vcf.gz"
	bqsr_snps="bqsr_snps.vcf.gz"
	bqsr_indels="bqsr_indels.vcf.gz"
	
	time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" SelectVariants \
		-R $ref \
		-V $mergeVcf \
		--select-type-to-include SNP \
		-O $raw_snps && echo "**SelectVariants  raw snps done **" \
	&&time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" SelectVariants \
			-R $ref \
			-V $mergeVcf \
			--select-type-to-include INDEL \
			-O $raw_indels && echo "**SelectVariants  raw indels done **" \
	&&time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" VariantFiltration \
			-R $ref \
			-V $raw_snps\
			-O $filtered_snps \
			-filter-name "QD_filter" -filter "QD < 2.0" \
			-filter-name "FS_filter" -filter "FS > 200.0" \
			-filter-name "MQ_filter" -filter "MQ < 40.0" \
			-filter-name "SOR_filter" -filter "SOR > 4.0" \
			--filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
			--filter-name "MQRSlt12" \
			--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
			--filter-name "RPRSlt8" \
			--cluster-size 3 \
			--cluster-window-size 20 &&echo "** varialtFiltration SNPs done **" \
	&& time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" VariantFiltration \
			-R $ref \
			-V $raw_indels \
			-O $filtered_indels \
			-filter-name "QD_filter" -filter "QD < 2.0" \
			-filter-name "FS_filter" -filter "FS > 200.0" \
			-filter-name "SOR_filter" -filter "SOR > 10.0" &&echo "** varialtFiltration indels done **" \
	&& time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" SelectVariants \
			--exclude-filtered \
			-V $filtered_snps \
			-O $bqsr_snps && echo "**  Exclude Filtered Variants snps done **" \
	&& time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" SelectVariants \
			--exclude-filtered \
			-V $filtered_indels \
			-O $bqsr_indels && echo "**  Exclude Filtered Variants indels done **" \
	&& echo "** filter done **"

	}
	
	runtime {
		docker: docker 
		cpu: "${NUM_THREAD}" 
		memory: "${MEMORY} GB"
		disk: "${DISK} GB"
	}

	output {
		File bqsr_snps="bqsr_snps.vcf.gz"
		File bqsr_snps_index="bqsr_snps.vcf.gz.tbi"
		File bqsr_indels="bqsr_indels.vcf.gz"
		File bqsr_indels_index="bqsr_indels.vcf.gz.tbi"
	}

}



# 以合并后的vcf做为起始
task run_VQSR {
	meta {
	description: "run VQSR for vcf files"
	}
	input{
	String docker 
	Int NUM_THREAD
	Int MEMORY
	Int DISK

	File RefFasta
	File RefIndex
	File RefDict
	File Refamb
	File Refann
	File Refbwt
	File Refpac
	File Refsa

	String species
	File merged_vcf
	File merged_vcf_index
	}

	command {
	organism=~{species};
	ref=~{RefFasta}
	filterDir="filterDir"
	if [ ! -d $filterDir ];then
	mkdir -p $filterDir;
	fi
	mergeVcf=~{merged_vcf}
	
	##  changed 
	rawSNP="raw_snps_recal.vcf.gz"
	rawINDEL="raw_indel_recal.vcf.gz"
	rawMAXED="All.raw.MAXED.vcf.gz"
	filteredSNPs="filtered_snps_final.vcf.gz"
	filteredINDELs="filtered_indels_final.vcf.gz"
	finalSNP="final_snps.vcf.gz"
	finalINDELs="final_indels.vcf.gz"
	finalVcf="final_both.vcf.gz"
	
	time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" SelectVariants \
			-R $ref \
			-V $mergeVcf\
			--select-type-to-include SNP \
			--output $rawSNP && echo "** Raw SNP  done **" \
	&& time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" SelectVariants \
			-R $ref \
			-V $mergeVcf \
			--select-type-to-include INDEL \
			-O $rawINDEL && echo "** Raw INDEL  done **" \
	&& time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" VariantFiltration \
			-R $ref \
			-O $filteredSNPs \
			-V $rawSNP \
			--filter-expression "FS > 60.0" \
			--filter-name "FS" \
			--filter-expression "MQ < 40.0" \
			--filter-name "MQ" \
			--filter-expression "QD < 2.0" \
			--filter-name "QD" \
			--filter-expression "SOR > 4.0" \
			--filter-name "SOR" \
			--filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
			--filter-name "MQRSlt12" \
			--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
			--filter-name "RPRSlt8" \
			--genotype-filter-expression "GQ < 13.0" \
			--genotype-filter-name "GQ13" \
			--cluster-size 3 --cluster-window-size 20 &&echo "** varialtFiltration SNP done **" \
	&&time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" VariantFiltration \
			-R $ref \
			-O $filteredINDELs \
			-V $rawINDEL \
			--filter-expression "FS > 200.0" \
			--filter-name "FS" \
			--filter-expression "QD < 2.0" \
			--filter-name "QD" \
			--filter-expression "SOR > 10.0" \
			--filter-name "SOR" \
			--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
			--filter-name "RPRSlt8" && echo "** varialtFiltration INDEL done **"  \
	&& time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" SelectVariants \
			--exclude-filtered \
			-V $filteredSNPs\
			--select-type-to-include SNP \
			--output $finalSNP && echo "** Raw SNP  done **" \
	&& time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" SelectVariants \
			--exclude-filtered \
			-V $filteredINDELs\
			--select-type-to-include SNP \
			--output $finalINDELs && echo "** Raw SNP  done **" \
	&& time gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" MergeVcfs \
			-I $finalSNP \
			-I $finalINDELs \
			-O $finalVcf && echo "** MergeVcfs done **" && echo "* VQSR done*"


	}
	
	runtime {
		docker: docker 
		cpu: "${NUM_THREAD}" 
		memory: "${MEMORY} GB"
		disk: "${DISK} GB"
	}

	output {
		File final_vcf="final_both.vcf.gz"
		File final_vcf_index="final_both.vcf.gz.tbi"
	}
}


## A1
task MergeGVCFs {
	input {
		Array[File] input_vcfs
		Array[File] input_vcfs_indexes
		String output_filename
		String docker
		Int NUM_THREAD
		Int MEMORY
		Int DISK
	
	}
	command {
	set -e
	gatk --java-options "-Xmx${MEMORY}G -XX:ParallelGCThreads=${NUM_THREAD}" MergeVcfs \
	--INPUT ${sep=' --INPUT ' input_vcfs} \
	--OUTPUT ${output_filename}.g.vcf.gz
	}
	runtime {
	docker: docker
	cpu: "${NUM_THREAD}"
	memory: "${MEMORY} GB"
	disk: "${DISK} GB"
	}
	output {
	File output_vcf="${output_filename}.g.vcf.gz"
	File output_vcf_index="${output_filename}.g.vcf.gz.tbi"
	}
}


# TASK DEFINITIONS
task Haplotypecaller {
	input {
		File input_bam
		File input_bam_index
		File interval_list
		String output_filename=basename(input_bam,".sorted.markup.bam")
		File ref_dict
		File ref_fasta
		File ref_fasta_index
		Float contamination=0.0
		String docker
		Int NUM_THREAD
		Int MEMORY
		Int DISK
	}
	
	command {
	set -e
	gatk --java-options "-Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" HaplotypeCaller \
	-R ${ref_fasta} \
	-I ${input_bam} \
	-L ${interval_list} \
	-O ${output_filename}.g.vcf.gz \
	-ERC GVCF \
	-stand-call-conf 20 \
	-contamination ${contamination}
	}
	runtime {
		docker: docker
		cpu: "${NUM_THREAD}"
		memory: "${MEMORY} GB"
		disk: "${DISK} GB"	
	}
	output {
		File output_gvcf="${output_filename}.g.vcf.gz"
		File output_gvcf_index="${output_filename}.g.vcf.gz.tbi"
	}
}


task ImportGVCFs {
	input {
		Array[File] gvcf_file_list
		Array[File] gvcf_file_list_index
		File interval
		File ref_fasta
		File ref_fasta_index
		File ref_dict
		String docker
		String output_combineGVCF_name="combineGVCF.g.vcf.gz"

		Int NUM_THREAD
		Int MEMORY
		Int DISK
	}
	command {
	set -euo pipefail
	gatk --java-options "-Xms${MEMORY}G -Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" CombineGVCFs \
	-R ${ref_fasta} \
	-L ${interval} \
	-O ${output_combineGVCF_name} \
	-V ${sep=' -V ' gvcf_file_list}
	}
	
	runtime {
	docker: docker
	cpu: "${NUM_THREAD}"
	memory: "${MEMORY} GB"
	disk: "${DISK} GB"
	}
	output {
	File output_combineGVCF="${output_combineGVCF_name}"
	File output_combineGVCF_index="${output_combineGVCF_name}.tbi"
	}
}

task GenotypeGVCFs {
	input {
		File conbineGVCF_file
		File conbineGVCF_file_index
		File interval
		 
		File ref_fasta
		File ref_fasta_index
		File ref_dict
		 
		 
		String docker
		Int NUM_THREAD
		Int MEMORY
		Int DISK
		String output_combineVCF_name="combineVCF.vcf.gz"

	}
	command {
	set -euo pipefail
	gatk --java-options "-Xms${MEMORY}G -Xmx${MEMORY}G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREAD}" GenotypeGVCFs \
	-R ${ref_fasta} \
	-V ${conbineGVCF_file} \
	-L ${interval} \
	-O ${output_combineVCF_name}
	}
	
	runtime {
	docker: docker
	cpu: "${NUM_THREAD}"
	memory: "${MEMORY} GB"
	disk: "${DISK} GB"
	}
	output {
	File output_vcf="${output_combineVCF_name}"
	File output_vcf_index="${output_combineVCF_name}.tbi"
	}
}

task MergeVCFs {
	input {
		Array[File] input_vcfs
		Array[File] input_vcfs_indexes
		String output_filename
		 
		String docker
		Int NUM_THREAD
		Int MEMORY
		Int DISK
		
	}
	command {
	set -e
	gatk --java-options "-Xmx${MEMORY}G -XX:ParallelGCThreads=${NUM_THREAD}" MergeVcfs \
	--INPUT ${sep=' --INPUT ' input_vcfs} \
	--OUTPUT ${output_filename}
	}
	runtime {
	docker: docker
	cpu: "${NUM_THREAD}"
	memory: "${MEMORY} GB"
	disk: "${DISK} GB"
	}
	output {
	File output_vcf="${output_filename}"
	File output_vcf_index="${output_filename}.tbi"
	}
}



## vcf annotation and imputation for per species
task annotation_and_imputation {
	input {
	File input_vcf
	String species
	String docker
	Int NUM_THREAD
	Int MEMORY
	Int DISK
	}
	command {
	#JavaMem=MEMORY*0.7
	java -Xmx30G -jar /software/bin/beagle.jar gt=${input_vcf} out=${species}_imputation
	#java -Xmx${MEMORY}G -jar /software/bin/snpEff/snpEff.jar -v ${species} ${species}_imputation.vcf.gz > ${species}.imputation.anno.vcf.gz
	} 
	runtime {
		docker: docker
		cpu: "${NUM_THREAD}"
		memory: "${MEMORY} GB"
		disk: "${DISK} GB"
	}
	output {
	File imputation_vcf = "${species}_imputation.vcf.gz"
	#File annotation_gzvcf = "${species}.imputation.anno.vcf.gz"
	}
}
