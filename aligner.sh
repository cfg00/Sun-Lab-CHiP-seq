#!/bin/bash

#important variables

fqDir=""
oDir=""
indexF=0
fqFlag=0
oFlag=0
align=0
readLength=50

threadFlag=0
threads=0

pFlag=0
paired=0
pair1="_1"
pair2="_2"





rsd=$(pwd) #directory where RNA seq commands are stored, may not be root of fastq reads

#usage of command line arguments
#-h - help //todo
#-f fastqdirectory
#-o output directory

#-r readLength
#-p signals paired ends (should have the same name and end in _1 and _2, otherwise type what files will end in)

#-index (compute index) usage of index with no other arguments will cause only indexing to occur. without this flag indexing will not occur //depricated
#-s slr file, specifications from graphical interface
#input command line arguments
for var in "$@"
do
	if [[ $fqFlag == 1 ]]
	then
		fqFlag=0
		fqDir="$var"
		quantify=2
		#next check for other flags before -
	elif [[ $oFlag == 1 ]]
	then
		oFlag=0
		oDir="$var"/
	
	
	
	
	elif [[ $threadFlag == 1 ]]
	then
		threadFlag=0
		threads=$var
	
	
	elif [[ "$var" == "-"* ]]
	then
		#some flags may have multiple phrases following them so the -option comes first in the else ifs. 
		#unflag these first
		if [[ $pFlag != 0 ]]
		then
			pFlag=0
		fi
		if [[ "$var" == *"-f"* ]]
		then
			#fast q files coming up
			fqFlag=1
		fi
		if [[ "$var" == *"-o"* ]]
		then
			oFlag=1
		fi
		
		if [[ "$var" == *"-t"* ]]
		then
			threadFlag=1
		fi
		
		
		
		if [[ "$var" == *"-p"* ]]
		then
			pFlag=1
			paired=1
		fi
		
		
		if [[ "$var" == "-"*"a"* ]]
		then
			align=2
		fi
		
		
		
		
	elif [[ $pFlag == 1 ]]
	then
		pFlag=2
		pair1="$var"
	elif [[ $pFlag == 2 ]]
	then
		pFlag=0
		pair2="$var"
	else
		#default choice to assume fqdir if fqdir is empty
		if [[ "$fqDir" == "" ]]
		then
			fqDir="$var"
		fi
	fi
done
	
    



    if [[ $paired == 1 ]]
	then
		cd "$fqDir"
		for filename in *.fastq.gz
		do
			mv $filename ${filename%.fastq.gz}.fq.gz
		done

		for filename in *.fq.gz
		do
			if [ ! -f "${filename%.*.*}.fq" ]
			then
				echo "unzipping $filename"
				gzip -dc $filename > "${filename%.*.*}.fq"
			fi
		done
		it=0
		
		for filename in *${pair1}.fq
		do
			if [ -f "${filename%${pair1}.fq}${pair2}.fq" ]
			then
				samples[$it]="${filename%${pair1}.fq}"
				it=$(expr $it + 1)
			else
				"ERROR: no matching pair for $filename ; will not be included"
			fi	
		done
		mkdir ${oDir%/}
		chmod -R 0777 ${oDir%/}
		mkdir "${oDir%/}/Alignment"
		chmod -R 0777 "${oDir%/}/Alignment"

		if [[ $align != 0 ]]
		then
			for i in "${!samples[@]}"
			do 
				base="${samples[$i]}"
				read1=${samples[$i]}${pair1}.fq
				read2=${samples[$i]}${pair2}.fq


				echo "aligning $base"

				STAROPTS="--outSAMattributes NH HI AS NM MD \
					--outFilterType BySJout \
					--outFilterMultimapNmax 20 \
					--outFilterMismatchNmax 999 \
					--outFilterMismatchNoverReadLmax 0.04 \
					--alignIntronMin 20 \
					--alignIntronMax 1000000 \
					--alignMatesGapMax 1000000 \
					--alignSJoverhangMin 8 \
					--alignSJDBoverhangMin 1 \
					--sjdbScore 1 \
					--limitBAMsortRAM 50000000000"
				
				 STAR --genomeDir ${rsd}/refGen/genome$readLength --readFilesIn ${read1} ${read2} --outFileNamePrefix "${oDir}Alignment/${base}" --runThreadN $threads --quantMode TranscriptomeSAM ${STAROPTS}
				 
				 echo "calculating expression of ${base}"
				 rsem-calculate-expression --num-threads $threads --paired-end --alignments "${oDir}Alignment/${base}Aligned.toTranscriptome.out.bam" ${rsd}/refGen/GCF_000001405.40_GRCh38.p14_genomic "${oDir}${base}"
			done
		fi
		#leave the raw read directory, where all alignments and quantifications have now been saved to
		#cd ../
		cd $rsd
	else
		#unpaired pipeline
		if [[ $align != 0 ]]
		then

			for filename in $fqDir/*.fq.gz 
			do
				echo "unzipping $filename"
				gzip -d ${filename}
			done

			for filename in $fqDir/*.fastq.gz 
			do
				echo "unzipping $filename"
				gzip -d ${filename}
			done

			mkdir ${oDir%/}
			chmod -R 0777 ${oDir%/}
			mkdir "${oDir%/}/$fqDir"
			chmod -R 0777 "${oDir%/}/$fqDir"

			for filename in $fqDir/*.fq 
			do
				echo "aligning $filename"
				STAR --genomeDir refGen/genome --readFilesIn ${filename} --outFileNamePrefix "$oDir${filename%.*}" --runThreadN $threads --quantMode TranscriptomeSAM
			done

			for filename in $fqDir/*.fastq 
			do
				echo "aligning $filename"
				STAR --genomeDir refGen/genome --readFilesIn ${filename} --outFileNamePrefix "$oDir${filename%.*}" --runThreadN $threads --quantMode TranscriptomeSAM
			done

			#quantifying gene expression
			#cd ${oDir%/}

			#it is very likely that the following for loop can be removed entirely
			#for filename in $oDir$fqDir/*.out.sam
			#do
			#	echo "converting $filename to "${filename%.*}.bam""
			#	samtools view --threads $threads -S -b $filename > "${filename%.*}.bam"
			#done

			for filename in $oDir$fqDir/*.toTranscriptome.out.bam
			do
				echo "calculating expression of ${filename}"
				sample_fname="${filename##*/}"
				sample_name="${sample_fname%.*}"
				rsem-calculate-expression --num-threads $threads --alignments "${filename%.*}.bam" refGen/GCF_000001405.40_GRCh38.p14_genomic "${sample_name}"
			done

			#cd ../
			#at this point we will have gene.results files?? and they can be used in the R language with Deseq2
		fi
	fi
fi