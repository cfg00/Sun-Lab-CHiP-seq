#!/bin/bash

#important variables

fqDir=""
oDir=""
indexF=0
fqFlag=0
oFlag=0
align=0
readLength=50
quantify=1
rFlag=0
pFlag=0
paired=0
pair1="_1"
pair2="_2"

hFlag = 0

threadFlag=0
threads=0
qReport=0
trim=0

rsd=$(pwd) #directory where RNA seq commands are stored, may not be root of fastq reads

#bash <file.sh> -f <folder where fastq files are> -o <folder where you want output> -t <number of threads>
#-r read length, matters if illumina seq/ alignment,-p paired reads or not (mostly everything is paired),
for var in "$@"
do
  if [[ $fqFlag == 1 ]]
	then
		fqFlag=0
		fqDir="$var"
		#quantify=2
	elif [[ $oFlag == 1 ]]
	then
		oFlag=0
		oDir="$var"/
	elif [[ $rFlag == 1 ]]
	then
		rFlag=0
		readLength=$var
	elif [[ $threadFlag == 1 ]]
	then
		threadFlag=0
		threads=$var
	elif [[ "$var" == "-"* ]]
	then
		if [[ "$var" == *"-a"* ]]
		then	
			align=1 
		fi

		if [[ "$var" == *"-f"* ]]
		then
			fqFlag=1
		fi

		if [[ "$var" == *"-o"* ]]
		then
			oFlag=1
		fi

		if [[ "$var" == *"-r"* ]]
		then
			rFlag=1
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

		if [[ "$var" == "-Q"* ]]
		then
			qReport=1
		fi
		
		if [[ "$var" == "-T"* ]]
		then
			trim=1
		fi

		#adding homer flag
		if [[ "$var" == "-h"* ]]
		then
			hFlag=1
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


if [[ $qReport != 0 ]]
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


		mkdir quality
		chmod -R 0777 quality
		
		for filename in *.fq
		do
			echo "assessing quality of $filename"
			fastqc -o quality $filename
		done
		
		cd $rsd	
fi

if [[ $trim != 0 ]]
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
	
	for filename in *${pair1}.fq
	do
		if [ -f "${filename%${pair1}.fq}${pair2}.fq" ]
		then
			echo "cleaning $filename and ${filename%${pair1}.fq}${pair2}.fq"
			fastp -i "${filename}" -o "${filename%${pair1}*}${pair1}_clean.fq" -I "${filename%${pair1}.fq}${pair2}.fq" -O "${filename%${pair1}*}${pair2}_clean.fq"
		else
			"ERROR: no matching pair for $filename ; will not be cleaned/included"
		fi	
	done
		
	for filename in *_clean.fq
	do	
		echo "removing ${filename%_clean.fq}.fq"
		rm "${filename%_clean.fq}.fq"
		
		echo "assessing quality of cleaned file $filename"
		fastqc -o quality $filename
		
		echo "renaming $filename to ${filename%_clean.fq}.fq"
		mv $filename "${filename%_clean.fq}.fq"
	done
fi

if [[ $align != 0 ]] 
then
#cd "$fqDir"
	#for prefix in Oct4 Klf4 Sox2 H3K27ac H3K4me2 input
	#do
    #bwa mem -t 6 /datasets/cs185-sp22-a00-public/genomes/GRCm38.fa /datasets/cs185-sp22-a00-public/lab4/${prefix}.esc.fastq | samtools view -bS > ${prefix}.bam
	#done


	#my colleague jordan wrote most of this.
	if ! [ -d refGen ] #checks if it a folder named this exits
	then
		mkdir refGen
	fi
	cd refGen
	if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.fna ]
	then
		if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.fna.gz ]
		then
			wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
		fi
		gzip -d GCF_000001405.40_GRCh38.p14_genomic.fna.gz
		#use mv command to change to .fa
	fi
#	if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.gtf ] might need this, check later
#	then
#		if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.gtf.gz ]
#		then
#			wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
#		fi
#		gzip -d GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
#	fi

	cd ../

	cd $fqDir

	mkdir ${oDir%/}
	chmod -R 0777 ${oDir%/}
	mkdir "${oDir%/}/Alignment"
	chmod -R 0777 "${oDir%/}/Alignment" #no crash
	for filename in *${pair1}.fq 
	do

		if [ -f "${filename%${pair1}.fq}${pair2}.fq" ]
		then
		base =  "${filename%${pair1}.fq}"

		read1 = $filename
		read2 = "${filename%${pair1}.fq}${pair2}.fq"

		bwa mem -t $threads ${rsd}/refGen/GCF_000001405.40_GRCh38.p14_genomic.fna $read1 $read2 | samtools view -bS > ${oDir%/}${base}.bam
		

    #bwa mem -t 6 /datasets/cs185-sp22-a00-public/genomes/GRCm38.fa /datasets/cs185-sp22-a00-public/lab4/${prefix}.esc.fastq | samtools view -bS > ${prefix}.bam
	done
	## next step is to add Peak finds, using homer

	


fi 

##how to do it: add command argument i.e -h (run homer) 
#also add bedGraph(done), findPeaks, PeakAnnotations (Q4), find motifs
#merge Peaks

if [[ $hFlag != 0 ]] 
then
	cd $fqDir
	cd $oDir
	#this loop makes the tag directories, using HOMER
	
	mkdir "tagDirectory"
	chmod -R 0777 "tagDirectory" #no crash
	for file in *.bam
	do
		echo "making tag directory for $file" 
		makeTagDirectory "tagDirectory/${file%.bam}" file

		echo "making the bedGraph for $file"
		makeUCSDfile "tagDirectory/${file}" -o auto

		



	done
	#TAG DIRECTORIES FINISHED
	#bedGraph FILES FINISHED :D


	##start bedGraph process:
	#for prefix in Oct4 Klf4 Sox2 H3K27ac H3K4me2
	#do
    #makeUCSCfile ~/lab4-spring22/tagdirs/${prefix} -o auto
	#done

	

fi