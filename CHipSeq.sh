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

threadFlag=0
threads=0
qReport=0
trim=0

rsd=$(pwd) #directory where RNA seq commands are stored, may not be root of fastq reads

for var in "$@"
do
  if [[ $fqFlag == 1 ]]
	then
		fqFlag=0
		fqDir="$var"
		quantify=2
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
			if [ ! -f "${filename.fq%.*}" ]
			then
				echo "unzipping $filename"
				gzip -dc $filename
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
		if [ ! -f "${filename.fq%.*}" ]
		then
			echo "unzipping $filename"
			gzip -dc $filename
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
