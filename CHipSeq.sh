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

#homer flag
hFlag = 0

#findPeaks flag
fpFlag = 0

#findPeaks Histone flag
histFlag = 0
histName = ""
#findPeaks TF flag
transFlag = 0

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
	#adding the varname to histName, this is the part that confused me
	elif [[ $histFlag == 1 ]]
	then
		histFlag=0
		histName=$var
		IFS=$'\n' read -d '' -r -a histones_list < $var

	elif [[ $transFlag == 1 ]]
	then
		transFlag=0
		histName=$var
		IFS=$'\n' read -d '' -r -a tf_list < $var
	
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
		if [[ "$var" == "-H"* ]]
		then
			hFlag=1
		fi

		#adding FindPeaks flag

		if [[ "$var" == "-FP"* ]]
		then
			fpFlag = 1
		fi

		#adding Findpeaks sub flags

		#transcription factors
		if [[ "$var" == "-tf"* ]]
		then
			transFlag = 1
		fi

		#histone modifications
		if [[ "$var" == "-hist"* ]]
		then
			histFlag = 1
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

		

		if [[$fpFlag != 0 ]]
		then
			##still confused about this flag loooool
			if [[$histFlag != 0 ]]
			then
			#this is where the findPeaks w/ hist will be run, make a list for histones.
			#IFS=$'\n' read -d '' -r -a histones_list < /etc/passwd
			
			for line in histones_list
			do
				#unsure wether to add i or not...
				#-i ~/lab4-spring22/tagdirs/input
				findPeaks ~/tagDirectory/${line} -i ~/tagDirectory/control -style histone -o auto
				

			done
			

			fi

			if [[$transFlag != 0 ]]
			then
			#this is where the findPeaks w/ TF will be run
			
			for line in tf_list
			do
				#again, not sure wether to add -i or not, ask later
				# 
				   findPeaks ~/tagDirectory/${line} -i ~/tagDirectory/control -style factor -o auto
				   names+=$line

			done


			fi

			if ! [ -d refGen ] #checks if it a folder named this exits, might change the name 
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

			#by this point you have a reference genome called whatever that is ^
			#Make a list with the variable names

			

		#	annotatePeaks.pl tss \
 		#	/refGen/GCF_000001405.40_GRCh38.p14_genomic.fna \
  		#	-size 8000 \
  		#	-hist 10 \
  		#	-d ~/lab4/tagdirs/Oct4 ~/lab4/tagdirs/Sox2 ~/lab4/tagdirs/Klf4 ~/lab4/tagdirs/H3K4me2 ~/lab4/tagdirs/H3K27ac \
  		#	-gtf /datasets/cs185-sp22-a00-public/genomes/GRCm38.75.gtf > ~/lab4/annotations/tss_histogram.txt

		#annotatePeaks.pl 

		cmdstring = ""
		buffer= "~/tagDirectory/"

		for ELEMENT in ${histones_list[@]}; do
  		cmdstring+="${placeholder}"
 	    cmdstring+="${ELEMENT} "
		done

		for line in ${tf_list[@]}; do
   		 cmdstring+="${placeholder}"
    	 buffer+="${line} "
		done
		
		#annotating peaks lol
		annotatePeaks.pl tss\
		/refGen/GCF_000001405.40_GRCh38.p14_genomic.fna \
		-size 8000 \
		-hist 10 \
		-d $cmdstring \
		-gtf /refGen/GCF_000001405.40_GRCh38.p14_genomic.gtf > ~/tagDirectory/annotations/annotated_output.txt

			
		#now we find motifs :D
		#prefix=Oct4
		#findMotifsGenome.pl \
 		#~/lab4/tagdirs/${prefix}/peaks.txt \
 		#/datasets/cs185-sp22-a00-public/genomes/GRCm38.fa \
 		#~/lab4/motifs/${prefix} \
 		#-mask -size 100

		for hist in ${histones_list[@]};
		do
   		findMotifsGenome.pl \
		~/tagDirectory/${hist}/annotated_output.txt \
		/refGen/GCF_000001405.40_GRCh38.p14_genomic.fna \
		/motifs/${hist} \
		-mask -size 100
		done
		
		for trf in ${tf_list[@]};
		do
   		findMotifsGenome.pl \
		~/tagDirectory/${trf}/annotated_output.txt \
		/refGen/GCF_000001405.40_GRCh38.p14_genomic.fna \
		/motifs/${trf} \
		-mask -size 100
		done

		fi



	done

	
	#TAG DIRECTORIES FINISHED
	#bedGraph FILES FINISHED 
	#Making FindPeaks thing.
	#done with FindPeaks 11/16/2022 i think.
	#added find motifs


	#loop through TF file and HM file


	##start bedGraph process:
	#for prefix in Oct4 Klf4 Sox2 H3K27ac H3K4me2
	#do
    #makeUCSCfile ~/lab4-spring22/tagdirs/${prefix} -o auto
	#done
fi
