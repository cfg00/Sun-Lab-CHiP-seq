#include <iostream>
using namespace std;

// gcc filename.cpp -o test.exe
#include <string>;
#include <cstdlib>;
std::string fromChar(char* s)
{
	int i = 0;
	std::string ret;
	while (s[i] != '\0')
	{
		ret.push_back(s[i]);
	}
	return ret;
}
std::string subString(std::string &original, int start = 0, int end = -1)
{
	//if end is -1, we use the end of the string
	std::string n;
	if (end == -1) end = original.length();
	for (int i = start; i < end; i++) n.push_back(original[i]);
	return n;
}
bool contains(std::string &text, std::string subtext)
{
	for (int i = 0; i < text.length(); i++)
	{
		bool match = true;
		for (int j = 0; j < subtext.length(); j++)
		{
			if (text[i] != subtext[j]) match = false;
			else j = subtext.length();
		}
		if (match) return true;
	}
	return false;
}
int main(int argc, char** argv)
{
	/*
	Command line arguments

	-f fastq directory, all files inside this will be aligned
	
	
	
	
	
	*/
	bool flagF = false; //true when the next command line argument is the fastq directory
	std::string fqDir;
	bool qReport = false;
	bool trim = false;
	bool align = false;
	int readLength = 150;
	bool readLengthFlag = false;
	int pairFlag = 0;
	std::string pair1;
	std::string pair2;
	bool flagO = false; 
	std::string oDir;
	for (int i = 0; i < argc; ++i)
	{
		std::string current = fromChar(argv[i]);
		if (flagF)
		{
			fqDir = current;
			flagF = false;
		}
		else if (readLengthFlag)
		{
			readLength = std::stoi(argv[i]);
			readLengthFlag = false;
		}if (pairFlag)
		{
			pair1 = fromChar(argv[i]);
			if(i!=argc)pair2 = fromChar(argv[++i]);
			pairFlag = 0;
		}if (flagO)
		{
			oDir = fromChar(argv[i]);
			flagO = false;
		}
		else {
			//a flag is not currently active so we check to see if a new flag will be made
			if (contains(current, "-f"))
			{
				flagF = true;
			}
			if (contains(current, "-Q"))
			{
				qReport = true;
			}
			if (contains(current, "-T"))
			{
				trim = true;
			}
			if (contains(current, "-a"))
			{
				align = true;
			}
			if (contains(current, "-r"))
			{
				readLengthFlag = true;
			}
			if (contains(current, "-p"))
			{
				pairFlag = 2;
			}
			if (contains(current, "-o"))
			{
				flagO = true;
			}
			
		}


		
		
	}

	if (qReport)
	{
		/*
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
		
		
		*/
		std::string command;
		command += "cd \"" + fqDir + "\";";
		command += "for filename in* .fastq.gz;";
		command += "do;";
		command += "\tmv $filename ${ filename % .fastq.gz }.fq.gz;";
		command += "done;";

		command += "for filename in *.fq.gz;";
		command += "do;";
		command += "\tif[!- f \"${filename%.*.*}.fq\"];";
		command += "\t\tthen;";
		command += "\t\t\techo \"unzipping $filename\";";
		command += "\t\t\tgzip - dc $filename > \"${filename%.*.*}.fq\";";
		command += "\t\tfi;";
		command += "\tdone;";

		command += "mkdir quality;";
		command += "chmod - R 0777 quality;";

		command += "for filename in* .fq;";
		command += "do;";
		command += "\techo \"assessing quality of $filename\";";
		command += "\tfastqc - o quality $filename;";
		command += "done;";

		std::system((const char*)command.c_str());

	}
	if (trim)
	{
		/*
				if [[ $trim != 0 ]]
		then
			#todo
			#${pair1}
	
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
		*/
		std::string command;
		command += "cd \"" + fqDir + "\";";
		command += "for filename in *.fastq.gz;";
		command += "do;";
		command += "\tmv $filename ${filename % .fastq.gz}.fq.gz; ";
		command += "done;";
		
		command += "for filename in* .fq.gz;";
		command += "do;";
		command += "\tif[!- f \"${filename%.*.*}.fq\"];";
		command += "\tthen;";
		command += "\t\techo \"unzipping $filename\";";
		command += "\t\tgzip - dc $filename > \"${filename%.*.*}.fq\";";
		command += "\tfi;";
		command += "done;";

		command += "for filename in* ${pair1}.fq; ";
		command += "do;";
		command += "\tif[-f \"${filename%${pair1}.fq}${pair2}.fq\"];";
		command += "\tthen;";
		command += "\t\techo \"cleaning $filename and ${filename%${pair1}.fq}${pair2}.fq\";";
		command += "\t\tfastp - i \"${filename}\" - o \"${filename%${pair1}*}${pair1}_clean.fq\" - I \"${filename%${pair1}.fq}${pair2}.fq\" - O \"${filename%${pair1}*}${pair2}_clean.fq\";";
		command += "\telse;";
		command += "\t\t\"ERROR: no matching pair for $filename ; will not be cleaned/included\";";
		command += "\tfi;";
		command += "done;";
		
		command += "for filename in* _clean.fq; ";
		command += "do;";
		command += "\techo \"removing ${filename%_clean.fq}.fq\";";
		command += "\trm \"${filename%_clean.fq}.fq\";";
		
		command += "\techo \"assessing quality of cleaned file $filename\";";
		command += "\tfastqc - o quality $filename;";

		command += "\techo \"renaming $filename to ${filename%_clean.fq}.fq\";";
		command += "\tmv $filename \"${filename%_clean.fq}.fq\";";
		command += "done;";

		std::system((const char*)command.c_str());

	}
	if (align)
	{
		std::string command;
		command += (std::string)"align.sh -r " + std::to_string(readLength) + (std::string)" -p " + pair1 + (std::string)" " + pair2 + " -o " + (std::string)oDir;
		/*
		Todo produce align.sh
		*/ 
		/*done
		
		*/
	}
	//make tag directories for every file presented
	if(true){
		//makeTagDirectory <output directory> <input BAM File>
	}

	//Visualizing the data w IGV:

	if(true){
		//i need arg 1 and arg 2 to be made available, so I can actually run this
	string command3 = "python makeUCSCfile.py " + arg1 + " " + arg2;
	const char * com3 = command3.c_str();
	system(com3);
	}

	//Identifying binding sites:
		//find peaks <tag directory> -i <control tag> - style
		//post2bed.pl peaks.txt > <>.peaks.bed

	//visualizing binding patterns
	if(true){
		string command4 = "findpeaks "+ tagDir+ "-i " + controlTag + "-style " + <h or f> + "-o auto";

		string command5 = "pos2bed.pl" + peaks.txt > oct4.peaks.bed;
	}


	

	

	return 0;
}
