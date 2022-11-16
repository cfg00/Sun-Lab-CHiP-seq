#this script just downloads the genome we need, this is more of a test
if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.fna ]
	then
		if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.fna.gz ]
		then
			wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
		fi
		gzip -d GCF_000001405.40_GRCh38.p14_genomic.fna.gz