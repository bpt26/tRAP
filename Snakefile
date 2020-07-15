rule all:
	input:
		"{root_name}tRNAScores.txt"
rule makeHiConfBed:
	input:
		"{root_name}-tRNAs-confidence-set.out"
		"{root_name}-tRNAs.bed"
		"{root_name}.chrom.sizes"
	shell:
	    python makeHiConfBed.py -i input[0] -b input[1] -l input[2]
rule tRNAFasta:
	input:
		"{root_name}.fa" # this is a whole genome fasta file
	shell:
		python tRNAFasta.py -b tRNAHiConf350.bed -g input[0]
rule getMFE:
	input:
		"{root_name}-tRNAs-confidence-set.ss"
	shell:
		python getMFE.py -s input[0] -b tRNAHiConf.bed
rule RNAfold:
	shell:
		RNAfold --noPS -C < tRNAFoldIn.txt > tRNA.mfe
rule classify:
	input:
		"{root_name}-tRNAs-confidence-set.out"
		"{root_name}.chrom.sizes"
		"{root_name}tRNAScores.txt"
	shell:
		python classifytRNAs.py -b tRNAHiConf.bed -e tRNAHiConf350.bed -t input[0] -m tRNA.mfe -f tRNAHiConf350.fa -l input[1] -d humanSimplifiedTrainingData.tsv -o input[2] -x True