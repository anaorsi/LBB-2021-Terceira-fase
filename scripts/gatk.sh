#!/usr/bin/bash

cd tools/gatk-workflows

# Indexação e confecção do dicionário
java -jar ./picard/build/libs/picard.jar CreateSequenceDictionary R=./inputs/ref_coconut.fasta O=./inputs/ref_coconut.dict
/usr/local/bin/samtools faidx ./inputs/ref_coconut.fasta

# Genome Generate STAR
mkdir ./inputs/STAR_ref
~/tools/STAR-2.7.9a/source/STAR --runMode genomeGenerate --genomeDir ./inputs/STAR_ref --genomeFastaFiles ./inputs/ref_coconut.fasta --sjdbGTFfile ./inputs/coconut.gtf --limitGenomeGenerateRAM 90607874144 --sjdbOverhang 74 --runThreadN 6


while IFS="" read -r p || [ -n "$p" ]
do

	# Obtenção dos dados
	NAME="$p"
	echo $NAME
	mkdir "./inputs/"$NAME
	cd "./inputs/"$NAME
        ~/tools/sra-tools/sratoolkit.2.9.6-ubuntu64/bin/prefetch $NAME
        ~/tools/sra-tools/sratoolkit.2.9.6-ubuntu64/bin/fasterq-dump $NAME
        cd ../..
	
	# Conversão fastq -> ubam
	java -Xmx8G -jar ./picard/build/libs/picard.jar FastqToSam -FASTQ "./inputs/"$NAME"/"$NAME"_1.fastq" -FASTQ2 "./inputs/"$NAME"/"$NAME"_2.fastq" -OUTPUT "./inputs/"$NAME"/"$NAME"_fastqtosam.bam" -READ_GROUP_NAME Test -SAMPLE_NAME $NAME -LIBRARY_NAME Novaseq -PLATFORM illumina

	# Alinhamento inicial
	mkdir "./inputs/"$NAME"/STAR_ini/"
	~/tools/STAR-2.7.9a/source/STAR --genomeDir ./inputs/STAR_ref --runThreadN 6 --readFilesIn "./inputs/"$NAME"/"$NAME"_1.fastq" "./inputs/"$NAME"/"$NAME"_2.fastq"  --sjdbOverhang 74 --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix "./inputs/"$NAME"/STAR_ini/"
	
	# Chamada de variantes simples, que servirá como input para o GATK
	/usr/local/bin/samtools index "./inputs/"$NAME"/STAR_ini/Aligned.sortedByCoord.out.bam" "./inputs/"$NAME"/STAR_ini/Aligned.sortedByCoord.out.bam.bai"

	echo "Mpileup"
	bcftools mpileup --max-depth 10000 --ff UNMAP,SECONDARY -I -Ou -f "./inputs/ref_coconut.fasta" "./inputs/"$NAME"/STAR_ini/Aligned.sortedByCoord.out.bam" | bcftools call -mv -Ob -V indels -m --ploidy X -o "./inputs/"$NAME"/"$NAME".bcf"
	echo "BCF->VCF"
	bcftools view -i '%QUAL>=50' -o "./inputs/"$NAME"/"$NAME".vcf" "./inputs/"$NAME"/"$NAME".bcf"
	bgzip -c "./inputs/"$NAME"/"$NAME".vcf" > "./inputs/"$NAME"/"$NAME".vcf.gz"
	tabix -p vcf "./inputs/"$NAME"/"$NAME".vcf.gz"
	
	# Limpeza de arquivos desnecessários (liberar espaço no servidor)
	rm -r "./inputs/"$NAME"/STAR_ini/"

	# Edição do arquivo de inputs
	awk '{sub("NAME", "'$NAME'/'$NAME'")}1' ./gatk4-rnaseq-germline-snps-indels/to_replace.json > ./gatk4-rnaseq-germline-snps-indels/new.inputs.json

	# Execução do workflow
	java -jar ./cromwell-33.1.jar run ./gatk4-rnaseq-germline-snps-indels/gatk4-rna-best-practices.wdl --inputs ./gatk4-rnaseq-germline-snps-indels/new.inputs.json
	
done < ./inputs/accession.txt
