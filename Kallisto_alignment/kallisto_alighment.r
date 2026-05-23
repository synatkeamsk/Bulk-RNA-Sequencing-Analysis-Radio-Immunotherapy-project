# Alighment of fastq file to reference genome using Kallisto

#' Building an index
cd kallisto/tests

kallisto index -i transcripts.idx transcripts.fasta.gz

#' create sample list file
sample1
sample2
sample3
sample4
sample5
sample6
sample7
sample8
sample9
sample10
sample11
sample12

#' run kallisto loop (corrected version)
mkdir -p logs

while read -r sample
do
    echo "Processing ${sample}"

    kallisto quant \
        -i transcripts.idx \
        -o kallisto_results/${sample} \
        -b 100 \
        fastq/${sample}_R1.fastq.gz \
        fastq/${sample}_R2.fastq.gz \
        > logs/${sample}.log 2>&1

done < samples.txt

