# read count in 1Mb bin from the alignment file
samtools view file.bam | cut -f 3,4 | sort -k1,1 -k2,2n | awk '{print $1 "\t" int($2 / 1000000) * 1000000}' | uniq -c | awk 'BEGIN{OFS="\t"}{$4=$2"_"$3}{print}' > 1Mb.file.bam..txt

# total bins in the genome (e.g. Jagger_1MbBin_mod.bed) as:
#chr1A   0       1000000 chr1A_0
#chr1A   1000000 2000000 chr1A_1000000

# any bins without reads should be zero
cat <(cut -f4 Jagger_1MbBin_mod.bed) <(cut -f 4 1Mb.file.bam..txt) | sort | uniq -c | awk '$1==1 {print $NF}' | awk 'BEGIN{FS="_"}{print "0\t"$1"\t"$2"\t"$1"_"$2}' | cat 1Mb.file.bam..txt - | sort -k2,2 -k3,3n > ga.1Mb.file.bam.txt

# coverage based on unique concordant mapped reads: (Jagger assembly size 14552150998)
cu_broken_bp=$(head -4 hisat2_mapping.log | tail -1 | awk '{print $1*302}')
cu_sample_cov=$(echo $cu_broken_bp | awk '{cov=sprintf("%.4f", $1/14552150998)}{print cov}')

# to print
total_pe_reads=$(head -1 hisat2_mapping.log | awk '{print $1}')
cu_reads=$(head -4 ${hisat2_mapping.log | tail -1 | awk '{print $1,$2}')

# normalize
awk -v CU_COV="$cu_sample_cov" -v SAMPLE="$sample_name" 'BEGIN{OFS="\t"}{$5=sprintf("%.4f", ($1*151)/(CU_COV*1000000)); $6=SAMPLE}{print $1,$2,$3,$5,$6}' ga.1Mb.file.bam.txt > final.ga.1Mb.file.bam.txt

echo -e "\n$sample_name:\nUnique concordant coverage: $cu_sample_cov\nTotal PE reads: $total_pe_reads\nConcordant unique PE reads: $cu_reads" >> hisat2_mapping.log

# reads distribution check:
cu_broken=$(head -4 ${hisat2_mapping.log | tail -1 | awk '{print $1*2}')
dis=$(awk '{sum+=$1}END{print sum}' final.ga.1Mb.file.bam.txt)
if [ "$cu_broken" -eq "$dis" ]; then
        echo -e "\nThe total concordant unique (PE*2) reads is EQUAL to reads distributed in the bins." >> hisat2_mapping.log
else
        echo -e "\nThe total concordant unique (PE*2) reads is NOT-EQUAL to reads distributed in the bins." >> hisat2_mapping.log
fi
