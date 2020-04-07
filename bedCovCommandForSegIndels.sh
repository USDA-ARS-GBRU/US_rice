module load bedtools/2.25.0
module load samtools/1.9

date
time samtools view rice1.sorted.bam -f 0x1 -F 0x800 -L ../long.bed -h |awk '$7 !~ /Ch.*/' | awk '$6 !~ /.*S/' |  samtools view -Sb - |bedtools coverage -a ../long.bed -b stdin > BedCovOut_long/rice1.out
time samtools view rice1.sorted.bam -f 0x1 -F 0x800 -L ../short.bed -h |awk '$7 !~ /Ch.*/' | awk '$6 !~ /.*S/' |  samtools view -Sb - |bedtools coverage -f 1 -a ../short.bed -b stdin > BedCovOut_short/rice1.out
cat  BedCovOut_long/301018_Blue.out BedCovOut_short/rice1.out |sed 's:^Chr::g'|sort -n -k1,2|sed 's:^:Chr:g'  > BedCovOut_merged/rice1.out
date
