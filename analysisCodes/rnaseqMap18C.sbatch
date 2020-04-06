#! /bin/bash
# run in hyperion:/home/whuang9/geneExpTemp/
# ============================================================

# map RNA-seq reads for 18 and 25 C females and males
# reads are located in 
# sequence/NCSU-GSL-0008 01 = 25 C female, 02 = male
# sequence/NCSU-GSL-0102 01 = 18 C female, 02 = male
# ============================================================

# request resource
# ============================================================

#SBATCH --cpus-per-task=16
#SBATCH --output=/home/whuang9/geneExpTemp/log/rnaseqMap18C.log

# add bowtie2 to path, necessary for tophat to function
# there is also a copy in /usr/local/bin
# run a print version before the command to make sure
# it's using the correct copy
# ============================================================

export PATH=/home/whuang9/software/bowtie2-2.2.3:$PATH
bowtie2 --version

# align sequences for 18 C samples
# ============================================================

/home/whuang9/software/tophat-2.0.13.Linux_x86_64/tophat --read-mismatches 6 --read-gap-length 6 --read-edit-dist 6 --read-realign-edit-dist 0 -o map/female18c --mate-inner-dist -20 --mate-std-dev 80 -i 20 -I 25000 --solexa1.3-quals -p 16 -g 20 --library-type fr-unstranded --segment-mismatches 2 --segment-length 25 --min-segment-intron 20 --max-segment-intron 25000 --no-coverage-search --transcriptome-index /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean/fb-r5.57.clean --rg-id female18c --rg-sample female18c --rg-library female18c --rg-platform illumina /home/whuang9/flybase/fb-r5.57/dmel.rmsk /home/whuang9/sequence/NCSU-GSL-0102/NCSU-GSL-0102_L1/130508_I1089_FCC1VKUACXX_L3_NCSU-GSL-0102-1/130508_I1089_FCC1VKUACXX_L3_NCSU-GSL-0102-1_1.fq.gz /home/whuang9/sequence/NCSU-GSL-0102/NCSU-GSL-0102_L1/130508_I1089_FCC1VKUACXX_L3_NCSU-GSL-0102-1/130508_I1089_FCC1VKUACXX_L3_NCSU-GSL-0102-1_2.fq.gz > /home/whuang9/geneExpTemp/log/female.18C.map.log 2>&1

/home/whuang9/software/tophat-2.0.13.Linux_x86_64/tophat --read-mismatches 6 --read-gap-length 6 --read-edit-dist 6 --read-realign-edit-dist 0 -o map/male18c --mate-inner-dist -20 --mate-std-dev 80 -i 20 -I 25000 --solexa1.3-quals -p 16 -g 20 --library-type fr-unstranded --segment-mismatches 2 --segment-length 25 --min-segment-intron 20 --max-segment-intron 25000 --no-coverage-search --transcriptome-index /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean/fb-r5.57.clean --rg-id male18c --rg-sample male18c --rg-library male18c --rg-platform illumina /home/whuang9/flybase/fb-r5.57/dmel.rmsk /home/whuang9/sequence/NCSU-GSL-0102/NCSU-GSL-0102_L1/130508_I1089_FCC1VKUACXX_L3_NCSU-GSL-0102-2/130508_I1089_FCC1VKUACXX_L3_NCSU-GSL-0102-2_1.fq.gz /home/whuang9/sequence/NCSU-GSL-0102/NCSU-GSL-0102_L1/130508_I1089_FCC1VKUACXX_L3_NCSU-GSL-0102-2/130508_I1089_FCC1VKUACXX_L3_NCSU-GSL-0102-2_2.fq.gz > /home/whuang9/geneExpTemp/log/male.18C.map.log 2>&1
