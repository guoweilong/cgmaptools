#!/usr/bin/env sh


chmod -w  RR.ATCGmap.gz  RR.bam WG.bam  genome.fa anno.refFlat
samtools index WG.bam
samtools index RR.bam
samtools faidx genome.fa

# ====================

cgmaptools convert bam2cgmap -b WG.bam -g genome.fa --rmOverlap -o WG
cgmaptools convert bam2cgmap -b RR.bam -g genome.fa --rmOverlap -o RR

cgmaptools convert atcgmap2atcgbz -c WG.ATCGmap.gz -b WG.ATCGbz
cgmaptools convert atcgbz2atcgmap -c WG2.ATCGmap.gz -b WG.ATCGbz

cgmaptools convert cgmap2cgbz -c RR.CGmap.gz -b RR.CGbz
cgmaptools convert cgbz2cgmap -c RR2.CGmap.gz -b RR.CGbz

cgmaptools fetch atcgbz -b WG.ATCGbz -C chr2 -L 90 -R 100
cgmaptools fetch cgbz -b RR.CGbz -C chr3 -L 2200 -R 2400

zcat RR2.CGmap.gz | gawk -F"\t" -vOFS="\t" '{$4="-"; $5="-"; print;}' | cgmaptools refill -g genome.fa -o RR3.CGmap.gz

cgmaptools intersect -1 WG.CGmap.gz -2 RR.CGmap.gz -C CG -o intersect_CG.gz

cgmaptools merge2 atcgmap -1 WG.ATCGmap.gz -2 RR.ATCGmap.gz | gzip > merge.ATCGmap.gz
cgmaptools merge2 cgmap -1 WG.CGmap.gz -2 RR.CGmap.gz | gzip > merge.CGmap.gz

zcat RR*.CGmap.gz WG.CGmap.gz | gawk '$8>=5' | cut -f1,3 | sort -u | cgmaptools sort -c 1 -p 2 > index
cgmaptools mergelist tomatrix -i index -f RR.CGmap.gz,RR2.CGmap.gz,WG.CGmap.gz -t RR,RR2,WG -c 5 -C 100 -o matrix.CG.gz

cgmaptools split -i WG.CGmap.gz -p WG -s CGmap.gz

for CHR in 1 2 3 4 5; do (for P in 1 2 3 4 5; do echo | gawk -vC=$CHR -vP=$P -vOFS="\t" '{print "chr"C, P*1000, P*1000+200, "+";}' ; done) ; done > region.bed
zcat WG.CGmap.gz | cgmaptools select region -r region.bed | head


gawk 'NR%100==50' index > site
cgmaptools select site -f RR.CGmap.gz -i site -o RR_select.CGmap.gz

cgmaptools snv -i WG.ATCGmap.gz -m bayes -v bayes.vcf -o bayes.snv --bayes-dynamicP

cgmaptools snv -i WG.ATCGmap.gz -m binom -o binom.snv


cgmaptools dms -i intersect_CG.gz -m 4 -M 100 -o DMS.gz -t fisher
cgmaptools dmr -i intersect_CG.gz -o DMR.gz


gawk '{if(/^#/){print}else{print "chr"$0;}}' bayes.vcf > bayes2.vcf
cgmaptools asm -r genome.fa -b WG.bam -l bayes2.vcf > WG.asm

zcat WG.CGmap.gz | cgmaptools mbed -b region.bed

cgmaptools mbin -i WG.CGmap.gz -B 500 -c 4 -f png -t WG -p WG > mbin.WG.data

cgmaptools mmbin -l WG.CGmap.gz,RR.CGmap.gz,RR2.CGmap.gz,merge.CGmap.gz -c 4 -B 2000 | gawk '{printf("%s:%s-%s", $1, $2, $3); for(i=4;i<=NF;i++){printf("\t%s", $i);} printf("\n");}' > mmbin

cgmaptools heatmap -i mmbin -c -o cluster.pdf -f pdf

for CHR in 1 2 3 4 5; do (for P in 1 2 3 4 5; do echo | gawk -vC=$CHR -vP=$P -vOFS="\t" '{print "chr"C, P*1000, P*1000+1000, "+";}' ; done) ; done | cgmaptools bed2fragreg -n 30 -F 50,50,50,50,50,50,50,50,50,50 -T 50,50,50,50,50,50,50,50,50,50 > fragreg.bed
gunzip -c WG.CGmap.gz | cgmaptools mfg -r fragreg.bed -c 2 -x CG > WG.mfg


cgmaptools mstat -i WG.CGmap.gz -c 4 -f png -p WG -t WG > WG.mstat.data

cgmaptools mtr -i WG.CGmap.gz -r region.bed -o WG.mtr.gz


cgmaptools oac bin -i WG.ATCGmap.gz -B 1000 -f png -p WG -t WG > WG.oac_bin.data

cgmaptools oac stat -i WG.ATCGmap.gz -p WG -f png > WG.oac_stat.data

cgmaptools mec bin -i WG.CGmap.gz -B 1000 -f png -p WG -t WG > WG.mec_bin.data

cgmaptools mec stat -i WG.CGmap.gz -p WG -f png > WG.mec_stat.data

cgmaptools lollipop -i matrix.CG.gz -a anno.refFlat -f pdf




cgmaptools mstat -i WG.CGmap.gz -c 4 -f png -p WG -t WG > WG.mstat.data

cgmaptools fragreg -i WG.mfg -f pdf -o WG_mfg.pdf

cgmaptools tanghulu -r genome.fa -b WG.bam -l chr1:2000-2400 -t CG

cgmaptools findCCGG -i genome.fa -o genome.ccgg


