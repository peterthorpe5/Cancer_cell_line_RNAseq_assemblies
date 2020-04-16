#!/bin/bash

#$ -cwd
#$ -V
#$ -N sge3bayes

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  AU565.vcf > AU565_Snpeff.vcf

mv snpEff_genes.txt AU565_snpEff_genes.txt

mv snpEff_summary.html AU565_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  HCC.vcf > HCC_Snpeff.vcf

mv snpEff_genes.txt HCC_snpEff_genes.txt

mv snpEff_summary.html HCC_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  index.sh.vcf > index.sh_Snpeff.vcf

mv snpEff_genes.txt index.sh_snpEff_genes.txt

mv snpEff_summary.html index.sh_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  MDAMB361.vcf > MDAMB361_Snpeff.vcf

mv snpEff_genes.txt MDAMB361_snpEff_genes.txt

mv snpEff_summary.html MDAMB361_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  nem.vcf > nem_Snpeff.vcf

mv snpEff_genes.txt nem_snpEff_genes.txt

mv snpEff_summary.html nem_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  ZR751.vcf > ZR751_Snpeff.vcf

mv snpEff_genes.txt ZR751_snpEff_genes.txt

mv snpEff_summary.html ZR751_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  example-concat.sh.vcf > example-concat.sh_Snpeff.vcf

mv snpEff_genes.txt example-concat.sh_snpEff_genes.txt

mv snpEff_summary.html example-concat.sh_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  HCC1937.vcf > HCC1937_Snpeff.vcf

mv snpEff_genes.txt HCC1937_snpEff_genes.txt

mv snpEff_summary.html HCC1937_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  MCF7.vcf > MCF7_Snpeff.vcf

mv snpEff_genes.txt MCF7_snpEff_genes.txt

mv snpEff_summary.html MCF7_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  MDAMB415.vcf > MDAMB415_Snpeff.vcf

mv snpEff_genes.txt MDAMB415_snpEff_genes.txt

mv snpEff_summary.html MDAMB415_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  SKBR3.vcf > SKBR3_Snpeff.vcf

mv snpEff_genes.txt SKBR3_snpEff_genes.txt

mv snpEff_summary.html SKBR3_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  ZR75B.vcf > ZR75B_Snpeff.vcf

mv snpEff_genes.txt ZR75B_snpEff_genes.txt

mv snpEff_summary.html ZR75B_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  HCC.vcf > HCC_Snpeff.vcf

mv snpEff_genes.txt HCC_snpEff_genes.txt

mv snpEff_summary.html HCC_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  HCC1954.vcf > HCC1954_Snpeff.vcf

mv snpEff_genes.txt HCC1954_snpEff_genes.txt

mv snpEff_summary.html HCC1954_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  MDAMB231.vcf > MDAMB231_Snpeff.vcf

mv snpEff_genes.txt MDAMB231_snpEff_genes.txt

mv snpEff_summary.html MDAMB231_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  MDAMB453.vcf > MDAMB453_Snpeff.vcf

mv snpEff_genes.txt MDAMB453_snpEff_genes.txt

mv snpEff_summary.html MDAMB453_snpEff_summary.html

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GRCh38.p7.RefSeq  T47D.vcf > T47D_Snpeff.vcf

mv snpEff_genes.txt T47D_snpEff_genes.txt

mv snpEff_summary.html T47D_snpEff_summary.html
