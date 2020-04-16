#!/bin/bash

#$ -cwd
#$ -V
#$ -N sge3bayes


# The path to Java 7 (or later).
JAVA=/shelf/apps/pjt6/conda/envs/trinity/bin/java

# The command takes two (required parameters):
# -d --dir      path to the directory containing vcf files
# -o --output   file to create containing the concatenated results
$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir output --output Bedale.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir Lindley --output Lindley.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir newt --output Newton.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir output_pa1 --output pa1.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir pa4 --output pa4.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir output_pa5 --output pa5.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir ../output --output Luffness.vcf

sed -i 's/unknown/Bedale/g' Bedale.vcf 
sed -i 's/unknown/Lindley/g' Lindley.vcf 
sed -i 's/unknown/Luffness/g' Luffness.vcf 
sed -i 's/unknown/Newton/g' Newton.vcf 
sed -i 's/unknown/pa1/g' pa1.vcf 
sed -i 's/unknown/pa4/g' pa4.vcf 
sed -i 's/unknown/pa5/g' pa5.vcf 


bgzip -c Bedale.vcf > Bedale.vcf.gz
tabix -p vcf Bedale.vcf.gz

bgzip -c Lindley.vcf > Lindley.vcf.gz
tabix -p vcf Lindley.vcf.gz

bgzip -c Luffness.vcf > Luffness.vcf.gz
tabix -p vcf Luffness.vcf.gz

bgzip -c Newton.vcf > Newton.vcf.gz
tabix -p vcf Newton.vcf.gz

bgzip -c pa1.vcf > pa1.vcf.gz
tabix -p vcf pa1.vcf.gz

bgzip -c pa4.vcf > pa4.vcf.gz
tabix -p vcf pa4.vcf.gz

bgzip -c pa5.vcf > pa5.vcf.gz
tabix -p vcf pa5.vcf.gz


/storage/home/users/pjt6/project/nem_pops/merged/vcftools/bin/vcf-merge Bedale.vcf.gz Lindley.vcf.gz Luffness.vcf.gz Newton.vcf.gz pa1.vcf.gz pa4.vcf.gz pa5.vcf.gz > all_nem_pops.vcf
# het

vcftools --vcf all_nem_pops.vcf --het --out all_nem_pops.het

# caluclate the stats on the vcf file

/storage/home/users/pjt6/rtg-tools-3.10.1/rtg vcfstats --allele-lengths all_nem_pops.vcf > all_nem_pops.vcf_SNP_stats.txt

