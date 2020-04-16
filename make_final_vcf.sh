#!/bin/bash

#$ -cwd
#$ -V
#$ -N sge3bayes


# The path to Java 7 (or later).
JAVA=/shelf/apps/pjt6/conda/envs/trinity/bin/java

# The command takes two (required parameters):
# -d --dir      path to the directory containing vcf files

cd /storage/home/users/pjt6/cancer_cell_lines/bam_files

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir AU565 --output AU565.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir HCC --output HCC.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir index.sh --output index.sh.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir MDAMB361 --output MDAMB361.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir nem --output nem.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir ZR751 --output ZR751.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir example-concat.sh --output example-concat.sh.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir HCC1937 --output HCC1937.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir MCF7 --output MCF7.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir MDAMB415 --output MDAMB415.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir SKBR3 --output SKBR3.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir ZR75B --output ZR75B.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir HCC --output HCC.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir HCC1954 --output HCC1954.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir MDAMB231 --output MDAMB231.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir MDAMB453 --output MDAMB453.vcf

$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.Concat --dir T47D --output T47D.vcf

