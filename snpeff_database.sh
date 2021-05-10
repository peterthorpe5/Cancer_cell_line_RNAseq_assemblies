# How to creat a snpeff database! The online doc sucks!
filter snps:http://www.ddocent.com/filtering/


step 1)


# this has a lot of the dependencies sorted. 
conda activate snpeffwrapper


cd /storage/home/users/pjt6/shelf_apps/apps/snpEff
# edit the config file
nano snpEffect.config 

# A.pisum N116 nmodels
a.pisum.N116.models.20180219.genome : A.pisum N116 PT

# A.pisum PSO1 nmodels
a.pisum.PSO1.models.20180219.genome : A.pisum PSO1 PT

e.g.

DZA.genome :  DZA


mkdir DZA (in the data folder)

snpEff.N116.genome <- "a.pisum.N116.models.20180219"
snpEff.PSO1.genome <- "a.pisum.PSO1.models.20180219"


# make a folder in the data folder in the snpeff software folder, 
mkdir /snpEff/data/a.pisum.N116.models.20180219
cp /genome/A.pisum_N116.genome.v1.0.fasta /snpEff/data/a.pisum.N116.models.20180219/sequences.fa
cp /genome/A.pisum_N116.models.20180219.gff /snpEff/data/a.pisum.N116.models.20180219/genes.gff
cp /genome/A.pisum_N116.models.20180219.aa.fasta /snpEff/data/a.pisum.N116.models.20180219/proteins.fa

java -jar ./snpEff.jar build -gff3 a.pisum.N116.models.20180219


#

copy all the relevent files to DZA.

then 

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx4g -jar snpEff.jar build -gff3 -v DZA

########################################################################################################
# example 2:

GCA_003473485.2_MtrunA17r5.0-ANR_genomic.genome : GCA_003473485.2_MtrunA17r5.0-ANR_genomic

cd data
mkdir GCA_003473485.2_MtrunA17r5.0-ANR_genomic

cd GCA_003473485.2_MtrunA17r5.0-ANR_genomic

cp /storage/home/users/pjt6/medicago/A17/GCA_003473485.2_MtrunA17r5.0-ANR_genomic.fna sequences.fa
cp /storage/home/users/pjt6/medicago/A17/GCA_003473485.2_MtrunA17r5.0-ANR_genomic.gff genes.gff
cp /storage/home/users/pjt6/medicago/A17/GCA_003473485.2_MtrunA17r5.0-ANR_protein.faa proteins.fa

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx4g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar build -gff3 -v GCA_003473485.2_MtrunA17r5.0-ANR_genomic


# first filter the snps a bit:
vcftools --vcf Medicago_V5_SNPs_freebayes.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --minDP 3  --out Medicago_V5_SNPs_freebayes.raw.g5mac3dp3.vcf
http://www.ddocent.com/filtering/

to run:

#java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar  Acyrthosiphon_pisum raw.g5mac3.recode.vcf > ACYPI_gene_PS01_N116_g5mac3.recode.snpef.vcf

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GCA_003473485.2_MtrunA17r5.0-ANR_genomic Medicago_V5_SNPs_freebayes.raw.g5mac3dp3.vcf.recode.vcf >  Medicago_V5_SNPs_freebayes.raw.g5mac3dp3.snp_effect.vcf






# first filter the snps a bit:
vcftools --vcf Medicago_V4_SNPs_freebayes.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --minDP 3  --out Medicago_V4_SNPs_freebayes.raw.g5mac3dp3.vcf
http://www.ddocent.com/filtering/

to run:

java -XX:-UsePerfData -Djava.io.tmpdir=/storage/home/users/pjt6/ -Xmx30g -jar /storage/home/users/pjt6/shelf_apps/apps/snpEff/snpEff.jar GCA_000219495.2_MedtrA17_4.0_genomic Medicago_V4_SNPs_freebayes.raw.g5mac3dp3.vcf.recode.vcf >  Medicago_V4_SNPs_freebayes.raw.g5mac3dp3.snp_effect.vcf




example 3:
I had to make the data folder this time. 
newtonv1.0

/mnt/shared/scratch/kv42609/common_variant_calling/SnpEff/data/newtonv1.0

cp Gp_Newton_haplotype1.fasta /mnt/shared/scratch/kv42609/common_variant_calling/SnpEff/data/newtonv1.0/sequences.fa
cp Gp_Newton_haplotype1.fasta /mnt/shared/scratch/kv42609/common_variant_calling/SnpEff/data/newtonv1.0/Gp_Newton_haplotype1.fasta

cp Gpal_newton_newton.gff3 /mnt/shared/scratch/kv42609/common_variant_calling/SnpEff/data/newtonv1.0/genes.gff

cp Gpal_newton_newton.proteins.fa  /mnt/shared/scratch/kv42609/common_variant_calling/SnpEff/data/newtonv1.0/proteins.fa


java -jar /snpEff/snpEff.jar build -gff3 newtonv1.0

