#!/bin/bash -l
#SBATCH -J diamond   #jobname
#SBATCH -N 1     #node
#SBATCH --ntasks-per-node=16
#SBATCH --threads-per-core=2
#SBATCH -p bigmem
#SBATCH --mem=60GB

cd /home/pjt6/scratch/cancer_genomes/assemblies


conda activate trinity2.9.1

diamond blastp  --evalue 1e-20 --query AU565_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out AU565_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query MDAMB361_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out MDAMB361_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query AU565_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out AU565_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query MDAMB415_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out MDAMB415_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query HCC_1395_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out HCC_1395_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query MDAMB415_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out MDAMB415_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query HCC_1395_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out HCC_1395_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query MDAMB453_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out MDAMB453_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query HCC_1806_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out HCC_1806_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query MDAMB453_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out MDAMB453_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query HCC_1806_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out HCC_1806_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query SKBR3_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out SKBR3_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query HCC1937_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out HCC1937_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query SKBR3_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out SKBR3_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query HCC1937_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out HCC1937_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query T47D_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out T47D_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query HCC1954_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out HCC1954_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query T47D_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out T47D_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query HCC1954_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out HCC1954_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query ZR751_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out ZR751_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query MCF7_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out MCF7_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query ZR751_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out ZR751_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query MCF7_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out MCF7_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query ZR75B_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out ZR75B_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query MDAMB231_GG.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out MDAMB231_GG_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query ZR75B_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out ZR75B_trinity.Trinity_vs_human.tab -p 16 
 
diamond blastp  --evalue 1e-20 --query MDAMB231_trinity.Trinity.fasta.transdecoder.pep -d h.dmnd --max-target-seqs 1 --outfmt 6 --out MDAMB231_trinity.Trinity_vs_human.tab -p 16 

