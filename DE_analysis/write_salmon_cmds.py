from collections import defaultdict

import os
out = open("salmon_PT.sh","w")
out.write("#$ -cwd\n")

fq = """AU565_paired_1.fq.gz     HCC1937_paired_1.fq.gz  MDAMB231_paired_1.fq.gz  MDAMB453_paired_1.fq.gz  T47D_paired_1.fq
HCC_1395_paired_1.fq.gz  HCC1954_paired_1.fq.gz  MDAMB361_paired_1.fq.gz  SKBR3_paired_1.fq        ZR751_paired_1.fq
HCC_1806_paired_1.fq.gz  MCF7_paired_1.fq.gz     MDAMB415_paired_1.fq.gz  SKBR3_paired_1.fq.gz     ZR75B_paired_1.fq""".split()


trin_path="/shelf/apps/pjt6/apps/trinityrnaseq-Trinity-v2.8.4/util/"
fasta = "nt_human.fa"

# iterate through the file system
count = 0
exp_list = []
for filename in fq:
    count = count + 1
    prefix = filename.split("_paired")[0]
    exp_list.append(prefix)
    cmd = "%s/align_and_estimate_abundance.pl --transcripts %s    --est_method salmon  --seqType fq  --thread_count    16 --output_dir  %s --left    %s_paired_R1.fq.gz    --right %s_paired_R2.fq.gz\n" %(trin_path, fasta, prefix, prefix, prefix)
    if count == 1:
        cmd = cmd.rstrip() + "  --prep_reference \n"
    out.write(cmd)
exp_folder = ""
for i in exp_list:
    temp = i + "/quant.sf "
    exp_folder = exp_folder + temp
est_abun = "%sabundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map  none --name_sample_by_basedir --out_prefix %s.genes.counts.matrix %s \n" % (trin_path, os.path.split(fasta)[-1], exp_folder)
out.write(est_abun)
out.close()
