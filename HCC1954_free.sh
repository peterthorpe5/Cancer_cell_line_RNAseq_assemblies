#!/bin/bash

#$ -cwd
#$ -V
#$ -N sge3bayes

# The job should be split into a number of sub tasks that is (ideally) at least
# equal to the number of CPUs available on the cluster. The only time you
# wouldn't do this is if the number of contigs to process is less. It's also
# worth putting the number much higher though (eg 500) as the cluster doesn't
# have equal-power CPUs across all nodes; with a higher number, nodes that
# finish their task sooner can then begin work on another segment, rather than
# having all the fast nodes finish their work then sit around doing nothing
# while the slower nodes catch up.

#$ -t 1-10


# The path to Java 7 (or later).
JAVA=/usr/bin/java

# Paths to samtools and freebayes. The --samtools and --freebayes options can
# omitted if these tools already exist on the system or user's path.

SAMTOOLS=/shelf/apps/pjt6/conda/envs/freebayes/bin/samtools
FREEBAYES=/shelf/apps/pjt6/conda/envs/freebayes/bin/freebayes


# Final folder for results. Must exist before the job is started.
OUTPUT=HCC1954


# Define a temporary working directory to write the results to. This specifies
# a folder on the node's scratch space named by jobID-taskID.
WRKDIR=$HOME/$JOB_ID-$SGE_TASK_ID
mkdir $WRKDIR

# Run FreeBayes in parallel.
# The program requires the -b option pointing to the BAM file. All FreeBayes
# options must follow sge3bayes options. Note how the -v output option for
# FreeBayes has the TASK_ID included in its name (so that each part will be
# named 1.vcf, 2.vcf...500.vcf). Remember FreeBayes requires a .fai index file
# and this MUST exist prior to starting the job, otherwise every task will try
# to make the same one and they'll all fail.
# "/usr/local/Modules/modulefiles/tools/freebayes/gitv1_8d2b3a0/bin/freebayes --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.4 " 
$JAVA -cp /shelf/apps/pjt6/apps/sge3bayes/lib/pfb.jar sge3bayes.ParallelRunner --samtools $SAMTOOLS \
--freebayes $FREEBAYES -b HCC1954*.bam -f Homo_sapiens.GRCh38.dna.primary_assembly.fa -v $WRKDIR/$SGE_TASK_ID.vcf \
 --no-population-priors --min-alternate-count 5 

# Once finished, copy the task's .vcf file back to the output folder.
mv $WRKDIR/* $OUTPUT/
rm -rf $WRKDIR