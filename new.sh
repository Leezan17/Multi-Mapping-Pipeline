#!/bin/bash
#SBATCH -J jellyfish_map
#SBATCH -p charlesworths
#SBATCH -A alarracu_lab
#SBATCH -t 05-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 10
#SBATCH --mem=24gb
#SBATCH --mail-user=tpham9@u.rochester.edu
#SBATCH --mail-type=ALL
#SBATCH --error=mapp.out
#SBATCH -o mapp.log

module load bedtools
module load jellyfish
module load python3
gff=$1
reference=$2
kmer=$3
#thread=$4
dir=$(pwd)
index=dmel_iso1.jf
mkdir parse_gff
jellyfish count -m $kmer -s 100M -t 10 $reference -o $index
awk 'BEGIN {counter=0}{printf "%s\n", $n>"parse_gff/"counter".gff"} {counter++}' $gff
#bedtools getfasta -fi /scratch/alarracu_lab/Genome/Dmelanogaster/dmel_scaffold2_plus0310.fasta -bed locus.gff >all.fasta
cd parse_gff
for fullname in *.gff; do
        name=${fullname%.gff}
        echo $name
        fasta=$name.fasta
        jf=$name.jf
        bedtools getfasta -fi /scratch/alarracu_lab/Genome/Dmelanogaster/dmel_scaffold2_plus0310.fasta -bed $fullname >$fasta
        echo $jf
        jellyfish count -m 31 -s 100M -t 10 $fasta -o $jf
        local=$name.count.1
        global=$name.count.2
        jellyfish query $jf -s $fasta >$local
        echo $fasta
        jellyfish query $dir/$index -s $fasta> $global
        python3 ../graph.py $local $global $fullname
done

