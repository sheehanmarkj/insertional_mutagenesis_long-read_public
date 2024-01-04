## NOT on the login node, please

seqtk=/path-to-program/seqtk/seqtk
bwa=/path-to-program/bwa-0.7.17/bwa

cp /path/to/host.dna.fa ./host_vector_hybrid.fasta

$seqtk seq /path/to/vector.fa | grep -A1 "^>contig-of-interest1$" >> ./host_vector_hybrid.fasta
$seqtk seq /path/to/vector.fa | grep -A1 "^>contig-of-interest2$" >> ./host_vector_hybrid.fasta
$seqtk seq fxn.fa | grep -A1 "^>contig-of-interest3$" >> ./host_vector_hybrid.fasta

$seqtk seq -l 60 ./host_vector_hybrid.fasta > /scratch/temp; mv /scratch/temp ./host_vector_hybrid.fasta

$bwa index ./host_vector_hybrid.fasta


