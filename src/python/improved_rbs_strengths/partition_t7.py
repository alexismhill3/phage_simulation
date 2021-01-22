# partions the T7 genome into ~10000 bp segments, with 100 bp overlap between segments
# preps input for the Salis RBS calculator (which only accepts sequences <10000 bp in length)

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

t7_record = SeqIO.read("../../../data/t7_genome_full.fasta", "fasta")
t7_sequence = t7_record.seq

lower_bound = 0
upper_bound = 0
genome_partitions = []
while upper_bound != len(t7_sequence):
    upper_bound = min(lower_bound + 10000, len(t7_sequence))
    seq_partition = t7_sequence[lower_bound:upper_bound]
    genome_partitions.append(SeqRecord(seq_partition, 
                                       id = f"{t7_record.id}-{lower_bound}-{upper_bound - 1}",
                                       name = "",
                                       description = f"{t7_record.description}, position {lower_bound}..{upper_bound - 1}")
                            )
    lower_bound += 9900

SeqIO.write(genome_partitions, "t7_genome_split.fasta", "fasta")