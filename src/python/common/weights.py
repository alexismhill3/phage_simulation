from Bio import Entrez, SeqIO

OPT_CODONS_E_COLI = {'A': ['GCT'],
                     'R': ['CGT', 'CGC'],
                     'N': ['AAC'],
                     'D': ['GAC'],
                     'C': ['TGC'],
                     'Q': ['CAG'],
                     'E': ['GAA'],
                     'G': ['GGT', 'GGC'],
                     'H': ['CAC'],
                     'I': ['ATC'],
                     'L': ['CTG'],
                     'F': ['TTC'],
                     'P': ['CCG'],
                     'S': ['TCT', 'TCC'],
                     'T': ['ACT', 'ACC'],
                     'Y': ['TAC'],
                     'V': ['GTT', 'GTA']}

def compute_cds_weights(record, feature, factor, weights):
    # Grab the gene name
    nuc_seq = feature.location.extract(record).seq
    aa_seq = feature.qualifiers["translation"][0]
    weight_sum = 0
    for index, nuc in enumerate(nuc_seq):
        aa_index = int(index / 3)
        codon_start = aa_index * 3
        codon = nuc_seq[codon_start:codon_start + 3]
        genome_index = feature.location.start + index
        if aa_index < len(aa_seq):
            if aa_seq[aa_index] in OPT_CODONS_E_COLI:
                if codon in OPT_CODONS_E_COLI[aa_seq[aa_index]]:
                    weights[genome_index] = factor
                    weight_sum += factor
                else:
                    weights[genome_index] = 1
                    weight_sum += 1
    return weights


def normalize_weights(weights):
    # Average over all CDSs, which will have non-zero weights
    non_zero = sum(1 if i != 0 else 0.0 for i in weights)
    mean_weight = sum(weights) / non_zero
    norm_weights = [i / mean_weight for i in weights]
    # Replace non-CDS weights with 1
    norm_weights = [1 if i == 0 else i for i in norm_weights]
    return norm_weights


def t7_weights(factor):
    # Download T7 wild-type genbank records
    Entrez.email = "alexis.hill@gmail.com"
    handle = Entrez.efetch(db="nuccore",
                           id=["NC_001604"],
                           rettype="gb",
                           retmode="text")

    record = SeqIO.read(handle, "genbank")
    weights = [0.0] * len(record.seq)
    for feature in record.features:
        if feature.type == "CDS":
            weights = compute_cds_weights(record, feature, factor, weights)
    
    norm_weights = normalize_weights(weights)
    return norm_weights