# Constructs a yaml representation of the wild type T7 genome from the corresponding genbank file

from Bio import Entrez, SeqIO
import yaml


PHI10_BIND = 1.82e8  # Binding constant for phi10

IGNORE_REGULATORY = ["E. coli promoter E[6]",
                     "T7 promoter phiOR",
                     "T7 promoter phiOL",
                     "E. coli promoter A0 (leftward)"]

RELABEL_GENES = {"gene 2": "gp-2",
                 "gene 1": "rnapol-1",
                 "gene 3.5": "lysozyme-3.5",
                 "gene 0.7": "protein_kinase-0.7"}


def get_promoter_interactions(name):
    '''
    Calculate promoter binding strengths. The relative strengths defined here
    come from 2012 Covert, et al paper.
    '''
    ecoli_strong = ["E. coli promoter A1",
                    "E. coli promoter A2",
                    "E. coli promoter A3"]
    ecoli_weak = ["E. coli B promoter",
                  "E. coli C promoter"]
    phi1_3 = ["T7 promoter phi1.1A",
              "T7 promoter phi1.1B",
              "T7 promoter phi1.3",
              "T7 promoter phi1.5",
              "T7 promoter phi1.6"]
    phi3_8 = ["T7 promoter phi2.5",
              "T7 promoter phi3.8",
              "T7 promoter phi4c",
              "T7 promoter phi4.3",
              "T7 promoter phi4.7"]
    phi6_5 = ["T7 promoter phi6.5"]
    phi9 = ["T7 promoter phi9"]
    phi10 = ["T7 promoter phi10"]
    phi13 = ["T7 promoter phi13",
             "T7 promoter phi17"]

    if name in ecoli_strong:
        return {'ecolipol': 10e4,
                'ecolipol-p': 3e4}
    elif name in ecoli_weak:
        return {'ecolipol': 1e4,
                'ecolipol-p': 0.3e4}
    elif name in phi1_3:
        return {'rnapol-1': PHI10_BIND * 0.01,
                'rnapol-3.5': PHI10_BIND * 0.01 * 0.5}
    elif name in phi3_8:
        return {'rnapol-1': PHI10_BIND * 0.005,
                'rnapol-3.5': PHI10_BIND * 0.005 * 0.5}
    elif name in phi6_5:
        return {'rnapol-1': PHI10_BIND * 0.2,
                'rnapol-3.5': PHI10_BIND * 0.2}
    elif name in phi9:
        return {'rnapol-1': PHI10_BIND * 0.2,
                'rnapol-3.5': PHI10_BIND * 0.2}
    elif name in phi10:
        return {'rnapol-1': PHI10_BIND,
                'rnapol-3.5': PHI10_BIND}
    elif name in phi13:
        return {'rnapol-1': PHI10_BIND * 0.05,
                'rnapol-3.5': PHI10_BIND * 0.05}
    else:
        raise ValueError(
            "Promoter strength for {0} not assigned.".format(name))


def get_terminator_interactions(name):
    '''
    Get terminator efficiencies.
    '''
    if name == "E. coli transcription terminator TE":
        return {'ecolipol': 1.0,
                'ecolipol-p': 1.0,
                'rnapol-1': 0.0,
                'rnapol-3.5': 0.0}
    elif name == "T7 transcription terminator Tphi":
        return {'rnapol-1': 0.85,
                'rnapol-3.5': 0.85}
    else:
        return {'name': 0.0}


def main():
    # initialize yaml dict
    phage_params = {"gene": [],
                    "promoter": [],
                    "terminator": [],
                    "rnase_site": []
    }

    # Download T7 wild-type genbank records
    Entrez.email = "alexismhill3@gmail.com"
    handle = Entrez.efetch(db="nuccore",
                           id=["NC_001604"],
                           rettype="gb",
                           retmode="text")

    record = SeqIO.read(handle, "genbank")
    phage_params["genome_length"] = len(record.seq)

    for feature in record.features:
        start = feature.location.start.position + 1
        stop = feature.location.end.position
        name = ''
        if "note" in feature.qualifiers:
            name = feature.qualifiers["note"][0]
        # Grab promoters and terminators
        if feature.type == "regulatory":
            if name in IGNORE_REGULATORY:
                continue
            # Construct promoter
            if "promoter" in feature.qualifiers["regulatory_class"]:
                interactions = get_promoter_interactions(name)
                phage_params["promoter"].append({"name": name, "start": start, "stop": stop, "interactions": interactions})
            # Construct terminator params
            if "terminator" in feature.qualifiers["regulatory_class"]:
                interactions = get_terminator_interactions(name)
                phage_params["terminator"].append({"name": name, "start": start, "stop": stop, "interactions": interactions})
        # Grab genes/CDSes
        if feature.type == "gene":
            #if name in IGNORE_GENES:
                #continue
            if name in RELABEL_GENES:
                name = RELABEL_GENES[name]
            # Construct CDS parameters for this gene
            # default rbs stregth is from Jack et al., 2019
            phage_params["gene"].append({"name": name, "start": start, "stop": stop, "rbs": -30, "rbs_strength": 1e7})
        # Rnase sites
        if feature.type == "misc_structure":
            # default rnase strength is from Jack et al., 2019
            phage_params["rnase_site"].append({"name": name, "start": start, "stop": stop, "rnase_strength": 1e-2})

    with open("../../yaml/t7_wt_original.yaml", 'w') as file:
        data = yaml.dump(phage_params, file)


if __name__ == "__main__":
    main()
