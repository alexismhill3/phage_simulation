import yaml

def main():
    # initialize yaml dict
    model_params = {"polymerases": [],
                    "ribosomes": [],
                    "species": [],
                    "reactions": []
    }

    # polymerases
    model_params["polymerases"].append({"name": "rnapol-1", "footprint": 35, "speed": 230, "copy_number": 0})
    model_params["polymerases"].append({"name": "rnapol-3.5", "footprint": 35, "speed": 230, "copy_number": 0})
    model_params["polymerases"].append({"name": "ecolipol", "footprint": 35, "speed": 45, "copy_number": 0})
    model_params["polymerases"].append({"name": "ecolipol-p", "footprint": 35, "speed": 45, "copy_number": 0})
    model_params["polymerases"].append({"name": "ecolipol-2", "footprint": 35, "speed": 45, "copy_number": 0})
    model_params["polymerases"].append({"name": "ecolipol-2-p", "footprint": 35, "speed": 45, "copy_number": 0})

    # ribosomes
    model_params["ribosomes"].append({"footprint": 30, "speed": 30, "copy_number": 0})

    # other molecular species
    model_params["species"].append({"name": "bound_ribosome", "copy_number": 10000})
    model_params["species"].append({"name": "bound_ecolipol", "copy_number": 1800})
    model_params["species"].append({"name": "bound_ecolipol_p", "copy_number": 0})
    model_params["species"].append({"name": "ecoli_genome", "copy_number": 0})
    model_params["species"].append({"name": "ecoli_transcript", "copy_number": 0})
    
    # species level reactions
    model_params["reactions"].append({"rate_constant": 1e6, "reactants": ["ecoli_transcript", "__ribosome"], "products": ["bound_ribosome"]})
    model_params["reactions"].append({"rate_constant": 0.04, "reactants": ["bound_ribosome"], "products": ["__ribosome", "ecoli_transcript"]})
    model_params["reactions"].append({"rate_constant": 0.001925, "reactants": ["ecoli_transcript"], "products": ["degraded_transcript"]})
    model_params["reactions"].append({"rate_constant": 1e7, "reactants": ["ecolipol", "ecoli_genome"], "products": ["bound_ecolipol"]})
    model_params["reactions"].append({"rate_constant": 0.3e7, "reactants": ["ecolipol-p", "ecoli_genome"], "products": ["bound_ecolipol_p"]})
    model_params["reactions"].append({"rate_constant": 0.04, "reactants": ["bound_ecolipol"], "products": ["ecolipol", "ecoli_genome", "ecoli_transcript"]})
    model_params["reactions"].append({"rate_constant": 0.04, "reactants": ["bound_ecolipol_p"], "products": ["ecolipol-p", "ecoli_genome", "ecoli_transcript"]})
    model_params["reactions"].append({"rate_constant": 3.8e7, "reactants": ["protein_kinase-0.7", "ecolipol"], "products": ["ecolipol-p", "protein_kinase-0.7"]})
    model_params["reactions"].append({"rate_constant": 3.8e7, "reactants": ["protein_kinase-0.7", "ecolipol-2"], "products": ["ecolipol-2-p", "protein_kinase-0.7"]})
    model_params["reactions"].append({"rate_constant": 3.8e7, "reactants": ["gp-2", "ecolipol"], "products": ["ecolipol-2"]})
    model_params["reactions"].append({"rate_constant": 3.8e7, "reactants": ["gp-2", "ecolipol-p"], "products": ["ecolipol-2-p"]})
    model_params["reactions"].append({"rate_constant": 1.1, "reactants": ["ecolipol-2-p"], "products": ["gp-2", "ecolipol-p"]})
    model_params["reactions"].append({"rate_constant": 1.1, "reactants": ["ecolipol-2"], "products": ["gp-2", "ecolipol"]})
    model_params["reactions"].append({"rate_constant": 3.8e9, "reactants": ["lysozyme-3.5", "rnapol-1"], "products": ["rnapol-3.5"]})
    model_params["reactions"].append({"rate_constant": 3.5, "reactants": ["rnapol-3.5"], "products": ["lysozyme-3.5", "rnapol-1"]})

    model_params["cell_volume"] = 1.1e-15

    with open("../../yaml/t7_model_params_original.yaml", 'w') as file:
        data = yaml.dump(model_params, file)


if __name__ == "__main__":
    main()