import load_sim_from_config

t7_conf = "../../yaml/t7_wt_degrade_jack_2019.yaml"
plasmid_conf = "../../yaml/three_genes_plasmid.yaml"
phage = load_sim_from_config.construct_genome(t7_conf)
plasmid = load_sim_from_config.construct_genome(plasmid_conf)
genomes = phage + plasmid
model = load_sim_from_config.construct_model_with_genomes(t7_conf, genomes)
model.simulate(time_limit=1000, time_step=5, output="two_different_genomes.tsv")