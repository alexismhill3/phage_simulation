import load_sim_from_config

conf = "../../yaml/t7_wt_degrade_jack_2019.yaml"
model = load_sim_from_config.construct_model(conf)
model.simulate(time_limit=1200, time_step=5, output="phage_degrade_test.tsv")