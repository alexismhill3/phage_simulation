import pinetree as pt
import phage_differential_evolution
import read_from_config
import tempfile

def run_simulation(x, *args):
    model_config, diff_evolution_config = *args

    # initalize genome from fixed paramaters
    genome = read_from_config.build_genome(model_config)
    # add variable parameters to fully specify genome
    genome_with_params = phage_differential_evolution.specify_model(genome, x, diff_evolution_config)
    # register genome and finishing specifying the simulation
    model = read_from_config.build_model(genome_with_params, model_config)

    # create a temporary file to write output to
    out = tempfile.NamedTemporaryFile()
    model.simulate(time_limit=diff_evolution_config["time-limit"], time_step=diff_evolution_config["time-step", output=out.name])
    results = phage_differential_evolution.parse_output(out.name)
    
    return phage_differential_evolution.calculate_rmse(results)