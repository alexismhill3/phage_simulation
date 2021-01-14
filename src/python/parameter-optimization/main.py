from scipy.optimize import differential_evolution
import minimize
import sys

model_config, diff_evolution_config = sys.argv
bounds = diff_evolution_config["bounds"]
result = differential_evolution(minimize.run_simulation, bounds, args=(model_config, diff_evolution_config))
print(result.x)