# main.py
#
# Simulating ninth planet and perform analysis on stability
#


# Import functions
import functions as func
import matplotlib.pyplot as plt
import rebound
import numpy as np


if __name__ == "__main__":
    
    # sim_start = rebound.Simulation("simulations/" + "ss_1990.bin")
    # ps_start = sim_start.particles
    # func.AddPlanet(sim_start, a=0.55, e=0, inc=0)
    # ps = sim_start.particles
    # r9 = ps[9]
    #func.Integrate(sim_start, f"Integration_1990_dt0.2_t50000_e{r9.e}_inc{r9.inc}_a{r9.a}", years=50000, step=10000)  
    
    #df = func.Unload("Data")
    #func.AnimatePositions(df, filename="Orbit_e3")

    #func.AnalyseConvergence()

    func.AnalyseStability()
    
    #func.AnalyseAccuracyErrors()

    #func.AnalyseSolarSystemPerturbation()

    #func.AnalyseArchive()


