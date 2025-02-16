# functions.py
#
# Contains all function regarding rebound simulation and analysis of ninth planet
#


# Imports
import rebound
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import numpy as np
import pandas as pd
from datetime import datetime


# Global variables
SIMULATIONS_FOLDER = "simulations/"
PLOTS_FOLDER = "plots/"
OUTPUT_FOLDER = "output/"
ANIMATIONS_FOLDER = "animations/"
LABELS = {
    0: "Sun",
    1: "Mercury",
    2: "Venus",
    3: "Earth",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "R9"
}


def InitSolarSytem():
    """Returns rebound simulation object that contains all planets"""

    print("SolarSystem...")

    # File and directory of solar system simulation object
    ss_file = "ss.bin"

    # Check if solarsystem simulation object is already in simulations folder
    if os.path.exists(SIMULATIONS_FOLDER + ss_file):

        sim = rebound.Simulation(SIMULATIONS_FOLDER + ss_file)
        return sim

    # Check if saving directory exists otherwise make one
    if not os.path.isdir(SIMULATIONS_FOLDER):

        os.mkdir(SIMULATIONS_FOLDER)
        print(f"NOTE: created this folder {SIMULATIONS_FOLDER}")
    
    # Create simulation object
    sim = rebound.Simulation()

    # Adding planets in our solar system
    ss_objects = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
    sim.add(ss_objects)

    # Save simulation object to file
    sim.save(SIMULATIONS_FOLDER + ss_file)

    return sim


def PlotOrbits(sim, filename="auto", show=False, overwrite=False):
    """Show orbits of all objects in rebound simulation object"""

    print("PlotOrbit...")

    # Create figure
    fig, ax1, ax2, ax3 = rebound.OrbitPlot(sim, 
        unitlabel="[AU]",  
        color=True, 
        slices=0.5, 
        xlim=[-10,10], 
        ylim=[-10,10])

    # If saving directory does not exists make one
    if not os.path.isdir(PLOTS_FOLDER):

        os.mkdir(PLOTS_FOLDER)
        print(f"NOTE: created this folder {PLOTS_FOLDER}")
    
    # Check if filename is already used
    if os.path.exists(PLOTS_FOLDER + filename + ".png") and not overwrite:

        print(f"WARNING: The file {filename}.png already exists, auto filename has been used")
        filename = "auto"

    # Create automatic filename
    if filename == "auto":

        # Get string if current date and time
        current_time = datetime.now().strftime("%Y-%m-%d_%H-%M")

        filename = "Orbits_" + current_time

    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False) 

    # Save figure
    fig.savefig(PLOTS_FOLDER + filename + ".png", dpi=200) 

    # Show figure
    if show:
        fig.show()

    # Clearing figure window for next figures
    fig.clf()


def PlotOrbits3D(sim, filename="auto", show=False, overwrite=False):
    """Show 3D orbits of all objects in rebound simulation object"""

    print("PlotOrbits3D...")

    # Setup plot variables
    fig = plt.figure(figsize=(10,10), dpi=500)
    ax = fig.gca(projection='3d')
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    ax.set_zlabel('z [AU] ')
    plt.title('Solar system with planet R9')

    # Plot particles
    ps = sim.particles
    ax.scatter(ps[0].x, ps[0].y, ps[0].z, color='k', marker='*')
    ax.scatter(ps[1].x, ps[1].y, ps[1].z, color='k', marker='o')
    ax.scatter(ps[2].x, ps[2].y, ps[2].z, color='k', marker='o')
    ax.scatter(ps[3].x, ps[3].y, ps[3].z, color='k', marker='o')
    ax.scatter(ps[4].x, ps[4].y, ps[4].z, color='k', marker='o')
    ax.scatter(ps[5].x, ps[5].y, ps[5].z, color='k', marker='o')
    ax.scatter(ps[6].x, ps[6].y, ps[6].z, color='k', marker='o')
    ax.scatter(ps[7].x, ps[7].y, ps[7].z, color='k', marker='o')
    ax.scatter(ps[8].x, ps[8].y, ps[8].z, color='k', marker='o')
    ax.scatter(ps[9].x, ps[9].y, ps[9].z, color='r', marker='o')

    # Loop through all planets
    for planet in ps[1:]:

        # Get orbit using kepler
        o = np.array(planet.sample_orbit())
        
        # Save orbit positions
        pos_x, pos_y, pos_z = [], [], []
        for x, y, z in o:
            pos_x.append(x)
            pos_y.append(y)
            pos_z.append(z)

        # Plot orbit
        plt.plot(pos_x, pos_y, pos_z, 'k-', alpha=0.2)

    # If saving directory does not exists make one
    if not os.path.isdir(PLOTS_FOLDER):

        os.mkdir(PLOTS_FOLDER)
        print(f"NOTE: created this folder {PLOTS_FOLDER}")
    
    # Check if filename is already used
    if os.path.exists(PLOTS_FOLDER + filename + ".png") and not overwrite:

        print(f"WARNING: The file {filename}.png already exists, auto filename has been used")
        filename = "auto"

    # Create automatic filename
    if filename == "auto":

        # Get string if current date and time
        current_time = datetime.now().strftime("%Y-%m-%d_%H-%M")

        filename = "Orbits3D_" + current_time

    # Save figure
    fig.savefig(PLOTS_FOLDER + filename + ".png") 

    # Show figure
    if show:
        fig.show()

    # Clearing figure window for next figures
    fig.clf()


def AnalyseStability():
    """Analysing difference between two integrations"""

    print("AnalyseStability...")

    # Open simulation archives
    sa_ss = rebound.SimulationArchive(f"{SIMULATIONS_FOLDER}Integration_1990_dt0.2_t50000.bin")
    sa_r9 = rebound.SimulationArchive(f"{SIMULATIONS_FOLDER}Integration_1990_dt0.2_t50000_e0_inc0_a0.55.bin")

    # Check if simulation archives have the same length
    if len(sa_ss) != len(sa_r9):
        print("WARNING: The simulation archives have not the same length")
        return

    # Initialize array
    time = []
    v_dif = [[],[],[],[],[],[],[],[],[],[]]
    inc_dif = [[],[],[],[],[],[],[],[],[],[]]
    e_dif = [[],[],[],[],[],[],[],[],[],[]]
    a_with_R9 = [[],[],[],[],[],[],[],[],[],[]]
    sim_start = sa_r9[0]
    ps_start = sim_start.particles

    # Iterate over two simulations
    for sim_ss, sim_r9 in zip(sa_ss, sa_r9):
    
        # Save time
        t = sim_ss.t /(2*np.pi)
        time.append(t)

        # Get list of particles without sun
        ps_ss = sim_ss.particles[1:9]
        ps_r9 = sim_r9.particles[1:9]

        # Set particle index
        i = 1
        # Loop over particles
        for p_ss, p_r9 in zip(ps_ss, ps_r9):

            # Get eccentricity difference
            p_inc_norm = p_ss.inc
            p_inc_new = p_r9.inc
            p_inc_dif = p_inc_norm - p_inc_new
            p_inc_dif = abs(p_inc_new - ps_start[i].inc)
            a_r9 = p_r9.a

            # Get speed difference
            p_v_norm = np.sqrt(p_ss.vx**2 + p_ss.vy**2 + p_ss.vz**2)
            p_v_new = np.sqrt(p_r9.vx**2 + p_r9.vy**2 + p_r9.vz**2)
            p_v_dif = p_v_norm - p_v_new

            # Get eccentricity
            p_e_norm = p_ss.e
            p_e_new = p_r9.e
            p_e_dif = p_e_norm - p_e_new
            p_e_dif = abs(p_e_new - ps_start[i].e)
            a_a0_dif = abs(a_r9 - ps.start[i].a)

            # Save information
            inc_dif[i].append(p_inc_dif)
            v_dif[i].append(p_v_dif)
            e_dif[i].append(p_e_dif)
            a_with_R9.append(a_a0_dif)
            
            # Increment particle index
            i += 1
    
    # Get max perturbation
    for i in range(1, 9):
        
        print("-----------------------------")
        print(LABELS[i])
        print(f"max_inc = {max(inc_dif[i])}")
        print(f"max_e = {max(e_dif[i])}")

    plt.title("Semimajor axis evolution")
    for i in range(1,9):
         plt.plot(time, a_with_R9[i], label=LABELS[i], alpha=.9, lw=.8)  
    plt.legend(loc='upper right', prop={'size': 6})
    plt.xlabel("time [years]")
    plt.ylabel("$|inc - inc_0|$")
    #plt.yscale("log")
    plt.show()
    # # Subplot
    # plt.subplot(2,1,1)
    # plt.title("dt0.2 t50000 e0 inc0 a0.55")
    # for i in range(1,9):
    #     plt.plot(time, inc_dif[i], label=LABELS[i], alpha=.9, lw=.8)    
    # plt.legend(loc='upper right', prop={'size': 6})
    # plt.yscale("log")
    # plt.xlabel("time [years]")
    # plt.ylabel("$|inc - inc_0|$")

    # # Subplot
    # plt.subplot(2,1,2)
    # # plt.title("Eccentricity difference")
    # for i in range(1,9):
    #     plt.plot(time, e_dif[i], label=LABELS[i], alpha=.9, lw=.8)    
    # plt.legend(loc='upper right', prop={'size': 6})
    # plt.yscale("log")
    # plt.xlabel("time [years]")
    # plt.ylabel("$|e - e_0|$")

    # # Subplot
    # # plt.subplot(3,1,3)
    # # # plt.title("Speed difference")
    # # for i in range(1,8):
    # #     plt.plot(time, v_dif[i], label=LABELS[i])      
    # # plt.legend(loc='upper right', prop={'size': 6})
    # # plt.xlabel("time [years]")
    # # plt.ylabel("d|v| [m/s]")

    # # Save figure
    # plt.tight_layout()
    # filename = "Stability_t50000"
    # plt.savefig(f"{PLOTS_FOLDER}{filename}.png", dpi=200)


def AddPlanet(sim, a=0.55, e=0, inc=0):
    """Add planet to simulation"""

    print("AddPlanet...")

    sim.add(m=3e-5, a=a, e=e, inc=inc)


def Integrate(sim, filename="auto", overwrite=False, dt=.2, years=10, step=2):
    """Integrate all objects in simulation, dt in days and t_end in years"""

    print("Integrate...")

    # Integration settings
    t_start =  0
    t_end = 2*np.pi * years
    N = int((years * 365) / dt)
    sim.dt = t_end / N 
    sim.integrator = "leapfrog"
    print(f"  dt = {dt} days")
    print(f"  t_end = {years} years")
    print(f"  N = {N} steps") 
    print(f"  sim.dt = {sim.dt}")
    estimated_time = N**1.1 / 175 * 2 / step 
    seconds = int(estimated_time % 60)
    minutes = int(estimated_time / 60)
    print(f"  estimated running time = {minutes} min and {seconds} s")

    # If saving directory does not exists make one
    if not os.path.isdir(OUTPUT_FOLDER):

        os.mkdir(OUTPUT_FOLDER)
        print(f"NOTE: created this folder {OUTPUT_FOLDER}") 

    # Check if filename is already used
    if os.path.exists(OUTPUT_FOLDER + filename + ".csv") and not overwrite:

        print(f"WARNING: The file {filename}.csv already exists, auto filename has been used")
        filename = "auto"

    # Create automatic filename
    if filename == "auto":

        ps = sim.particles

        # Check if a nine planet is added
        if len(ps) < 10:
            filename = f"Integration_dt{dt}_t{years}"

        else:
            r9 = ps[9]
            filename = f"Integration_dt{dt}_t{years}_e{r9.e}_inc{r9.inc}_a{r9.a}"

    # Save simulation objects in archive, this might be usefull later on
    sim.automateSimulationArchive(f"{SIMULATIONS_FOLDER}{filename}.bin", step=step , deletefile=True)

    # Setup structure for saving data per timestamp
    # data_dict = {}
    # data_dict["time"] = np.array([])

    # ps = sim.particles
    # for i in range(len(ps)):

    #     data_dict[f"{i}_x"] = np.array([])
    #     data_dict[f"{i}_y"] = np.array([])
    #     data_dict[f"{i}_z"] = np.array([])
    #     data_dict[f"{i}_vx"] = np.array([])
    #     data_dict[f"{i}_vy"] = np.array([])
    #     data_dict[f"{i}_vz"] = np.array([])
    #     data_dict[f"{i}_inc"] = np.array([])
    #     data_dict[f"{i}_a"] = np.array([])   
    #     data_dict[f"{i}_e"] = np.array([]) 

    # Running integration over time
    prog_frac = N/100
    time = np.linspace(t_start, t_end, N)   
    for iteration, t in enumerate(time):

        # Print status
        if int(iteration % prog_frac) == 0:
            print(f"  Progress = {int(iteration/prog_frac)}%", end='\r')

        # Integrate simulation
        sim.integrate(t)
        # data_dict["time"] = np.append(data_dict["time"], t)

    #     # Saving attributes per particle in simulation
    #     for i, p in enumerate(ps):

    #         # Get x, y and z position of particle
    #         pos = p.xyz
    #         x, y, z = pos[0], pos[1], pos[2]

    #         # Get vx, vy and vz position of particle
    #         speed = p.vxyz
    #         vx, vy, vz = speed[0], speed[1], speed[2]

    #         # Get inclination, eccentricity and semi major axis
    #         if i > 0:
                
    #             inc, e, a = p.inc, p.e, p.a
            
    #         else:

    #             inc, e, a = 0, 0, 0

    #         # Save all attributes into dictionary
    #         data_dict[f"{i}_x"] = np.append(data_dict[f"{i}_x"], x)
    #         data_dict[f"{i}_y"] = np.append(data_dict[f"{i}_y"], y)
    #         data_dict[f"{i}_z"] = np.append(data_dict[f"{i}_z"], z)
    #         data_dict[f"{i}_vx"] = np.append(data_dict[f"{i}_vx"], vx)
    #         data_dict[f"{i}_vy"] = np.append(data_dict[f"{i}_vy"], vy)
    #         data_dict[f"{i}_vz"] = np.append(data_dict[f"{i}_vz"], vz)
    #         data_dict[f"{i}_inc"] = np.append(data_dict[f"{i}_inc"], inc)
    #         data_dict[f"{i}_a"] = np.append(data_dict[f"{i}_a"], a)
    #         data_dict[f"{i}_e"] = np.append(data_dict[f"{i}_e"], e)
        
    print("")
    # # Write data to csv
    # df = pd.DataFrame(data_dict)
    # df.to_csv(OUTPUT_FOLDER + filename + ".csv", index=False)


def AnalyseConvergence():

    # Make integrations
    dts = [50, 10, 5, 2, 1, .5, .2, .1, .05, .01]
    # for dt in dts:
    #     sim = InitSolarSytem()
    #     # func.AddPlanet(sim)
    #     Integrate(sim, filename=f"Norm_dt{dt}", overwrite=True, dt=dt, years=1)

    # Make subplot
    plt.subplot(2,1,1)
    index = 1

    # Get archives
    for dt in dts:
        sa = rebound.SimulationArchive(f"simulations/Norm_dt{dt}.bin")

        # Initialize arrays
        t = []
        E = []
        a = []
        E_0 = sa[0].calculate_energy()
        a_0 = sa[0].particles[index].a

        # Loop over archive
        for sim in sa:
        
            if sim.t > 0:
                t.append(sim.t/(2*np.pi))
            else:
                t.append(sim.t)

            E.append(abs((sim.calculate_energy() - E_0)/E_0))   
            a.append(abs((sim.particles[index].a - a_0)/a_0))

        # plt.plot(t, E, alpha=.9, lw=.5, label=f"dt = {dt} days")
        plt.plot(t, a, alpha=.9, lw=.5, label=f"dt = {dt} days")
        
    # plt.title("Difference total energy of solar system")
    plt.plot([0, 1], [1e-5, 1e-5], linestyle='--', color="r", lw=0.5, label="Variation = 1e-5")
    plt.title(f"Difference a of {LABELS[index]}")
    plt.yscale("log")
    # plt.ylim(-1.7e-9, 1.7e-9)
    plt.legend(loc='upper right', prop={'size': 6})
    plt.xlabel("t [years]")
    # plt.ylabel("|dE|/E_0")
    plt.ylabel("|da|/a_0 [AU]")

    # Get archives
    plt.subplot(2,1,2)
    dE = []
    da = []
    for dt in dts:
        sa = rebound.SimulationArchive(f"simulations/Norm_dt{dt}.bin")

        # Initialize arrays
        a = []
        E = []
        a_0 = sa[0].particles[index].a
        E_0 = sa[0].calculate_energy()

        # Loop over archive
        for sim in sa:
        
            if sim.t > 0:
                t.append(sim.t/(2*np.pi))
            else:
                t.append(sim.t)

            E.append(abs((sim.calculate_energy() - E_0)/E_0))   
            a.append(abs((sim.particles[index].a - a_0)/a_0)) 
        
        dE.append(max(E))
        da.append(max(a))
        plt.plot(dt, max(a), marker="o", alpha=.9, label=f"dt = {dt} days")
        # plt.plot(dt, max(E), marker="o", alpha=.9, label=f"dt = {dt} days")
   
    # Make subplot
    # plt.title("Max difference total energy of solar system")
    plt.title(f"Max difference a of {LABELS[index]}")
    plt.xlim(max(dts), min(dts))
    # plt.plot(dts, dE, linestyle='--', marker='o', color='b', lw=0.5)
    plt.plot(dts, da, linestyle='--', color='k', lw=0.5)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("dt [days]")
    # plt.ylabel("|dE|/E_0 max")
    plt.ylabel("|da|/a_0 max [AU]")

    # Make sure subplots dont overlap
    plt.tight_layout()

    # Save figure
    filename = f"Convergence_A_{LABELS[index]}"
    plt.savefig(f"{PLOTS_FOLDER}{filename}.png", dpi=200)


def AnalyseSolarSystemPerturbation():
    """Analyse accuracy of rebound"""

    print("AnalyseSolarSystemPerturbation...")

    # Loop over archive and save planet orbit properties
    time_int = []
    incs_int = [[],[],[],[],[],[],[],[],[],[]]
    es_int = [[],[],[],[],[],[],[],[],[],[]]
    sa = rebound.SimulationArchive(f"{SIMULATIONS_FOLDER}Integration_1990_dt0.2_t50000.bin")
    sim_start = sa[0]
    ps_start = sim_start.particles

    # Loop over archive
    for sim in sa:

        # Get planets
        i = 1
        ps = sim.particles[i:]

        # Loop over planets
        for p in ps:

            # Save perturbation
            incs_int[i].append(abs(p.inc - ps_start[i].inc))
            es_int[i].append(abs(p.e - ps_start[i].e))

            i += 1 

        # Add time 
        time_int.append(sim.t/(2*np.pi))

    # Get max perturbation
    for i in range(1, 9):
        
        print("-----------------------------")
        print(LABELS[i])
        print(f"max_inc = {max(incs_int[i])}")
        print(f"max_e = {max(es_int[i])}")

    # # Subplot
    # plt.subplot(2,1,1)
    # plt.title("Solar system perturbations")
    # for i in range(1,9):  
    #     plt.plot(time_int, incs_int[i], color=f"C{i}", label=f"{LABELS[i]}", alpha=.7, lw=.7)      
    # plt.plot([0, 50000], [2.22e-7, 2.22e-7], color="k", linestyle="--", label="Rebound error")
    # plt.xlabel("time [years]")
    # plt.ylabel("$|inc - inc_0|$")
    # plt.yscale("log")
    # plt.legend(loc='upper right', prop={'size': 6})

    # # Subplot
    # plt.subplot(2,1,2)
    # for i in range(1,9): 
    #     plt.plot(time_int, es_int[i], color=f"C{i}", label=f"{LABELS[i]}", alpha=.7, lw=.7)    
    # plt.plot([0, 50000], [5.55e-5, 5.55e-5], color="k", linestyle="--", label="Rebound error") 
    # plt.xlabel("time [years]")
    # plt.ylabel("$|e - e_0|$")
    # plt.yscale("log")
    # plt.legend(loc='upper right', prop={'size': 6})

    # plt.tight_layout()
    # plt.savefig(PLOTS_FOLDER + "AnalysePerturabtions_t50000_1.png", dpi=200)


def AnalyseAccuracyErrors():
    """Calculate rebound integration error"""

    print("AnalyseAccuracyErrors...")

    # Initialize
    years = range(1990, 2020, 1)
    time_real = []
    dincs = [[],[],[],[],[],[],[],[],[]]
    des = [[],[],[],[],[],[],[],[],[]]

    # Define archive
    sa = rebound.SimulationArchive(f"{SIMULATIONS_FOLDER}Integration_1990_dt0.2_t30.bin")

    # Open simulations and save orbital parameters
    for year in years:      

        # Get planets
        i = 1
        sim_real = rebound.Simulation(SIMULATIONS_FOLDER + f"ss_{year}.bin")
        sim_int = sa.getSimulation(2 * np.pi * (year-years[0]))
        
        # Get particles
        ps_real = sim_real.particles[i:]
        ps_int = sim_int.particles[i:]
        
        # Save properties of planets
        for p_real, p_int in zip(ps_real, ps_int):

            dincs[i].append(abs(p_real.inc - p_int.inc))
            des[i].append(abs(p_real.e - p_int.e))

            i += 1

        # Add time 
        time_real.append(year - years[0])

    # Get max values
    errbar_inc = 0
    for dinc in dincs[1:]:
        maxim = max(dinc)

        if maxim > errbar_inc:
            errbar_inc = maxim

    # Get max values
    errbar_e = 0
    for de in des[1:]:
        maxim = max(de)

        if maxim > errbar_e:
            errbar_e = maxim

    print(f"  max error inc = {errbar_inc}")
    print(f"  max error e = {errbar_e}")

    # Subplot
    plt.subplot(2,1,1)
    plt.title("Rebound integration accuracy")
    for i in range(1,9):
        # plt.plot(time_real, incs_real[i], color=f"C{i}", linestyle="--", lw=.6)   
        plt.plot(time_real, dincs[i], color=f"C{i}", label=f"{LABELS[i]}", alpha=.5)
        plt.plot([0, 29], [errbar_inc, errbar_inc], linestyle="--", color="k")         
    plt.xlabel("time [years]")
    plt.ylabel("|dinc|")
    plt.yscale("log")
    plt.legend(loc='upper right', prop={'size': 6})

    # Subplot
    plt.subplot(2,1,2)
    for i in range(1,9):
        # plt.plot(time_real, es_real[i], color=f"C{i}", linestyle="--", lw=.6)   
        plt.plot(time_real, des[i], color=f"C{i}", label=f"{LABELS[i]}", alpha=.5) 
        plt.plot([0, 29], [errbar_e, errbar_e], linestyle="--", color="k")    
    plt.xlabel("time [years]")
    plt.ylabel("|de|")
    plt.yscale("log")
    plt.legend(loc='upper right', prop={'size': 6})

    plt.tight_layout()
    plt.savefig(PLOTS_FOLDER + "AccuracyError_1.png", dpi=200)    


def AnalyseArchive():
    """Calculate rebound integration error"""

    print("AnalyseArchive...")

    # Open archive
    sa = rebound.SimulationArchive(f"{SIMULATIONS_FOLDER}Integration_1990_dt0.2_t50000.bin")

    # Intialize
    a = [[],[],[],[],[],[],[],[],[],[]]
    e = [[],[],[],[],[],[],[],[],[],[]]
    inc = [[],[],[],[],[],[],[],[],[],[]]

    # Get info
    for sim in sa:
        
        # Get planets
        i = 1
        ps = sim.particles[i:]

        # Save info
        for p in ps:

            a[i].append(p.a)
            e[i].append(p.e)
            inc[i].append(p.inc)

            i += 1

    orbit_param = "e"

    if orbit_param == "e":
        for i in range(1, 9):
            plt.plot(a[i][0], e[i][0], marker=">", markersize=3, color=f"C{i}", label=LABELS[i])
            plt.plot(a[i][-1], e[i][-1], marker="<", markersize=3, color=f"C{i}")
            plt.plot(a[i], e[i], marker=".", markersize=0.001, color=f"C{i}", alpha=.4)

        plt.legend()
        plt.xscale("log")
        plt.xlabel("a [AU]")
        plt.ylabel("e")
        plt.savefig(f"{PLOTS_FOLDER}Analyse_e.png", dpi=200)

    if orbit_param == "inc":
        for i in range(1, 9):
            plt.plot(a[i][0], inc[i][0], marker=">", markersize=3, color=f"C{i}", label=LABELS[i])
            plt.plot(a[i][-1], inc[i][-1], marker="<", markersize=3, color=f"C{i}")
            plt.plot(a[i], inc[i], marker=".", markersize=.001, color=f"C{i}", alpha=.4)

        plt.legend()
        plt.xscale("log")
        plt.xlabel("a [AU]")
        plt.ylabel("inc")
        # plt.savefig(f"{PLOTS_FOLDER}Analyse_inc.png", dpi=200)

def Unload(filename):
    """Unload csv into pandas dataframe"""
  
    print("Unload...")

    # Check if file exists
    if not os.path.exists(OUTPUT_FOLDER + filename + ".csv"):

        print(f"WARNING: The file {filename}.csv cannot be found")
        return 

    df = pd.read_csv(f"output/{filename}.csv")

    return df    


def AnimatePositions(df, filename="auto", dpi=100, fps=30):
    """Make animation of positions"""

    print("AnimatePositions...")

    # Check if dataframe is not empty
    if df is None:

        print("WARNING: The dataframe you passed is empty")
        return

    # If saving directory does not exists make one
    if not os.path.isdir(ANIMATIONS_FOLDER):

        os.mkdir(ANIMATIONS_FOLDER)
        print(f"NOTE: created this folder {ANIMATIONS_FOLDER}")

    # Check if you have the desired animation write
    if "ffmpeg" not in animation.writers.list():

        print("WARNING: You do not have the ffmpeg writer, ask Simon :p")
        return

    # Setup plot variables
    fig = plt.figure(figsize=(10,10), dpi=dpi)
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.xlabel('x [AU]')
    plt.ylabel('y [AU]')
    plt.title('Animation of solar system with planet R9')

    # Initialize plots for all particles pt = point and ln = line
    pt0, = plt.plot([], [], 'k*', markersize=4, animated=True)
    ln0, = plt.plot([], [], 'k-', markersize=1, alpha=0.5, animated=True)
    pt1, = plt.plot([], [], 'ko', markersize=2, animated=True)
    ln1, = plt.plot([], [], 'k-', markersize=1, alpha=0.5, animated=True)
    pt2, = plt.plot([], [], 'ko', markersize=2, animated=True)
    ln2, = plt.plot([], [], 'k-', markersize=1, alpha=0.5, animated=True)
    pt3, = plt.plot([], [], 'ko', markersize=2, animated=True)
    ln3, = plt.plot([], [], 'k-', markersize=1, alpha=0.5, animated=True)
    pt4, = plt.plot([], [], 'ko', markersize=2, animated=True)
    ln4, = plt.plot([], [], 'k-', markersize=1, alpha=0.5, animated=True)
    pt5, = plt.plot([], [], 'ko', markersize=2, animated=True)
    ln5, = plt.plot([], [], 'k-', markersize=1, alpha=0.5, animated=True)
    pt6, = plt.plot([], [], 'ko', markersize=2, animated=True)
    ln6, = plt.plot([], [], 'k-', markersize=1, alpha=0.5, animated=True)
    pt7, = plt.plot([], [], 'ko', markersize=2, animated=True)
    ln7, = plt.plot([], [], 'k-', markersize=1, alpha=0.5, animated=True)
    pt8, = plt.plot([], [], 'ko', markersize=2, animated=True)
    ln8, = plt.plot([], [], 'k-', markersize=1, alpha=0.5, animated=True)
    pt9, = plt.plot([], [], 'ro', markersize=4, animated=True)
    ln9, = plt.plot([], [], 'r-', markersize=1, alpha=0.5, animated=True)

    def animate(i, N):

        print(f"  Generating animation frame {i}/{N}", end='\r')

        # Select data range
        pts = df.loc[[i]] 
        lns = df.loc[0 if i-100 < 0 else i-100:i]

        # Set data
        pt0.set_data(pts["0_x"].to_numpy(), pts["0_y"].to_numpy())
        ln0.set_data(lns["0_x"].to_numpy(), lns["0_y"].to_numpy())
        pt1.set_data(pts["1_x"].to_numpy(), pts["1_y"].to_numpy())
        ln1.set_data(lns["1_x"].to_numpy(), lns["1_y"].to_numpy())
        pt2.set_data(pts["2_x"].to_numpy(), pts["2_y"].to_numpy())
        ln2.set_data(lns["2_x"].to_numpy(), lns["2_y"].to_numpy())
        pt3.set_data(pts["3_x"].to_numpy(), pts["3_y"].to_numpy())
        ln3.set_data(lns["3_x"].to_numpy(), lns["3_y"].to_numpy())
        pt4.set_data(pts["4_x"].to_numpy(), pts["4_y"].to_numpy())
        ln4.set_data(lns["4_x"].to_numpy(), lns["4_y"].to_numpy())
        pt5.set_data(pts["5_x"].to_numpy(), pts["5_y"].to_numpy())
        ln5.set_data(lns["5_x"].to_numpy(), lns["5_y"].to_numpy())
        pt6.set_data(pts["6_x"].to_numpy(), pts["6_y"].to_numpy())
        ln6.set_data(lns["6_x"].to_numpy(), lns["6_y"].to_numpy())
        pt7.set_data(pts["7_x"].to_numpy(), pts["7_y"].to_numpy())
        ln7.set_data(lns["7_x"].to_numpy(), lns["7_y"].to_numpy())
        pt8.set_data(pts["8_x"].to_numpy(), pts["8_y"].to_numpy())
        ln8.set_data(lns["8_x"].to_numpy(), lns["8_y"].to_numpy())
        pt9.set_data(pts["9_x"].to_numpy(), pts["9_y"].to_numpy())
        ln9.set_data(lns["9_x"].to_numpy(), lns["9_y"].to_numpy())
        
    # Setup saving animation
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=fps, metadata=dict(artist='Simon van Eeden'), bitrate=1800)

    # Animation variables
    N = len(df)

    # Make antimation
    ani = matplotlib.animation.FuncAnimation(fig, 
        func=animate,
        fargs=(N,),
        frames=N, 
        repeat=True)

    # Create automatic filename
    if filename == "auto":
        
        # Get string if current date and time
        current_time = datetime.now().strftime("%Y-%m-%d_%H-%M")

        filename = "Positions_" + current_time

    ani.save(f"{ANIMATIONS_FOLDER}{filename}.mp4", writer=writer, dpi=dpi)


def SaveSim(sim, filename="auto", overwrite=False):
    """Save simulation to /simulations folder"""

    print("SaveSim...")

    # Check if saving directory exists otherwise make one
    if not os.path.isdir(SIMULATIONS_FOLDER):

        print(f"NOTE: created this folder {SIMULATIONS_FOLDER}")
        os.mkdir(SIMULATIONS_FOLDER)

    # Check if filename is already used
    if os.path.exists(SIMULATIONS_FOLDER + filename + ".bin") and not overwrite:

        print(f"WARNING: The file {filename}.bin already exists, auto filename has been used")
        filename = "auto"

    # Create automatic filename
    if filename == "auto":
        
        # Get string if current date and time
        current_time = datetime.now().strftime("%Y-%m-%d_%H-%M")

        filename = "Sim_" + current_time

    # Save simulation object to file
    sim.save(SIMULATIONS_FOLDER + filename + ".bin")
