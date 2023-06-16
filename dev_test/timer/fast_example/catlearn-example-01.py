from ase.build import fcc100, add_adsorbate
from ase.calculators.emt import EMT
from ase.io import read
from ase.constraints import FixAtoms
from ase.neb import NEB
from ase.optimize import BFGS, MDMin
import matplotlib.pyplot as plt
from catlearn.optimize.mlneb import MLNEB
from ase.neb import NEBTools
from catlearn.optimize.tools import plotneb
from datetime import datetime
import time
import sys
from tabulate import tabulate
import numpy as np
import os


""" 
    Toy model for the diffusion of a Au atom on an Al(111) surface.  
    This example contains: 
    1. Optimization of the initial and final end-points of the reaction path. 
    2.A. NEB optimization using CI-NEB as implemented in ASE. 
    2.B. NEB optimization using our machine-learning surrogate model.
    3. Comparison between the ASE NEB and our ML-NEB algorithm.
"""


# Calculator
def ase_calculator():
    return EMT()


# Define number of images:
n_images = 7

repeats = 1

info = np.zeros((repeats, 3))

for i in range(repeats):
    neb_catlearn = MLNEB(
        start="initial.traj",
        end="final.traj",
        ase_calc=ase_calculator(),
        n_images=n_images,
        interpolation="linear",
        restart=False,
    )

    neb_catlearn.run(
        fmax=1,
        trajectory="ML-NEB.traj",
        # full_output=True,
    )

    neb_catlearn.timer.report()
    neb_catlearn.timer.summarize(
        [
            "MLNEB Run",
            "Train GP Model",
            "ML NEB optimization",
            "Get results",
            "Acquisition function",
            "Add & evaluate",
            "Store results",
            "Check convergence",
        ]
    )

    ct_obj = neb_catlearn.timer.data_summary
    info[i, 0] = ct_obj["Train GP Model"]["sum"]
    info[i, 1] = ct_obj["MLNEB Run"]["sum"]
    info[i, 2] = ct_obj["Add & evaluate"]["sum"]

thread = os.environ["OMP_NUM_THREADS"]
np.savetxt(f"info-{thread}.txt", info, delimiter=",")

# from pprint import pprint

# pprint(neb_catlearn.timer.data_summary)

# sys.exit()
