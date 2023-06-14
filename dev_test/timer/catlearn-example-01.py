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


# Timestamp
def timestamp():
    return datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")


# Local timer
class Timer:
    def __init__(self):
        self.timings = {}

    @staticmethod
    def timestamp():
        return (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), time.time())

    def start(self, name):
        assert name not in self.timings.keys()
        dt, t = self.timestamp()
        self.timings[name] = {
            "title": name,
            "start_time": t,
            "start_time_str": dt,
            "end_time": None,
            "end_time_str": None,
            "duration": None,
            "duration_str": None,
        }

        print(f"[ {dt} ] START  | {name} ...")

    def stop(self, name):
        dt, t = self.timestamp()

        t_obj = self.timings[name]
        t_obj["end_time"] = t
        t_obj["end_time_str"] = dt
        t_obj["duration"] = t_obj["end_time"] - t_obj["start_time"]
        t_obj["duration_str"] = f"{t_obj['duration']:.6f}"

        print(f"[ {dt} ] STOP   | {name} ... {t_obj['duration']} seconds")

    def report(self):
        headers = ["Title", "Start Time", "End Time", "Duration"]
        table = []
        for k, v in self.timings.items():
            row = [
                v["title"],
                v["start_time_str"],
                v["end_time_str"],
                v["duration"],
            ]
            table.append(row)
        kwargs = dict(
            floatfmt=".6f",
            tablefmt="grid",
            numalign="decimal",
            # stralign="left",
        )
        print(tabulate(table, headers, **kwargs))


# Define number of images:
n_images = 7

print(f"{timestamp()} Optimizing NEB images path (CatLearn METHOD):")

T = Timer()

T.start("Catlearn NEB initialization")
neb_catlearn = MLNEB(
    start="initial.traj",
    end="final.traj",
    ase_calc=ase_calculator(),
    n_images=n_images,
    interpolation="linear",
    restart=False,
)
T.stop("Catlearn NEB initialization")

T.start("Catlearn NEB optimization")
neb_catlearn.run(
    fmax=0.01,
    trajectory="ML-NEB.traj",
    # full_output=True,
)
T.stop("Catlearn NEB optimization")

T.report()
