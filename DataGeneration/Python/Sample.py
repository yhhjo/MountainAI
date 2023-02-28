import matplotlib.pyplot as plt
import torch 
import numpy as np
import random
import os.path 
from GeoHistory import generate_exhumation_history, generate_lithology_and_heat

# CONSTANTS
EPSILON = .01

class Sample():
    def __init__(self, params):
        self.duration = params.duration
        self.n_horizons = params.n_horizons
        self.surf_temp = params.surf_temp
        self.horizon_dep = params.horizon_dep 

        # Exhumation history hyperparameters
        self.num_events = params.num_events
        self.uplift_rate_mean = params.uplift_rate_mean
        self.uplift_rate_std = params.uplift_rate_std
        
        # Lithology and heat hyperparameters
        self.b_cond_mean = params.b_cond_mean
        self.surf_hp_mean = params.surf_hp_mean
        self.hp_dep_mean = params.hp_dep_mean
        self.surf_hf_mean = params.surf_hf_mean
        self.std_scale = params.std_scale 
        
        # Generate exhumation history
        p = [self.duration, self.num_events, self.uplift_rate_mean, self.uplift_rate_std, EPSILON]
        self.exhumation_history, self.exhumation_rates = generate_exhumation_history(*p)

        # Generate lithology and thermal history
        p = [self.b_cond_mean, self.surf_hp_mean, self.hp_dep_mean, self.surf_hf_mean, self.std_scale]
        r = generate_lithology_and_heat(*p)
        self.basement_cond, self.surface_hp, self.hp_depth, self.surf_hf = r
        
    def summary(self):
        print(f"duration = {self.duration}")
        print(f"n_horizons = {self.n_horizons}")
        print(f"horizon_dep = {self.horizon_dep}")
        print(f"num_events = {self.num_events}")
        print(f"uplift_rate_mean = {self.uplift_rate_mean}")
        print(f"uplift_rate_std = {self.uplift_rate_std}")
        print(f"exhumation_history = {self.exhumation_history}")
        print(f"exhumation_rates = {self.exhumation_rates}")
        print(f"basement_cond = {self.basement_cond}")
        print(f"surface_hp = {self.surface_hp}")
        print(f"hp_depth = {self.hp_depth}")
        print(f"surf_hf = {self.surf_hf}")
        print(f"surf_temp = {self.surf_temp}")

    # Creates output file readable by tqtec
    def toOutputFile(self, filename):
        with open(filename, 'w') as f:
            f.write(f"{filename}\n")
            f.write(f"{self.duration:10.0f}{2:10d}{self.surf_temp:10.4f}{self.surf_hf:10.4f}{self.basement_cond:10.4f}{self.surface_hp:10.4f}{self.hp_depth:10.4f}\n")
            f.write(f"{0:10}\n")
            f.write(" ".join(f"{e:10.4f}" for e in self.horizon_dep) + "\n")
            f.write(f"{0:10}\n")
            f.write(f"{self.num_events:10d}\n")
            for r in self.exhumation_history:
                for e in r:
                    f.write(f"{e:10.4f}")
                f.write(f"\n")
            f.write(f"{0:10}\n{0:10}")
            