import matplotlib.pyplot as plt
import torch 
import numpy as np
import random

# CONSTANTS FOR CONTINENT COLLISION
CC_DUR_MEAN = 25
CC_DUR_STD = 5
N_HORIZONS = 10 
HORIZON_DEP = [5, 8, 11, 13, 15, 17, 19, 21, 23, 25]
MIN_EVENTS = 2
MAX_EVENTS = 5
CC_UPLIFT_RATE_MEAN = 1.5
CC_UPLIFT_RATE_STD = 1
CC_B_COND_MEAN = 2.5
CC_HP_MEAN = 2
CC_HP_DEP_MEAN = 45
CC_SURF_HF_MEAN = 100
STD_SCALE = .1
SURF_TEMP_MEAN = 5
SURF_TEMP_STD = 5

# Sample model parameters with means specified in get_params
class Continental_Collision():
    def __init__(self):
        self.duration = round(random.normalvariate(CC_DUR_MEAN, CC_DUR_STD), 0)
        self.n_horizons = N_HORIZONS
        self.horizon_dep = HORIZON_DEP
        self.num_events = round(random.uniform(MIN_EVENTS, MAX_EVENTS))
        self.uplift_rate_mean = CC_UPLIFT_RATE_MEAN
        self.uplift_rate_std = CC_UPLIFT_RATE_STD
        self.b_cond_mean = CC_B_COND_MEAN
        self.surf_hp_mean = CC_HP_MEAN
        self.hp_dep_mean = CC_HP_DEP_MEAN
        self.surf_hf_mean = CC_SURF_HF_MEAN
        self.std_scale = STD_SCALE
        self.surf_temp = round(random.normalvariate(SURF_TEMP_MEAN, SURF_TEMP_STD), 0)