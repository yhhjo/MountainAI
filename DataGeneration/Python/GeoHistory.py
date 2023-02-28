import matplotlib.pyplot as plt
import torch 
import numpy as np
import random



def generate_exhumation_history(model_duration, num_events, UPLIFT_RATE_MEAN, UPLIFT_RATE_STD, EPSILON):
    MAX_OFFSET = .4 #Longest model can wait for first event 
    MAX_DURATION = .7 #Each event can use, at most, 70% of remaining duration
    DURATION_STD_SCALE = .25
    
    events = []
    exhumation_rates = []
        
    # Generate the first event start time 
    offset = round(random.uniform(0, model_duration*MAX_OFFSET), 1)
    remaining_time = model_duration - offset

    mean = remaining_time/num_events
    
    for i in range(0, num_events):
        event_duration = abs(round(random.normalvariate(mean, mean*DURATION_STD_SCALE), 1))

        #  Update params
        event_start_time = model_duration - remaining_time
        remaining_time -= event_duration

        uplift = round(abs(random.normalvariate(UPLIFT_RATE_MEAN, UPLIFT_RATE_STD))*event_duration, 1)
        events.append([event_start_time, event_duration, uplift])
        exhumation_rates.append(round(uplift/(event_duration + EPSILON), 2))
    return np.array(events), exhumation_rates


# Arguments represent means. Std_scale will draw samples from a normal distribution where std = mean*std_scale.
def generate_lithology_and_heat(b_cond, surf_hp, hp_dep, surf_hf, std_scale):
    basement_cond = round(abs(random.normalvariate(b_cond, b_cond*std_scale)), 2)
    surface_heat_prod = round(abs(random.normalvariate(surf_hp, surf_hp*std_scale)), 2)
    heat_prod_depth = round(abs(random.normalvariate(hp_dep, hp_dep*std_scale)), 0)
    surf_hf = round(abs(random.normalvariate(surf_hf, surf_hf*std_scale)), 1)
    
    # surf_hf must be greater than heat generated in crust
    if surf_hf <= heat_prod_depth*surface_heat_prod: 
        surf_hf = round(heat_prod_depth*surface_heat_prod + abs(random.normalvariate(3,1)), 3)
    return basement_cond, surface_heat_prod, heat_prod_depth, surf_hf
