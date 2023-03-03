import numpy as np
import Sample

def process_input(x)->np.array([]):
    """Takes input of shape H_thermal and concatenates the normalized vectors (age, AHe, Aft, ZHe, Zft)
    Output is shape d_input. 
    """
    x = np.array([*x])
    std = np.std(x, axis=1, keepdims=True)
    std[std==0] = .0001
    x = ((x - np.mean(x, axis=1, keepdims=True)) / std).flatten() # Normalize
    return x

def process_output(y: Sample, d_output: int) -> np.array([]):
    """y: Sample object. Returns a concatenated, padded array of size d_output"""
    y = [y.duration, y.num_events, y.basement_cond, y.surface_hp, y.hp_depth, y.surf_hf, y.surf_temp, 
    *y.exhumation_history.flatten(), *y.exhumation_rates]
    y = np.pad(y, (0, d_output-len(y)), constant_values=0) #Pad
    return y

def formattedLoader(d_input, d_output, dataset_path):
    """Returns data formatted specifically for transformer"""
    dataset = np.load(dataset_path,allow_pickle=True)
    formatted_data = []
    for x, y in dataset:
        x = process_input(x)
        y = process_output(y, d_output)
        formatted_data.append((x,y))
    return formatted_data

def rawLoader(dataset_path: str)->np.array([]):
    return np.load(dataset_path,allow_pickle=True)