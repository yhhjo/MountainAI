import numpy as np

def process_input(dataset, normalize=True):
    """Pulls inputs from dataset and returns an array (dataset_size, 10, 5) normalized across the 1st dimension"""
    y = np.array([*dataset])[:,0]
    y = np.array([np.array(z) for z in y])

    if normalize:
        std = np.std(y, axis=1, keepdims=True)
        mean = np.mean(y, axis=1, keepdims=True)
        std[std==0] = .0001
        y = (y - mean) / std
    y = np.transpose(y, (0, 2, 1))
    return y
    
def process_output(dataset, d_output: int) -> np.array([]):
    """Input: dataset. Returns the concatenated values of Sample object 0-padded to be size d_output"""
    Y = np.array([*dataset])[:,1]
    ret = np.zeros((Y.shape[0], d_output))
    for i, y in enumerate(Y):
        t = [y.duration, y.num_events, y.basement_cond, y.surface_hp, y.hp_depth, y.surf_hf, y.surf_temp, 
        *y.exhumation_history.flatten(), *y.exhumation_rates]
        t = np.pad(t, (0, d_output-len(t)), constant_values=0) #Pad
        ret[i] = t
    return ret

def formattedLoader(d_output, dataset_path):
    """Returns data formatted specifically for transformer"""
    dataset = np.load(dataset_path,allow_pickle=True)
    X = process_input(dataset)
    Y = process_output(dataset, d_output)
    formatted_data = []
    for (x, y) in zip(X, Y):
        formatted_data.append((x, y))

    return formatted_data

def rawLoader(dataset_path: str)->np.array([]):
    return np.load(dataset_path,allow_pickle=True)