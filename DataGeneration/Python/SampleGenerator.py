import matplotlib.pyplot as plt
import torch 
import numpy as np
import os 
import subprocess

from Parser import parse_tensor
import Sample
import Orogens
import concurrent.futures

# Constants
FNAME = "sample.in" 
OUTNAME = "sample.out" 
CLOSURE_NAME = "TC"
TEMPS = ["60", "120", "180", "240"]
PATH_TO_TQTEC = "../Fortran/"

def cleanup(i, names):
    for f in names:
        if os.path.exists(f):
            os.remove(f)

# Multithreaded function that executes tqtec and rdtqtec, storing the temporary input and output files in 'tmp' before 
# deleting them. 
def process_sample_thread(i, orogen, tmp_dir):
    i = str(i)
    fname=tmp_dir+'/'+FNAME+i
    outname=tmp_dir+'/'+OUTNAME+i
    closure_name=tmp_dir+'/'+CLOSURE_NAME+i
    names = [fname, outname, closure_name, TEMPS]

    tqtec = [PATH_TO_TQTEC+"tqtec", "-f", fname, "-o", outname]
    rdtqtec = [PATH_TO_TQTEC+"rdtqtec", outname, "-closure", closure_name, *TEMPS]

    if (orogen=="CC"):
        p = Orogens.Continental_Collision()

    s = Sample.Sample(p)
    try:
        s.toOutputFile(fname)
        subprocess.run(tqtec, check=True)
        subprocess.run(rdtqtec, check=True)
        tc = parse_tensor(closure_name, plot=False)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e} idx: {i}")
        cleanup(i, names[:-1])
        return 0
    cleanup(i, names[:-1])
    return tc, s

# Gathers temperature-time-elevation data from rdtqtec and stores the results in a 2d array of 'input data' and 'closure temperature' results. Every 'chunk_size' iterations, it saves its progress to 'output dir' as a numpy array. 
def generate_models(n, chunk_size, orogen, output_dir, tmp_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
        
    results = []
    for chunk in range(chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(process_sample_thread, i, orogen, tmp_dir) for i in range(int(n/chunk_size))]
            for idx, future in enumerate(concurrent.futures.as_completed(futures)):
                result = future.result()
                if result:
                    results.append(result)
        np.save(f"{output_dir}/results_CHUNK_{chunk}.npy", np.array(results, dtype=object))
    return results

# Returns chunks as a list of numpy arrays
def read_chunk(file):
    return np.load(file, allow_pickle=True)