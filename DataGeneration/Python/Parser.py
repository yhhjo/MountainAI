import argparse
import matplotlib.pyplot as plt
import torch 
import numpy as np

def parse_tensor(fname, plot):
    ahe = [] 
    aft = [] 
    zhe = [] 
    zft = [] 
    ages = []
    chronometer = ''
    
    with open(fname, 'r') as f:
        for line in f:
            if '>' in line:
                if '60' in line:
                    chronometer = '60'
                elif '120' in line:
                    chronometer = '120'
                elif '180' in line:
                    chronometer = '180'
                elif '240' in line:
                    chronometer = '240'
                continue
            
            cols = line.split()
            if len(cols) != 2:
                continue
            
            ages.append(eval(cols[1]))
            elevation = eval(cols[0])

            if chronometer == '60':
                ahe.append(elevation)
            elif chronometer == '120':
                aft.append(elevation)
            elif chronometer == '180':
                zhe.append(elevation)
            elif chronometer == '240':
                zft.append(elevation)
    ages = ages[:max(len(ahe), len(aft), len(zhe), len(zft))]
    ages = -(ages - np.max(ages))
    if plot:
        plt.plot(ages, aft)
        plt.plot(ages, ahe)
        plt.plot(ages, zhe)
        plt.plot(ages, zft)
    return (ages, ahe, aft, zhe, zft)