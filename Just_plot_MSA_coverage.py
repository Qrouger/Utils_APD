#This python script is made for plot MSA with coverage, LDDT and PAE from Alphapulldown 1.0.4 outputs. It's a modification of a Alphafold #plot script create by Sébastien Lemal. https://blog.biostrand.ai/explained-how-to-plot-the-prediction-quality-metrics-with-alphafold2
import os
import glob
import pickle
import json
import numpy as np
import matplotlib.pyplot as plt
 
class ARG:
    def __init__(self, repo):
        self.input_dir = repo
        self.output_dir = repo
        self.name = pickle_file
        self.pickle_file = pickle_file

print("Enter PATH to file")
repo = [str(input())]

print ("Enter name of .pkl") 
pickle_file = input()

for r in repo:
    args = ARG(r)
    pre_feature_dict = pickle.load(open(f'{args.input_dir}/{args.pickle_file}.pkl','rb'))
    feature_dict = pre_feature_dict.feature_dict
    
def generate_output_images(feature_dict, out_dir, name):
    msa = feature_dict['msa']
    seqid = (np.array(msa[0] == msa).mean(-1))
    seqid_sort = seqid.argsort()
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan
    final = non_gaps[seqid_sort] * seqid[seqid_sort, None]

    ###################### PLOT MSA WITH COVERAGE ####################
    
    plt.figure(figsize=(14, 4), dpi=100)
    plt.subplot(1, 2, 1)
    plt.title(f"Sequence coverage ({name})")
    plt.imshow(final,
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    plt.plot((msa != 21).sum(0), color='black')
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label="Sequence identity to query", )
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}coverage.pdf")


if __name__ == "__main__" :
    generate_output_images(feature_dict,args.output_dir if args.output_dir else args.input_dir, args.name)
