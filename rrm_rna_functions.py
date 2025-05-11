import os
import json
import glob
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from itertools import product
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from sklearn.neighbors import KernelDensity
from collections import Counter

import pathlib
abs_path = str(pathlib.Path(__file__).parent.resolve())

class RNAScoring():
    def __init__(self):
        self.res_list = ['R', 'K', 'H', 'D', 'E', 'Q', 'N', 'S', 'T', 'Y', 'F',
                    'W', 'M', 'L', 'I', 'V', 'A', 'G', 'P', 'C'
                    ]
        self.nt_list = ["A", "G", "C", "U"]

        # Protein bias
        fasta_rna_aln = abs_path + \
                "/alignment_data/RNAs_aligned_cluster_0_nt_full_refined.fasta"
        rna_seqs = list(SeqIO.parse(open(fasta_rna_aln), 'fasta'))
        self.rna_seqs_dict_full = {i.id: str(i.seq) for i in rna_seqs}
        self.bias_dict = Counter(
            [i.id.split("_")[0] + "_" + i.id.split("_")[1] for i in rna_seqs])

        self.fasta_prot_aln = abs_path + \
          "/alignment_data/rrm_bound_domains_aligned_processed_squeezed.fasta"
          
        prot_seqs = list(SeqIO.parse(open(self.fasta_prot_aln), 'fasta'))
        self.prot_seqs_dict_full = {"_".join(i.id.split("_")[:2]): str(i.seq) for i in prot_seqs}

        # Input data
        master_aln = abs_path + \
            '/alignment_data/RRM_master_interaction_info_with_nt_unbiased_cluster_0.json'
        self.aln_data = json.load(open(master_aln))
        rna_aln = abs_path + "/alignment_data/neighbours_rna_cluster_0.json"
        self.rna_data = json.load(open(rna_aln))
        self.entries_list = list(self.rna_seqs_dict_full.keys())

        self.hmm_to_master_pos = {
            29: 77, 31: 79, 33: 81, 34: 82, 37: 90, 58: 122, 60: 124, 62: 126,
            66: 135, 68: 137, 70: 139, 72: 141, 94: 166, 97: 169, 99: 171,
            101: 173, 102: 174, 103: 175, 104: 176, 105: 177
        }


    def load_scores(self, scores_folder=abs_path + "/alignment_data/precalculated_scores/*.pkl"):
        self.all_df_dict = {}
        for file in glob.glob(scores_folder):
            df = pd.read_pickle(file)
            res_num = int(file.split("/")[-1].split(".")[0].split("_")[0])
            nuc_num = int(file.split("/")[-1].split(".")[0].split("_")[1])
            self.all_df_dict[(res_num, nuc_num)] = df

        return self.all_df_dict

    def score_out_seq(self, rrm_seq, rna_seq, rna_pos_range):
        window_size = rna_pos_range[1] - rna_pos_range[0]

        # Now we'll store all windows and scores in order
        all_windows = []
        all_scores = []
        
        for i in range(len(rna_seq) - (window_size - 1)):
            rna_window = rna_seq[i: i + window_size]
            scores_by_pos = []

            for selected_pos, df in self.load_scores().items():
                if rna_pos_range[0] <= selected_pos[1] < rna_pos_range[1]:
                    nuc = rna_window[selected_pos[1] - rna_pos_range[0]]
                    res = rrm_seq[selected_pos[0]].upper()
                    if res != "-" and nuc != "-":
                        scores_by_pos.append(df.loc[res][nuc])

            score = np.nanmean(scores_by_pos)
            all_windows.append(rna_window)
            all_scores.append(score)
        
        # Create a dictionary for compatibility with existing code
        self.scores_dict = dict(zip(all_windows, all_scores))
        
        # But also return a list of (window, score) tuples to preserve order
        return self.scores_dict, list(zip(all_windows, all_scores))

    def plot_rna_kde(self, rna_seq, scores_dict, window_size):
        pos_file = open(
            "alignment_data/leave-one-out_positive_scores_avg_perc_20.txt",
            "r")
        pos_scores = [float(i.strip()) for i in pos_file.readlines()]
        neg_file = open(
            "alignment_data/leave-one-out_negative_scores_avg_perc_20.txt",
            "r")
        neg_scores = [float(i.strip()) for i in neg_file.readlines()]
        pos_kde = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(
            np.array(pos_scores).reshape(-1, 1))
        neg_kde = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(
            np.array(neg_scores).reshape(-1, 1))

        # Match the scores with the RNA fragments (Single RNA fragment)
        labels = []
        score_list = []
        kde_score = []

        for i in range(len(rna_seq) - (
                window_size - 1)):  # All the windows starting positions
            score = scores_dict[rna_seq[
                                i: i + window_size]] + 0.89  # Value from ROC-AUC
            score_list.append(score)
            labels.append(rna_seq[i: i + window_size])
            kde_score.append(float(
                np.exp(pos_kde.score_samples(np.array(score).reshape(-1, 1))) / \
                np.exp(neg_kde.score_samples(np.array(score).reshape(-1, 1)))))

        min_kde = float(
            np.exp(
                pos_kde.score_samples(np.array(-1.7 + 0.89).reshape(-1, 1))) / \
            np.exp(
                neg_kde.score_samples(np.array(-1.7 + 0.89).reshape(-1, 1))))
        max_kde = float(
            np.exp(pos_kde.score_samples(np.array(0 + 0.89).reshape(-1, 1))) / \
            np.exp(neg_kde.score_samples(np.array(0 + 0.89).reshape(-1, 1))))

        # Normalize the data
        kde_normalized = [
            (x - min_kde) / (max_kde - min_kde) for x in kde_score]

        my_cmap = plt.cm.get_cmap('RdYlGn')
        colors = my_cmap(kde_normalized)

        fig, ax = plt.subplots()
        y_pos = np.arange(len(score_list))
        ax.barh(y=y_pos, width=score_list, align='center', color=colors)

        sm = ScalarMappable(cmap=my_cmap)
        cbar = plt.colorbar(sm, pad=0.15)
        cbar.set_label('Confidence score', rotation=270, labelpad=25)

        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_yticks(y_pos)
        # ax.set_xlim(-1.8, 0)
        ax.set_yticklabels(labels)
        plt.gca().invert_yaxis()
        plt.xlabel("Scores")
        plt.savefig("rna_scores_plot.png")
        # plt.savefig("{}.png".format(outfile))

    def find_best_rna(self, rrm_seq):

        # To find the best 3-mer RNA for a rrm_seq
        nt_list = [["A", "C", "G", "U"]]
        mer3_list = list(product(*nt_list * 3))
        # Calculate the scores for all the 3-mer RNA fragments

        mer3_scores = {}
        for rna_seq in mer3_list:
            mer3_scores.update(
                self.score_out_seq(rrm_seq, rna_seq=rna_seq,
                                   rna_pos_range=(3, 6)))

        sort_scores = sorted(mer3_scores.items(), key=lambda x: x[1],
                             reverse=True)

        # Add the first and last nucleotides
        rna_ends = list(product(*[["A", "C", "G", "U"], ["A", "C", "G", "U"]]))

        mer5_list = []
        # We take the top-5 3-mers and add all the possible ends to get the 5-mer
        for rna_seq in sort_scores[:5]:
            for rna_end in rna_ends:
                mer5_list.append(rna_end[0] + "".join(rna_seq[0]) + rna_end[1])

        # Calculate the scores for all the 5-mer RNA fragments
        mer5_scores = {}
        for rna_seq in mer5_list:
            mer5_scores.update(
                self.score_out_seq(rrm_seq, rna_seq=rna_seq,
                                   rna_pos_range=(2, 7)))

        sort_scores = sorted(mer5_scores.items(), key=lambda x: x[1],
                             reverse=True)

        # To print the top scoring RNAs
        for rna_seq, score in sort_scores:
            print(rna_seq, score)

    def align_RRM_to_hmm(self, fasta_file):

        fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
        ids_to_align = [fasta.id for fasta in fasta_sequences]

        bashCmd = ["hmmalign", "--mapali",
                   abs_path + "/alignment_data/rrm_bound_domains_aligned_processed_squeezed.fasta",
                   "-o", abs_path + "/output_files/out_test_mapali.stockholm",
                   abs_path + "/alignment_data/rrm_bound.hmm",
                   "{}".format(fasta_file)]

        process = subprocess.Popen(bashCmd, stdout=subprocess.PIPE)
        process.communicate()
        if process.returncode != 0:
            exit("hmmalign failure")

        align = AlignIO.read(abs_path +\
                    "/output_files/out_test_mapali.stockholm", "stockholm")

        seqs_dict = {
            aln_seq.id: str(aln_seq.seq) for aln_seq in align if aln_seq.id in ids_to_align}

        if len(list(seqs_dict.values())[0]) > 243:
            # Get the new mapping manually #TODO: MUST be a better way...
            original_alignment = SeqIO.parse(
                open(abs_path + "/alignment_data/rrm_bound_domains_aligned_processed_squeezed.fasta"), 'fasta')

            orig_seq_to_map = list([str(i.seq).upper() for i in original_alignment][0])
            new_seq_to_map = list([str(i.seq).upper() for i in align][0])

            aln_length = len(orig_seq_to_map)

            new_mapping = {}
            new_gaps = 0
            for i in range(aln_length):
                if new_seq_to_map[:i+1] == orig_seq_to_map[:i+1]:
                    new_mapping[i] = i + new_gaps
                else:
                    while new_seq_to_map[:i+1] != orig_seq_to_map[:i+1]:
                        new_gaps += 1
                        new_seq_to_map.pop(i)
                    new_mapping[i] = i + new_gaps

            # Update the seqs_dict to keep the right size and relevant part
            # IT IS NOT THE REAL SEQUENCE

            for seqid, seq in seqs_dict.items():
                updated_seq = "".join([seq[new_mapping[pos]] for pos in new_mapping.keys()])
                seqs_dict[seqid] = updated_seq

        return seqs_dict
