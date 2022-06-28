# Python wrapper to have available all the scoring classes in one place and
# merge the results

import sys
import argparse
from rrm_rna_functions import RNAScoring


class Manager():
    # This is the general input manager of the scoring framework
    def __init__(self, usr_input):
        rna_scoring.__init__()
        self.usr_input = usr_input

    def input_handler(self):
        if self.usr_input.fasta_file:
            seqs_dict = rna_scoring.align_RRM_to_hmm(fasta_file=self.usr_input.fasta_file)

        elif self.usr_input.RRM_id:
            seqs_dict = {self.usr_input.RRM_id:
                             rna_scoring.prot_seqs_dict_full[self.usr_input.RRM_id]}

        for seq_id, seq in seqs_dict.items():
            if self.usr_input.top:
                rna_scoring.find_best_rna(rrm_seq=seq)

            elif self.usr_input.rna_seq:
                rna_scoring.score_out_seq(
                    rrm_seq=seq, rna_seq=self.usr_input.rna_seq,
                    rna_pos_range=self.usr_input.rna_pos_range)

                if self.usr_input.plot:
                    rna_scoring.plot_rna_kde(rna_seq=self.usr_input.rna_seq,
                                      scores_dict=rna_scoring.scores_dict,
                                      window_size=self.usr_input.window_size)
                else:
                    print(seq_id)
                    for key, score in rna_scoring.scores_dict.items():
                        print(key, score)

class UserInput():
    def __init__(self):
                # Input files parser
        args = sys.argv[1:]
        parser = argparse.ArgumentParser(
            description='RRM-RNA scoring framework')

        input_arg = parser.add_mutually_exclusive_group(required=True)
        input_arg.add_argument('-RRM', '--RRM_identifier',
                            help='')
        input_arg.add_argument('-fasta', '--fasta_file',
                            help='')

        feat_arg = parser.add_mutually_exclusive_group(required=True)
        feat_arg.add_argument('-RNA', '--RNA_sequence',
                            help='')
        feat_arg.add_argument('-top', '--top_RNAs', action="store_true",
                            help='')

        parser.add_argument('-ws', '--window_size', required=True,
                            help='')
        parser.add_argument('-plot', '--plot_RNAs', action="store_true",
                            help='')

        self.input_files = parser.parse_args(args)

        self.fasta_file = self.input_files.fasta_file
        self.RRM_id = self.input_files.RRM_identifier

        self.rna_seq = self.input_files.RNA_sequence
        self.top = self.input_files.top_RNAs
        self.plot = self.input_files.plot_RNAs

        self.window_size = int(self.input_files.window_size)

        if self.window_size == 3:
            self.rna_pos_range = (3, 6)

        elif self.window_size == 5:
            self.rna_pos_range = (2, 7)


if __name__ == '__main__':
    rna_scoring = RNAScoring()
    usr_input = UserInput()
    manager = Manager(usr_input)
    manager.input_handler()


