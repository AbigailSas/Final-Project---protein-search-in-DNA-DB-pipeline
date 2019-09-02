import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict

IRRELEVANT_LINES = 15
INDEX_OF_SEQ_DIRECTION = 43
MONTH_ABBREVIATION = 3
EVAL_LOC = 4
SEQ_NAME_LOC = 1
SCORE_LOC = 2
MONTH_TO_NUMBER_DICT = {"jan":1,"feb":2,"mar":3,"apr":4,"may":5,"jun":6,"jul":7,"aug":8,
                                       "sep":9,"oct":10,"nov":11,"dec":12}

RPS2 = "small subunit ribosomal protein S2 RPS2 consensus\n"
RPL1 = "large subunit ribosomal protein L1 consensus\n"
IF_2 = "translation initiation factor IF-2 consensus\n"
RPL22 = "large subunit ribosomal protein L22\n"

class Parser:

    def __init__(self, directory, e_val_threshold):
        self.directory = directory
        self.e_val_threshold = e_val_threshold


    def parse(self):
        """
        Read all the files in directory - chlF files and chlA files are parsed differently. Keep the number
        of detected chls and return as 2 dicts.
        Files in directory are BLAST result files, from blasting ChlG and PsbA4+PsbA1 against Dor reservoir
        Assume directory is a directory with 2 sub directories called chlA and chlF. In each one there
        are 23 files with names in the format of "dec3_chlf.tx
        E value is the parameter to choose which sequences are a good enough match to be considered chlA
        :return: 2 dictionaries, one for chlF and one for chlA where the key is the month (as int) and the value
        is a list of the values of number of chl found
        """
        # initiate dicts with lists
        chlf_results, chla_results, single_copy_results = defaultdict(list), defaultdict(list), defaultdict(list)
        for file_ in os.listdir(self.directory):
            full_path = self.directory + "\\" + file_
            for file in os.listdir(full_path):
                month_int_repr = MONTH_TO_NUMBER_DICT[file[:MONTH_ABBREVIATION]]
                year_identifier = int(file[MONTH_ABBREVIATION])
                if file_ == "ChlA":
                    chla_results[month_int_repr].append((year_identifier, self.read_file_chl_A(full_path + "\\" + file)))
                elif file_ == "ChlF":
                    chlf_results[month_int_repr].append((year_identifier,
                                                         self.read_file_chl_f(full_path + "\\" + file)))
                elif file_ == "singleCopy":
                    single_copy_results[month_int_repr].append\
                        ((year_identifier,
                          self.read_file_single_copy_genes(full_path + "\\" + file)))
        return chlf_results, chla_results, single_copy_results


    def read_file_chl_A(self,file_name):
        """
        reads file, and counts dna seqs that match chlA with an e value smaller or equal to e_value param
        :param file_name: file name in the format starting with <month><number>
        :param e_value: maximal e value score to be considered as a good hit for homology matches to chlA
        :return: number of sequences found
        """
        f = open(file_name)
        unique = []
        for line in f.readlines()[IRRELEVANT_LINES:]:
            if "NB501373" in line and ">" not in line:
                seq_line_features = line.split("  ")
                dna_seq = seq_line_features[SEQ_NAME_LOC]
                e_val = float(seq_line_features[EVAL_LOC][:-1])
                if (dna_seq not in unique) and (e_val <= self.e_val_threshold):
                    unique.append(dna_seq)
        # identify paired end reads
        unique = set([t[:INDEX_OF_SEQ_DIRECTION] for t in unique])
        return len(unique)

    def read_file_chl_f(self,file_name):
        """
        reads file and counts number of chlF reads, not including paired end reads and reads that are PsbA1
        :param file_name: file name in the format starting with <month><number>
        :param e_value: maximal e value score to be considered as a good hit for homology matches to chlA
        :return: number of sequences found
        """
        f = open(file_name)
        prot_name = ""
        results = defaultdict(list)
        for line in f.readlines()[IRRELEVANT_LINES:]:
            if line.startswith("Query= "):
                prot_name = line[7:]
            if "NB501373" in line and ">" not in line:
                seq_line_features = line.split("  ")
                dna_seq = seq_line_features[SEQ_NAME_LOC]
                score = float(seq_line_features[SCORE_LOC])
                e_val = float(seq_line_features[EVAL_LOC])
                if prot_name != "":
                    results[dna_seq].append([prot_name, score, e_val])

        # decide if read is chlf or psba1
        final = []
        for read in results:
            chlf_vals, psba_vals, f_e_vals = [], [], []
            for read_name, read_score, read_e_val in results[read]:
                if "chlorophyll f" in read_name:
                    chlf_vals.append(read_score)
                    f_e_vals.append(read_e_val)
                elif "Photosystem II" in read_name:
                    psba_vals.append(read_score)
            chlf_average, psba_average, e_val_average = np.average(chlf_vals), np.average(psba_vals), np.average(f_e_vals)
            if chlf_average > psba_average and e_val_average<=float(1e-10):
                final.append(read)

        # identify paired end reads
        unique = set([t[:INDEX_OF_SEQ_DIRECTION] for t in final])
        return len(unique)

    def read_file_single_copy_genes(self, file_name):
        f = open(file_name)
        unique = {RPS2: [], RPL1:[], IF_2:[]}
        gene_name = ""
        results = defaultdict(list)
        for line in f.readlines()[IRRELEVANT_LINES:]:
            if line.startswith("Query= "):
                gene_name = line[7:]
            if "NB501373" in line and ">" not in line and gene_name!=RPL22:
                seq_line_features = line.split("  ")
                dna_seq = seq_line_features[SEQ_NAME_LOC]
                e_val = float(seq_line_features[EVAL_LOC][:-1])
                if (dna_seq not in unique[gene_name]) and (e_val <= self.e_val_threshold):
                    unique[gene_name].append(dna_seq)
        # identify paired end reads
        for gene in unique:
            unique[gene] = len(set([t[:INDEX_OF_SEQ_DIRECTION] for t in unique[gene]]))

        return unique

