#!/usr/bin/env python
# coding: utf-8

# In[10]:


import gzip

def analyze_fastq_gz(file_path):
    read_num = 0
    tot_lenght = 0
    gc_tot = 0

    max_len = 150
    sum_score_pos = [0] * max_len
    pos_count = [0] * max_len

    print(f"I'm analyzing the compress file: {file_path}...")

    try:
        # We use 'gzip.open' instead of the regular 'open' because the file is a .gz file
        with gzip.open(file_path, 'rt') as f:
            while True:
                header = f.readline()
                if not header:
                    break

                sequence = f.readline().strip()
                separator = f.readline()
                quality = f.readline().strip()

                read_num += 1

                # 1. Length calculation
                lenght = len(sequence)
                tot_lenght += lenght

                # 2. GC content calculation
                g_num = sequence.count('G') + sequence.count('g')
                c_num = sequence.count('C') + sequence.count('c')
                if lenght > 0:
                    gc_tot += ((g_num + c_num) / lenght) * 100

                # 3.  Conversion of Phred Scores and Quality by Position
                phred_scores = [ord(char) - 33 for char in quality]

                for i, q_score in enumerate(phred_scores):
                    if i < max_len:
                        sum_score_pos[i] += q_score
                        pos_count[i] += 1

            # Print results
            print("\n--- Analysis of the results ---")
            print(f"Total readings analyzed: {read_num}")
            if read_num > 0:
                print(f"Average lenght of readings: {tot_lenght / read_num:.2f} bp")
                print(f"Average GC content: {gc_tot / read_num:.2f}%")

                print("\nAverage quality for the top 10 positions:")
                for i in range(10):
                    if pos_count[i] > 0:
                        average_pos = sum_score_pos[i] / pos_count[i]
                        print(f" Position{i+1}: Q = {average_pos:.2f}")

    except Exception as e:
        print(f"Error during processing: {e}")

if __name__ == "__main__":
    analyze_fastq_gz("1M_SRR9336468_1.fastq.gz")


# Comparison Test with BioPython
import gzip
from Bio import SeqIO

def analyze_with_biopython(file_path):
    read_num = 0
    tot_lenght = 0
    tot_gc = 0

    print(f"Analysis of the file using Biopython: {file_path} ...")
    with gzip.open(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_num += 1
            seq_len = len(record.seq)
            tot_lenght += seq_len

            g_count = record.seq.count("G") + record.seq.count("g")
            c_count = record.seq.count("C") + record.seq.count("c")
            if seq_len > 0:
                tot_gc += ((g_count + c_count) / seq_len) * 100

    print("\n--- Results using Byopython ---")
    print(f"Total number of readings analyzed: {read_num}")
    if read_num > 0: 
        print(f"Average lenght: {tot_lenght / read_num:.2f} bp")
        print(f"Average GC content: {tot_gc / read_num:.2f}%")

if __name__ == "__main__":
    analyze_with_biopython("1M_SRR9336468_1.fastq.gz")

# In[ ]:




