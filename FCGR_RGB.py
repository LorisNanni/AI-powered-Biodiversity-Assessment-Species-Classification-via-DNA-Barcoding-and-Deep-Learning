"""
Created on Mon Mar 27 14:56:14 2023

@author: daniela
"""
import collections
from matplotlib import pyplot as plt
import math
import cv2
import numpy as np
import os
from Bio import SeqIO


def base_color(base):
    if base == "A":
        return np.array([1, 0, 0])  # Red
    elif base == "C":
        return np.array([0, 1, 0])  # Green
    elif base == "G":
        return np.array([0, 0, 1])  # Blue
    elif base == "T":
        return np.array([1, 1, 0])  # Yellow
    else:
        return np.array([0, 0, 0])  # Black for other characters
def custom_color(kmer):
    base_counts = collections.Counter(kmer)
    total_bases = len(kmer)
    color = np.zeros(3)
    
    for base, count in base_counts.items():
        color += (count / total_bases) * base_color(base)
    
    return color

def count_kmers(sequence, k, data):
    d = collections.defaultdict(int)
    for i in range(len(data)-(k-1)):
        d[sequence[i:i+k]] +=1
    keys_to_remove = []
    for key in d.keys():
        if "N" in key:
            keys_to_remove.append(key)
    for key in keys_to_remove:
        del d[key]
    return d

 
def probabilities(kmer_count, k, data):
    probabilities = collections.defaultdict(float)
    N = len(data)
    for key, value in kmer_count.items():
        probabilities[key] = float(value) / (N - k + 1)
    return probabilities
 
def chaos_game_representation(probabilities, k):
    array_size = int(math.sqrt(4**k))
    chaos = np.zeros((array_size, array_size, 3))
 
    maxx = array_size
    maxy = array_size
    posx = 1
    posy = 1
    for key, value in probabilities.items():
        for char in key:
            if char == "T":
                posx += maxx / 2
            elif char == "C":
                posy += maxy / 2
            elif char == "G":
                posx += maxx / 2
                posy += maxy / 2
            maxx = maxx / 2
            maxy /= 2
        
        # Assegna il colore personalizzato al k-mer
        color = custom_color(key)
        chaos[int(posy-1), int(posx-1)] = np.array(color) * value
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
    
    return chaos
 
# Function to process all records in a FASTA file and save them with numbered filenames
def process_fasta_file(fasta_file, output_folder, k=2):
    sequence_number = 1  # Start numbering from 1
    for record in SeqIO.parse(fasta_file, "fasta"):
        data = str(record.seq)
        kmer_count = count_kmers(data, k, data)
        kmer_prob = probabilities(kmer_count, k, data)
        chaos = chaos_game_representation(kmer_prob, k)
        
        # Resize the image
        chaos_resized = cv2.resize(chaos, (227, 227), interpolation=cv2.INTER_NEAREST)
        
        # Normalize the color values to range [0, 1]
        chaos_normalized = chaos_resized / np.max(chaos_resized)
        
        # Save the image with a numbered filename
        output_path = os.path.join(output_folder, f"{sequence_number}.png")
        plt.imsave(output_path, chaos_normalized, cmap='viridis')
        #plt.imsave(output_path, chaos_resized, cmap=cm.gray_r)
        sequence_number += 1  # Increment the counter for the next sequence

# Process the single FASTA file with multiple entries
input_base = "/Users/dc3000/Desktop/PhD/projects/other/vec2mat/forBLOG/NonAlign/"  

files = ["ATest.fas", "ATrain.fas", "batsTest.fas", "batsTrain.fas", "CypraeidaeTest.fas", "CypraeidaeTrain.fas", "Drosophilatest.fas", "Drosophilatrain.fas", "fishesTest.fas", "fishesTrain.fas", "IngaTest.fas", "IngaTrain.fas"]


k_values = [5,6]
for k in k_values:
    output_base = "FCGR_output/output_k=" + str(k) + "/"
    for file in files:
        fasta_file =  input_base + file
        output_folder = output_base + file 
        os.makedirs(output_folder, exist_ok=True)
        process_fasta_file(fasta_file, output_folder, k=k)