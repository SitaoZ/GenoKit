# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq

import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.motifs import create, jaspar
from Bio.motifs.logo import calculate_information_content
from io import StringIO
import logomaker
import pandas as pd

class MotifFinder:
    def __init__(self, fasta_file, motif_length=6, min_occurrences=5):
        """
        Initialize the MotifFinder with a FASTA file
        
        Args:
            fasta_file (str): Path to input FASTA file
            motif_length (int): Length of motifs to search for
            min_occurrences (int): Minimum number of occurrences to consider a motif
        """
        self.fasta_file = fasta_file
        self.motif_length = motif_length
        self.min_occurrences = min_occurrences
        self.sequences = self._load_sequences()
        self.motifs = []
        
    def _load_sequences(self):
        """Load sequences from FASTA file"""
        return [str(record.seq) for record in SeqIO.parse(self.fasta_file, "fasta")]
    
    def _find_all_kmers(self):
        """Find all possible k-mers in the sequences"""
        kmers = defaultdict(int)
        for seq in self.sequences:
            for i in range(len(seq) - self.motif_length + 1):
                kmer = seq[i:i+self.motif_length]
                if 'N' not in kmer:  # Skip ambiguous bases
                    kmers[kmer] += 1
        return kmers
    
    def find_motifs(self):
        """Find overrepresented motifs in the sequences"""
        kmers = self._find_all_kmers()
        
        # Filter by minimum occurrences
        filtered_kmers = {k: v for k, v in kmers.items() if v >= self.min_occurrences}
        
        # Convert to BioPython motif format
        for motif_seq, count in filtered_kmers.items():
            instances = [Seq(motif_seq) for _ in range(count)]
            motif = create(instances)
            self.motifs.append({
                'sequence': motif_seq,
                'count': count,
                'motif_object': motif,
                'ic': calculate_information_content(motif)
            })
        
        # Sort by information content
        self.motifs.sort(key=lambda x: x['ic'], reverse=True)
        return self.motifs
    
    def visualize_motif(self, motif_index=0, output_file=None):
        """
        Visualize a motif using sequence logo
        
        Args:
            motif_index (int): Index of motif to visualize
            output_file (str): Optional file path to save the visualization
        """
        if not self.motifs:
            self.find_motifs()
            
        if motif_index >= len(self.motifs):
            raise ValueError(f"Motif index {motif_index} out of range")
            
        motif = self.motifs[motif_index]['motif_object']
        
        # Convert to PWM matrix for logomaker
        pwm = []
        for position in motif.pwm:
            pwm.append({
                'A': position['A'],
                'C': position['C'],
                'G': position['G'],
                'T': position['T']
            })
        df = pd.DataFrame(pwm)
        
        # Create sequence logo
        plt.figure(figsize=(10, 3))
        logo = logomaker.Logo(df, font_name='Arial')
        
        # Style the logo
        logo.style_spines(visible=False)
        logo.style_spines(spines=['left', 'bottom'], visible=True)
        logo.ax.set_ylabel('Information content', labelpad=-1)
        logo.ax.set_xlabel('Position', labelpad=-1)
        plt.title(f"Motif: {self.motifs[motif_index]['sequence']}\n"
                 f"Count: {self.motifs[motif_index]['count']} "
                 f"IC: {self.motifs[motif_index]['ic']:.2f}")
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
    
    def save_motifs(self, output_file):
        """Save found motifs to a file in JASPAR format"""
        if not self.motifs:
            self.find_motifs()
            
        with open(output_file, 'w') as f:
            for i, motif in enumerate(self.motifs):
                jaspar.write(motif['motif_object'], f)
                f.write(f"\n# Motif {i+1}: {motif['sequence']}\n")
                f.write(f"# Count: {motif['count']}\n")
                f.write(f"# Information content: {motif['ic']:.2f}\n\n")

    def get_top_motifs(self, n=5):
        """Get top n motifs by information content"""
        if not self.motifs:
            self.find_motifs()
        return self.motifs[:n]


def get_motif(args):
    """
    Parameters:
    """
    pass
# Example usage
if __name__ == "__main__":
    # Initialize with your FASTA file
    finder = MotifFinder("example.fasta", motif_length=6, min_occurrences=5)
    
    # Find motifs
    motifs = finder.find_motifs()
    print(f"Found {len(motifs)} motifs")
    
    # Visualize the top motif
    finder.visualize_motif(0, output_file="top_motif.png")
    
    # Save all motifs
    finder.save_motifs("found_motifs.jaspar")
    
    # Get top 3 motifs
    top_motifs = finder.get_top_motifs(3)
    for i, motif in enumerate(top_motifs, 1):
        print(f"\nTop motif #{i}:")
        print(f"Sequence: {motif['sequence']}")
        print(f"Count: {motif['count']}")
        print(f"Information content: {motif['ic']:.2f}")
