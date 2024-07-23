### This script evaluates the quality of a multiple sequence alignment of coding regions by counting the number of stop codons and incomplete codons in each sequence

from Bio import SeqIO
import pandas as pd
import argparse

if __name__ == "__main__":

    argparser = argparse.ArgumentParser(description="Evaluate the quality of a multiple sequence alignment of coding regions by counting the number of stop codons and incomplete codons in each sequence")
    argparser.add_argument("--alignments", help="Input alignment files to be evaluated", nargs="+")
    argparser.add_argument("--output", help="Output file (csv) containing the total number of stop codons and incomplete codons in the alignment")

    args = argparser.parse_args()

    alignment_files = args.alignments
    output_file = args.output

    # For each alignment file, for each sequence in the alignment, count the number of stop codons and incomplete codons
    # Output: dataframe with three columns: alignment file, number of sequences, length of the sequences, number of stop codons, number of incomplete codons
    
    results_df = pd.DataFrame(columns=["alignment", "nr_sequences", "sequence_length", "avg_gaps", "stop_codons", "incomplete_codons_ends", "incomplete_codons_internal"])
    
    
    for alignment_file in alignment_files:
        alignment = list(SeqIO.parse(alignment_file, "fasta"))
        
        nr_sequences = len(alignment)
        sequence_length = len(alignment[0].seq)
        stop_codons_total = 0
        incomplete_codons_ends_total = 0
        incomplete_codons_internal_total = 0
        gaps_alignment = 0
          
        for record in alignment:
            # Check if all sequences within the alignment have the same length, otherwise the alignment is not valid
            if len(record.seq) != sequence_length:
                print(f"Sequence {record.id} in alignment {alignment_file} has a different length than the first sequence in the alignment (counting gaps). Make sure the sequences are aligned.")
                break
            
            # Alignments are nucleotide alignments, not AA alignments
            # Stop codons are TGA, TAA, TAG
            # Incomplete codons are sequences that contain 1 or 2 gap characters (e.g., "AT-", "A-G", "--T", "-A-""), but not 3 ("---", that would just be a gap in the sequence)
            # Differentiate between incomplete codons at the start and end of the sequence and in the middle of the sequence
            
            # Walk through the sequence in steps of 3
            stop_codons_seq = 0
            incomplete_codons_ends = 0
            incomplete_codons_internal = 0
            sequence = record.seq
            for i in range(0, len(sequence), 3):
                codon = sequence[i:i+3]
                if codon in ["TGA", "TAA", "TAG", "tga", "taa", "tag"]:
                    stop_codons_seq += 1
                    print(f"Stop codon {codon} found starting at nucleotide position {i} in sequence {record.id} in alignment {alignment_file}.")
                elif "-" in codon and codon.count("-") < 3:
                    # Check whether the incomplete codon is at the start or end of the sequence, i.e., only gaps precede or follow the codon, respectively
                    if sequence[:i].count("-") == i or sequence[i+3:].count("-") == len(sequence) - i - 3:
                        incomplete_codons_ends += 1
                    else:
                        incomplete_codons_internal += 1
                        print(f"Internal incomplete codon {codon} found starting at nucleotide position {i} in sequence {record.id} in alignment {alignment_file}.")
            # Add the counts for the sequence to the alignment counts
            stop_codons_total += stop_codons_seq
            incomplete_codons_ends_total += incomplete_codons_ends
            incomplete_codons_internal_total += incomplete_codons_internal
            
            # Count the number of gaps in the sequence and add it to the alignment count
            gaps_seq = sequence.count("-")
            gaps_alignment += gaps_seq
            
        # Calculate the average number of gaps per sequence
        avg_gaps = round(gaps_alignment / nr_sequences, 2)

        # Add the counts for the alignment to the dataframe using pd.concat
        alignment_name = alignment_file.split("/")[-1] # Get the name of the alignment file
        # Remove .fasta extension if present
        if alignment_name.endswith(".fasta"):
            alignment_name = alignment_name[:-6]
        alignment_results_df = pd.DataFrame({"alignment": [alignment_name], "nr_sequences": nr_sequences, "sequence_length": sequence_length, "avg_gaps": avg_gaps, "stop_codons": [stop_codons_total], "incomplete_codons_ends": [incomplete_codons_ends_total], "incomplete_codons_internal": [incomplete_codons_internal_total]})
        results_df = pd.concat([results_df, alignment_results_df], ignore_index=True)
        
    # Write the results to a csv file
    results_df.to_csv(output_file, index=False)
                    
        
