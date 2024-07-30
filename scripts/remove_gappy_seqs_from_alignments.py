### This script removes sequences with more than the given proportion of gaps/Ns from the alignment.

from Bio import SeqIO
import argparse

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Remove sequences with a proportion of gaps higher than max_gap_fraction from the alignment.")
    argparser.add_argument("--alignments", help="Input alignment file(s) in fasta format.", nargs="+")
    argparser.add_argument("--max_gap_fraction_gene", help="Maximum fraction of gaps allowed in a gene sequence.", type=float)
    argparser.add_argument("--max_gap_fraction_utr", help="Maximum fraction of gaps allowed in the untranslated regions (5'UTR and 3'UTR).", type=float)
    argparser.add_argument("--output_dir", help="Output directory for the filtered alignments.")
    
    args = argparser.parse_args()
    
    # Read the alignment
    alignments = {}
    for alignment_file in args.alignments:
        name = alignment_file.split("/")[-1].split("_")[0]
        alignment = list(SeqIO.parse(alignment_file, "fasta"))
        # Add the alignment to the dictionary
        alignments[name] = alignment

    # Filter alignments based on the maximum gap fraction
    
    # Distinguish between gene and UTR sequences    
    filtered_alignments = {}
    for name, alignment in alignments.items():
        if "UTR" in name:
            max_gap_fraction = args.max_gap_fraction_utr
        else:
            max_gap_fraction = args.max_gap_fraction_gene
        
        # Filter the alignment based on Ns and gaps
        filtered_alignment = []
        for record in alignment:
            sequence = str(record.seq)
            
            # Convert ambiguous nucleotides (R, Y, etc.) to N
            valid_chars = set("ACGTNatcgn-+")  # "+" denotes masked positions --> convert to N after filtering
            sequence = "".join([char if char in valid_chars else "N" for char in sequence])
            
            gap_count = sequence.count("-") + sequence.count("N") + sequence.count("n")
            gap_fraction = gap_count / len(sequence)
            if gap_fraction <= max_gap_fraction:
                # Convert "+" to "N" for masked positions
                record.seq = sequence.replace("+", "N")
                filtered_alignment.append(record)
            else:
                print(f"Removed {record.id} from {name} alignment due to gap fraction of {gap_fraction} (max_gap_fraction: {max_gap_fraction}).")
                
        filtered_alignments[name] = filtered_alignment
    
    # Write the filtered alignments to the output directory
    for name, alignment in filtered_alignments.items():
        output_file = f"{args.output_dir}/{name}_filtered.fasta"
        SeqIO.write(alignment, output_file, "fasta")