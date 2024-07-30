### This script moves accession numbers in the name of a sequence (seperated by "|") to the back of a sequence's name to make them sortable by serotype in an alignment viewer.

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Move accession numbers to the back of a sequence's name to make them sortable by name.")
    parser.add_argument("--input", help="Input fasta file.")
    parser.add_argument("--output", help="Output fasta file.")
    args = parser.parse_args()

    with open(args.input, "r") as f:
        lines = f.readlines()
        

    with open(args.output, "w") as f:
        for line in lines:
            if line.startswith(">"):            
                parts = line[1:].split(" |")  # Remove the initial ">"
                accession = parts[0]
                rest = " |".join(parts[1:]).strip()
                f.write(f">{rest}|{accession}\n")
            else: 
                f.write(line)
                