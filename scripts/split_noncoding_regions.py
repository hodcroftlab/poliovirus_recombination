
### This script takes a fasta file of a full genome sequence alignment as an input and cuts-out the 5'UTR and 3'UTR
### Input: a fasta file of the full genome sequence alignment (KEEPING insertions)
### Outputs: 1) a fasta file of the 5'UTR; 2) for the 3'UTR sequences; 3) for the coding regions

from Bio import SeqIO
from BCBio import GFF
import argparse

if __name__ == "__main__":
    
    # Parse arguments
    argparser = argparse.ArgumentParser(description="Split the sequences into 5'UTR, 3'UTR, and coding regions.")
    argparser.add_argument("--alignment", help="Aligned fasta file of full genome sequences.")
    argparser.add_argument("--reference", help="Gene annotation file for the reference sequence in .gff format.")
    argparser.add_argument("--output_coding_region", help="Output file for the coding region.")
    argparser.add_argument("--output_5UTR", help="Output file for the 5'UTR.")
    argparser.add_argument("--output_3UTR", help="Output file for the 3'UTR.")
    
    args = argparser.parse_args()
    

    # Read reference .gff file, extract reference sequence ID and CDS positions

    def extract_cds_positions(gff_file):
        first_cds_start = None
        last_cds_end = None

        with open(gff_file) as in_handle:
            for rec in GFF.parse(in_handle):
                reference_id = rec.id
                for feature in rec.features:
                    if feature.type == "CDS":
                        start = feature.location.start + 1  # Convert to 1-based index
                        end = feature.location.end
                        
                        # Check for first CDS
                        if first_cds_start is None or start < first_cds_start:
                            first_cds_start = start
                        
                        # Check for last CDS
                        if last_cds_end is None or end > last_cds_end:
                            last_cds_end = end
        
        return reference_id, first_cds_start, last_cds_end
    
    reference_id, first_cds_start, last_cds_end = extract_cds_positions(args.reference)
    
    # Read input fasta
    alignment = list(SeqIO.parse(args.alignment, "fasta"))
    
    # Get reference sequence
    aligned_reference_seq = None
    for record in alignment:
        if record.id == reference_id:
            aligned_reference_seq = record.seq
            break
        
    if aligned_reference_seq is None:
        raise ValueError(f"Reference sequence with ID '{reference_id}' not found in the alignment.")
        
    # Map the original start and end positions of the reference sequence to its aligned counterpart
    def map_positions_to_alignment(aligned_seq, original_start, original_end):
        count = 0
        aligned_start = None
        aligned_end = None
        
        for i, base in enumerate(aligned_seq):
            if base != "-":
                count += 1
                if count == original_start:
                    aligned_start = i
                if count == original_end:
                    aligned_end = i
                    break

        return aligned_start, aligned_end
    
    cutoff_5utr, cutoff_3utr = map_positions_to_alignment(aligned_reference_seq, first_cds_start, last_cds_end)
    
    # Split the sequences
    five_utr_seqs = []
    three_utr_seqs = []
    coding_region_seqs = []
        
    for record in alignment:
        # 5'UTR
        five_utr = record[:cutoff_5utr]
        five_utr_seqs.append(five_utr)
            
        # 3'UTR
        three_utr = record[cutoff_3utr+1:]
        three_utr_seqs.append(three_utr)
            
        # Coding region
        coding_region = record[cutoff_5utr:cutoff_3utr+1]
        coding_region_seqs.append(coding_region)
            
    # Write the sequences to separate files
    SeqIO.write(five_utr_seqs, args.output_5UTR, "fasta")
    SeqIO.write(coding_region_seqs, args.output_coding_region, "fasta")
    SeqIO.write(three_utr_seqs, args.output_3UTR, "fasta")
            