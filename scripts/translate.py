### This script translates DNA sequences in a fasta file to an amino acid sequence

from Bio import SeqIO

# Read input fasta
fasta = SeqIO.parse(snakemake.input.coding_regions_fasta, "fasta")

# Translate the sequences
def translate(fasta):
    translated_seqs = []
    for record in fasta:

        # Remove gaps
        record.seq = record.seq.replace("-", "")
        
        # Check if sequences are a multiple of 3
        if len(record.seq) % 3 != 0:
            print(f"Sequence {record.id} is not a multiple of 3")
            
        # Remove incomplete codons (trailing bases) from the start
        # Create dictionary of sequence IDs and number of nucleotides that need to be removed from the start
        # Five sequences start with a gap, four of those have nucleotides that should be removed:
        trailing_bases = {
            "OK570202.1": 2,
            "OK570231.1": 1,
            "OK570237.1": 1,
            "MF990297.1": 1
        }
        
        for seq_id, num_bases in trailing_bases.items():
            if record.id == seq_id:
                record.seq = record.seq[num_bases:]
        
        # Assumption: the remaining sequences that are not a multiple of 3 arise due to extra nucleotides at the end
        # Therefore, remove the extra nucleotides from the end
        record.seq = record.seq[:len(record.seq) - len(record.seq) % 3]
        
        # Translate
        record.seq = record.seq.translate()
        translated_seqs.append(record)
        
    return translated_seqs
    
translated_coding_regions = translate(fasta)

# Write the translated sequences to a file
SeqIO.write(translated_coding_regions, snakemake.output.translated_coding_regions, "fasta")
