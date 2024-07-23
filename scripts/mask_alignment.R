
library(argparse)
library(DECIPHER)

parser <- ArgumentParser(description = "Mask highly variable regions of an alignment using DECIPHER's MaskAlignment.")
parser$add_argument("--alignment", help = "Path to the alignment file in fasta format.", required = TRUE)
parser$add_argument("--window_size", help = "Size of the window to calculate the fraction of gaps.", type = "integer", default = 5)
parser$add_argument("--threshold", help = "Threshold to mask a region if the fraction of gaps is above this value.", type = "numeric", default = 1)
parser$add_argument("--max_fraction_gaps", help = "Maximum fraction of gaps allowed in a region to not be masked.", type = "numeric", default = 0.2)
parser$add_argument("--output", help = "Path to the output file in fasta format.", required = TRUE)

args <- parser$parse_args()

# Read alignment as DNAStringSet
alignment <- Biostrings::readDNAStringSet(args$alignment)

# Run MaskAlignment with provided parameters
masked_alignment <- DECIPHER::MaskAlignment(
    alignment,
    type = "sequences",
    windowSize = args$window_size,
    threshold = args$threshold,
    maxFractionGaps = args$max_fraction_gaps,
    showPlot = TRUE
)

# Print number of masked columns as difference between the number of columns in the original and masked alignment
cat("Number of masked columns:", width(alignment[1]) - width(as(masked_alignment, "DNAStringSet")[1]), "\n")

# Display the complete DNA sequence set including the mask
masks <- lapply(width(colmask(masked_alignment)), rep, x="N")
masks <- unlist(lapply(masks, paste, collapse=""))
masked_alignment <- replaceAt(alignment, at=IRanges(colmask(masked_alignment)), value=masks)

# Output has type "DNAMultipleAlignment" which is not directly compatible with Biostrings::writeXStringSet
# Convert to DNAStringSet
masked_alignment <- as(masked_alignment, "DNAStringSet")

# Write masked alignment as .fasta file
Biostrings::writeXStringSet(masked_alignment, file=args$output, format="fasta")
