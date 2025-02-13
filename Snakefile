rule files:
    input:
        evc_sequences = "data/sequences/all_evc_sequences.fasta",
        metadata = "data/metadata/all_evc_metadata.csv",
        evc_annotation_gff = "data/gene_annotation/annotation.gff",
        evc_annotation_tsv = "data/gene_annotation/annotation.tsv",
        coding_region_annotation_gff = "data/gene_annotation/annotation_coding_region.gff",
        coding_region_annotation_tsv = "data/gene_annotation/annotation_coding_region.tsv",
        rivm_genotyping_full = "data/metadata/rivm_genotyping_results_fullgenome.csv",
        colors = "data/config/colors.tsv",
        vp1_clades = "data/config/vp1_clades.tsv",
        sites_of_interest = "data/config/sites_of_interest.tsv",
        auspice_config = "data/config/auspice_config.json"

files = rules.files.input

# Prepare data: filter out questionable sequences (from labs, patents, etc.), try to extract serotypes from metadata
# Safe unidentified sequences and metadata in extra files

rule prepare_data:
    input:
        sequences = files.evc_sequences,    
        metadata = files.metadata
    output:
        metadata_filtered = "data/metadata/all_evc_metadata_filtered.csv",
        sequences_filtered = "data/sequences/all_evc_sequences_filtered.fasta",
        metadata_unidentified = "data/metadata/unidentified_metadata.csv",
        sequences_unidentified = "data/sequences/unidentified_sequences.fasta"
    params:
        min_length = 300
    run:
        shell(
        "python scripts/prepare_data.py \
        --sequences {input.sequences} \
        --metadata {input.metadata} \
        --output_sequences {output.sequences_filtered} \
        --output_metadata {output.metadata_filtered} \
        --output_unidentified_sequences {output.sequences_unidentified} \
        --output_unidentified_metadata {output.metadata_unidentified} \
        --min_length {params.min_length}"
        )

rule filter_all_full_genomes:
    input:
        fasta = rules.prepare_data.output.sequences_filtered,
        metadata = rules.prepare_data.output.metadata_filtered
    output:
        filtered_fasta = "data/sequences/all_evc_full_genomes.fasta",
        filtered_metadata = "data/metadata/evc_full_genomes_metadata.csv"
    params:
        min_length = 6000
    run:
        shell(
            "python scripts/filter_fasta_by_length.py \
            --fasta {input.fasta} \
            --metadata {input.metadata} \
            --output_fasta {output.filtered_fasta} \
            --output_metadata {output.filtered_metadata} \
            --min_length {params.min_length}"
            )

rule rivm_genotyping:
    input: 
        metadata = rules.filter_all_full_genomes.output.filtered_metadata,
        genotyping = files.rivm_genotyping_full
    output:
        updated_metadata = "data/metadata/evc_full_genomes_metadata_rivm.csv"
    run:
        shell(
            "python scripts/add_rivm_genotyping.py \
            --metadata {input.metadata} \
            --genotypes {input.genotyping} \
            --output {output.updated_metadata}"
            )

rule plot_metadata_stats:
    input:
        metadata = rules.rivm_genotyping.output.updated_metadata,
    params:
        plot_dir = "plots/metadata_stats"
    run:
        shell(
            "python scripts/plot_metadata_stats.py \
            --metadata {input.metadata} \
            --colors {files.colors} \
            --outdir {params.plot_dir}"
        )

rule divide_fasta_into_ref_and_query:
    input:
        fasta = rules.filter_all_full_genomes.output.filtered_fasta
    output:
        reference_fasta = "data/sequences/ref_full_genomes.fasta",
        query_fasta = "data/sequences/query_full_genomes.fasta"
    params:
        reference_accn = "NC_002058.3"
    run:
        shell(
            "python scripts/divide_fasta_into_ref_and_query.py \
            --fasta {input.fasta} \
            --reference_accn {params.reference_accn} \
            --output_ref {output.reference_fasta} \
            --output_query {output.query_fasta}"
            )


### Try multiple alignment options and compare --> start with full genome alignments 

# Option 1: mafft
rule mafft_full_genomes:
    input:
        fasta_query = rules.divide_fasta_into_ref_and_query.output.query_fasta,
        fasta_ref = rules.divide_fasta_into_ref_and_query.output.reference_fasta
    output:
        aligned_fasta = "data/sequences/alignments/full_genomes/mafft/mafft_aligned_full_genomes.fasta"
    run:
        shell("mafft --auto --keeplength --add {input.fasta_query} {input.fasta_ref} > {output.aligned_fasta}")

# Option 2: nextclade align --> codon-aware alignment (with increased gap penalites)
rule nextalign_full_genomes:
    input:
        fasta_query = rules.divide_fasta_into_ref_and_query.output.query_fasta,
        fasta_ref = rules.divide_fasta_into_ref_and_query.output.reference_fasta,
        gff = files.evc_annotation_gff
    output:
        aligned_fasta = "data/sequences/alignments/full_genome/nextalign/nextalign_full_genomes.fasta"
    params:
        min_seed_cover = 0,
        penalty_gap_open = 24,
        penalty_gap_open_in_frame = 50,
        penalty_gap_open_out_of_frame = 60
    run:
        shell(
            "nextclade run \
            {input.fasta_query} \
            --verbose \
            --min-seed-cover {params.min_seed_cover} \
            --penalty-gap-open {params.penalty_gap_open} \
            --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
            --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
            --input-ref {input.fasta_ref} \
            --include-reference \
            --input-annotation {input.gff} \
            --output-fasta {output.aligned_fasta} \
            --output-translations='data/sequences/alignments/full_genome/nextalign/nextalign_translations_{{cds}}.fasta'"
            )



### Try codon-aware alignments for the coding regions only

## Requirement: separate coding regions from non-coding regions

# First step: run mafft nucleotide alignment keeping insertions
rule mafft_keep_insertions:
    input:
        fasta_query = rules.divide_fasta_into_ref_and_query.output.query_fasta,
        fasta_ref = rules.divide_fasta_into_ref_and_query.output.reference_fasta
    output:
        aligned_fasta = "data/sequences/alignments/full_genome/mafft/mafft_aligned_full_genomes_with_insertions.fasta"
    run:
        shell("mafft --auto --add {input.fasta_query} {input.fasta_ref} > {output.aligned_fasta}")

## Find 5'UTR
# Manually/visually identify start codon in alignment (AliView)
# Are start codons nicely aligned? --> yes, but some sequences start are missing the first few amino acids of the protein coding region
# Cut-off point in the sequence alignment: bp 1099 (A of ATG for Methionin)

## Find 3'UTR
# Search for sequence of the reference genome's 3'UTR in the alignment; begins with stop codon (TAG in reference genome)
# Are last base pairs of the protein coding region nicely aligned? --> yes, but some sequence completely lack the 3'UTR
# Cut-off point in the sequence alignment: bp 7948 (last nucleotide of the protein coding region, T in majority of sequences)

## Use these cut-off points to split off the non-coding regions 

rule split_noncoding_regions:
    input: 
        aligned_fasta = rules.mafft_keep_insertions.output.aligned_fasta,
        reference = files.evc_annotation_gff
    output:
        five_prime_utr = "data/sequences/alignments/noncoding_regions/5UTR.fasta",
        coding_region = "data/sequences/alignments/coding_region/coding_region.fasta",
        three_prime_utr = "data/sequences/alignments/noncoding_regions/3UTR.fasta"
    run:
        shell(
            "python scripts/split_noncoding_regions.py \
            --alignment {input.aligned_fasta} \
            --reference {input.reference} \
            --output_coding_region {output.coding_region} \
            --output_5UTR {output.five_prime_utr} \
            --output_3UTR {output.three_prime_utr}"
            )


## Try multiple alignment options for the coding region only
# Firstly, remove gaps from alignment again
rule ungap_coding_region:
    input:
        coding_region = rules.split_noncoding_regions.output.coding_region
    output:
        coding_region_ungapped = "data/sequences/alignments/coding_region/coding_region_nogaps.fasta"
    run:
        shell("python scripts/remove_gaps.py --input {input.coding_region} --output {output.coding_region_ungapped}")

# Secondly, divide into reference and query again
rule divide_coding_region_into_ref_and_query:
    input:
        fasta = rules.ungap_coding_region.output.coding_region_ungapped
    output:
        reference_fasta = "data/sequences/alignments/coding_region/coding_region_ref.fasta",
        query_fasta = "data/sequences/alignments/coding_region/coding_region_query.fasta"
    params:
        reference_accn = "NC_002058.3"
    run:
        shell(
            "python scripts/divide_fasta_into_ref_and_query.py \
            --fasta {input.fasta} \
            --reference_accn {params.reference_accn} \
            --output_ref {output.reference_fasta} \
            --output_query {output.query_fasta}"
            )

# Option 1: re-run mafft on the coding region only, removing insertions relative to the reference (--keeplength)
rule mafft_coding_region:
    input:
        fasta_query = rules.divide_coding_region_into_ref_and_query.output.query_fasta,
        fasta_ref = rules.divide_coding_region_into_ref_and_query.output.reference_fasta
    output:
        aligned_fasta = "data/sequences/alignments/coding_region/mafft/mafft_coding_region.fasta"
    benchmark: 
        "data/alignment_evaluation/benchmarks/mafft.txt"
    run:
        shell("mafft --auto --keeplength --add {input.fasta_query} {input.fasta_ref} > {output.aligned_fasta}")

# Option 2: nextalign applied to coding regions only
rule nextalign_coding_region:
    input:
        fasta_query = rules.divide_coding_region_into_ref_and_query.output.query_fasta,
        fasta_ref = rules.divide_coding_region_into_ref_and_query.output.reference_fasta,
        gff = files.coding_region_annotation_gff
    output:
        aligned_fasta = "data/sequences/alignments/coding_region/nextalign/nextalign_coding_region.fasta"
    params:
        min_seed_cover = 0,
        penalty_gap_open = 24,
        penalty_gap_open_in_frame = 60,
        penalty_gap_open_out_of_frame = 70
    benchmark:
        "data/alignment_evaluation/benchmarks/nextalign.txt"
    run:
        shell(
            "nextclade run \
            {input.fasta_query} \
            --verbose \
            --min-seed-cover {params.min_seed_cover} \
            --penalty-gap-open {params.penalty_gap_open} \
            --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
            --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
            --input-ref {input.fasta_ref} \
            --include-reference \
            --input-annotation {input.gff} \
            --output-fasta {output.aligned_fasta} \
            --output-translations='data/sequences/alignments/coding_region/nextalign/translation_{{cds}}.fasta' \
            --output-tsv='data/sequences/alignments/coding_region/nextalign/nextalign_output.tsv'"
            )


# Option 3: MACSE v2 (using enrichAlignment mode against the reference sequence, but this probably doesn't matter)
rule macse_align:
    input:
        query = rules.divide_coding_region_into_ref_and_query.output.query_fasta,
        ref = rules.divide_coding_region_into_ref_and_query.output.reference_fasta
    output:
        nuc_alignment = "data/sequences/alignments/coding_region/macse/macse_nucleotides.fasta",
        aa_alignment = "data/sequences/alignments/coding_region/macse/macse_aminoacids.fasta"
    benchmark:
        "data/alignment_evaluation/benchmarks/macse.txt"
    run:
        shell(
            "macse -prog enrichAlignment \
            -align {input.ref} \
            -seq {input.query} \
            -out_AA {output.aa_alignment} \
            -out_NT {output.nuc_alignment}"
            )

# MACSE denotes frameshift mutations (i.e., gaps in the middle of a codon) as an exclamation mark ("!") in the nucleotide alignment
# Convert these exclamation marks to regular gaps ("-") for downstream comparison

rule convert_macse_frameshifts:
    input:
        alignment = rules.macse_align.output.nuc_alignment
    output:
        alignment_converted = "data/sequences/alignments/coding_region/macse/macse_nucleotides_fs_converted.fasta"
    run:
        shell("sed -e 's/!/-/g' {input.alignment} > {output.alignment_converted}")


# Option 4: VIRULIGN
rule virulign:
    input:
        query = rules.divide_coding_region_into_ref_and_query.output.query_fasta,
        ref = rules.divide_coding_region_into_ref_and_query.output.reference_fasta
    output:
        nuc_alignment = "data/sequences/alignments/coding_region/virulign/virulign_nucleotides.fasta"
    benchmark:
        "data/alignment_evaluation/benchmarks/virulign.txt"
    run:
        shell(
            "virulign {input.ref} {input.query} \
            --exportAlphabet Nucleotides \
            --exportKind GlobalAlignment \
            --exportWithInsertions no \
            --exportReferenceSequence yes \
            --progress yes > {output.nuc_alignment}"
            )


# Option 5: kc-align

# kc-align run directly on the fasta produces an error after codon alignment, because the sequences are not all the same length
# Possibly, kc-align seems to expect the sequences to be a multiple of three in length?
# The intermediate output still looks good; the sequences are aligned and the codons are in frame
# Try: remove trailing nucleotides from the sequences and run kc-align again --> this works
# Note: kcalign does not remove remove insertions relative to the reference sequence 

rule remove_trailing_nucleotides:
    input:
        fasta = rules.divide_coding_region_into_ref_and_query.output.query_fasta,
        trailing_starts = "data/sequences/alignments/coding_region/kcalign/trailing_start_nucleotides.json"
    output:
        fasta_trimmed = "data/sequences/alignments/coding_region/kcalign/coding_region_query_trimmed.fasta"
    run:
        shell("python scripts/remove_trailing_nucleotides.py --fasta {input.fasta} --trailing_starts {input.trailing_starts} --output {output.fasta_trimmed}")

rule kcalign:
    input:
        query = rules.remove_trailing_nucleotides.output.fasta_trimmed,
        ref = rules.divide_coding_region_into_ref_and_query.output.reference_fasta,
    output:
        fasta = "data/sequences/alignments/coding_region/kcalign/kcalign_nucleotides.fasta",
        clustal = "data/sequences/alignments/coding_region/kcalign/kcalign_nucleotides.clustal"
    benchmark:
        "data/alignment_evaluation/benchmarks/kcalign.txt"
    run:
        shell("kc-align --mode gene --reference {input.ref} --sequences {input.query}")
        shell("mv codon_align.fasta data/sequences/alignments/coding_region/kcalign/kcalign_nucleotides.fasta")  # kc-align produces files in the current directory, there is no option for specifying output location --> move
        shell("mv codon_align.clustal data/sequences/alignments/coding_region/kcalign/kcalign_nucleotides.clustal")
        

# The kc-align script sometimes seems to have issues: no messages are printed to the console and the output is full of gaps
# Rerunning the exact same command sometimes produces a good alignment
# Possibly happens when I don't delete the clustal file before rerunning the command?
# Also: kcalign throws out sequences suspected of containing frameshits, early stop codons, or having too many Ns
# While this is probably OK, I would rather keep all the sequences at this stage
# For some reason, the reference sequence appears twice in the output, although it's not in the query file; remove the duplicate (it's the first sequence in the file)

rule remove_kcalign_duplicates:
    input:
        alignment = rules.kcalign.output.fasta
    output:
        alignment_new = "data/sequences/alignments/coding_region/kcalign/kcalign_nucleotides_no_duplicates.fasta"
    run:
        shell("python scripts/remove_duplicate_sequences_from_alignment.py --alignment {input.alignment} --output {output.alignment_new}")
        

# For MACSE and kc-align, insertions relative to the reference are not removed; do this in an extra step for downstream analysis and comparability

rule remove_insertions_relative_to_reference:
    input:
        kcalign = rules.remove_kcalign_duplicates.output.alignment_new,
        macse = rules.convert_macse_frameshifts.output.alignment_converted
    output:
        kcalign_trimmed = "data/sequences/alignments/coding_region/kcalign/kcalign_nucleotides_trimmed.fasta",
        macse_trimmed = "data/sequences/alignments/coding_region/macse/macse_nucleotides_trimmed.fasta"
    params:
        ref = "NC_002058.3"
    run:
        shell("python scripts/remove_insertions_relative_to_reference.py --alignment {input.kcalign} --reference {params.ref} --output {output.kcalign_trimmed}")
        shell("python scripts/remove_insertions_relative_to_reference.py --alignment {input.macse} --reference {params.ref} --output {output.macse_trimmed}")

# Gather codon-aware alignment files for comparison
def get_coding_region_alignments(wildcards):
    return [
        rules.mafft_coding_region.output.aligned_fasta,
        rules.virulign.output.nuc_alignment,
        rules.remove_insertions_relative_to_reference.output.macse_trimmed,
        rules.remove_insertions_relative_to_reference.output.kcalign_trimmed,
        rules.nextalign_coding_region.output.aligned_fasta
    ]

# Compare the results of the different codon-aware alignment tools by computing number of incomplete codons and stop codons
rule evaluate_codon_alignments:
    input:
        alignments = get_coding_region_alignments
    output:
        csv = "data/alignment_evaluation/alignment_evaluation.csv",
        plot = "plots/alignment_evaluation.png"
    run:
        shell("python scripts/evaluate_codon_alignments.py --alignments {input.alignments} --output {output.csv} --plot {output.plot}")

# Plot the number of gaps and mutations across the alignments
rule plot_gaps_mutations_pre_masking:
    input:
        alignments = get_coding_region_alignments,
        gene_annotation = files.coding_region_annotation_tsv
    params:
        reference_accn = "NC_002058.3",
        output_csv = "data/alignment_evaluation/gaps_mutations_pre_masking"
    output:
        output_plot = "plots/gaps_mutations_pre_masking/gaps_mutations_pre_masking_combined.png"
    run:
        shell(
            "python scripts/plot_gaps_mutations_across_alignment.py \
            --alignments {input.alignments} \
            --reference_accn {params.reference_accn} \
            --gene_annotation {input.gene_annotation} \
            --output_csv {params.output_csv} \
            --output_plot {output.output_plot}"
            )

# Choose which alignment to move forward with --> Nextalign
# Reorder sequence names in the alignment: move accession number to the back of the name, to make the sequences sortable by serotype in alignment viewer
rule reorder_sequence_names:
    input:
        alignment = rules.nextalign_coding_region.output.aligned_fasta
    output:
        alignment_reordered_names = "data/sequences/alignments/coding_region/nextalign/nextalign_coding_region_reordered_names.fasta"
    run:
        shell("python scripts/reorder_sequence_names.py --input {input.alignment} --output {output.alignment_reordered_names}")

# Mask coding region alignment using MaskAlignment function implemented in DECIPHER R to hide badly aligned regions
rule mask_alignment:
    input:
        alignment = rules.nextalign_coding_region.output.aligned_fasta
    output:
        masked_alignment = "data/sequences/alignments/coding_region/nextalign/nextalign_coding_region_masked_ws15.fasta",
        plot = "plots/nextalign_coding_region_masking_ws15.pdf"
    params:
        window_size = 15,
        threshold = 1,
        max_fraction_gaps = 0.4
    run:
        shell("Rscript scripts/mask_alignment.R --alignment {input.alignment} --output {output.masked_alignment} --window_size {params.window_size} --threshold {params.threshold} --max_fraction_gaps {params.max_fraction_gaps}")
        shell("mv Rplots.pdf {output.plot}")

# Get percentage of masked sites per gene
rule get_masking_stats:
    input:
        masked_alignment = rules.mask_alignment.output.masked_alignment,
        gene_annotation = files.coding_region_annotation_tsv
    params:
        mask_char = "+"
    output:
        csv = "data/alignment_evaluation/masking_stats.csv",
        plot = "plots/masking_stats.png"
    run:
        shell("python scripts/get_masking_stats.py --alignment {input.masked_alignment} --mask_char {params.mask_char} --gene_annotation {input.gene_annotation} --outcsv {output.csv} --outplot {output.plot}")

# Re-plot the number of gaps and mutations across the alignments
rule plot_gaps_mutations_post_masking:
    input:
        alignments = rules.mask_alignment.output.masked_alignment,
        gene_annotation = files.coding_region_annotation_tsv
    params:
        reference_accn = "NC_002058.3",
        output_csv = "data/alignment_evaluation/gaps_mutations_post_masking"
    output:
        output_plot = "plots/gaps_mutations_post_masking/gaps_mutations_post_masking_nextalign.png"
    run:
        shell(
            "python scripts/plot_gaps_mutations_across_alignment.py \
            --alignments {input.alignments} \
            --reference_accn {params.reference_accn} \
            --gene_annotation {input.gene_annotation} \
            --output_csv {params.output_csv} \
            --output_plot {output.output_plot}"
            )

genes = ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]

# Split the alignment by gene
rule split_coding_region_by_gene:
    input:
        alignment = rules.mask_alignment.output.masked_alignment,
        gene_annotation = files.coding_region_annotation_tsv
    output:
        expand("data/sequences/alignments/genes/{gene}_aligned.fasta", gene=genes)
    run:
        shell("python scripts/split_alignment_by_gene.py --alignment {input.alignment} --gene_annotation {input.gene_annotation} --output_dir data/sequences/alignments/genes")


## Re-align non-coding regions with mafft

# Again, ungap and split into query and reference sequences
rule ungap_noncoding_regions:
    input:
        five_prime = rules.split_noncoding_regions.output.five_prime_utr,
        three_prime = rules.split_noncoding_regions.output.three_prime_utr
    output:
        five_prime_ungapped = "data/sequences/alignments/noncoding_regions/5UTR_ungapped.fasta",
        three_prime_ungapped = "data/sequences/alignments/noncoding_regions/3UTR_ungapped.fasta"
    run:
        shell("python scripts/remove_gaps.py --input {input.five_prime} --output {output.five_prime_ungapped}")
        shell("python scripts/remove_gaps.py --input {input.three_prime} --output {output.three_prime_ungapped}")

rule divide_noncoding_regions_into_ref_and_query:
    input:
        five_prime = rules.ungap_noncoding_regions.output.five_prime_ungapped,
        three_prime = rules.ungap_noncoding_regions.output.three_prime_ungapped
    output:
        five_prime_ref = "data/sequences/alignments/noncoding_regions/5UTR_ref.fasta",
        five_prime_query = "data/sequences/alignments/noncoding_regions/5UTR_query.fasta",
        three_prime_ref = "data/sequences/alignments/noncoding_regions/3UTR_ref.fasta",
        three_prime_query = "data/sequences/alignments/noncoding_regions/3UTR_query.fasta"
    params:
        reference_accn = "NC_002058.3"
    run:
        shell("python scripts/divide_fasta_into_ref_and_query.py --fasta {input.five_prime} --reference_accn {params.reference_accn} --output_ref {output.five_prime_ref} --output_query {output.five_prime_query}")
        shell("python scripts/divide_fasta_into_ref_and_query.py --fasta {input.three_prime} --reference_accn {params.reference_accn} --output_ref {output.three_prime_ref} --output_query {output.three_prime_query}")

rule mafft_noncoding_regions:
    input:
        five_prime_query = rules.divide_noncoding_regions_into_ref_and_query.output.five_prime_query,
        five_prime_ref = rules.divide_noncoding_regions_into_ref_and_query.output.five_prime_ref,
        three_prime_query = rules.divide_noncoding_regions_into_ref_and_query.output.three_prime_query,
        three_prime_ref = rules.divide_noncoding_regions_into_ref_and_query.output.three_prime_ref
    output:
        five_prime_aligned = "data/sequences/alignments/noncoding_regions/5UTR_mafft_aligned.fasta",
        three_prime_aligned = "data/sequences/alignments/noncoding_regions/3UTR_mafft_aligned.fasta"
    run:
        shell("mafft --auto --keeplength --add {input.five_prime_query} {input.five_prime_ref} > {output.five_prime_aligned}")
        shell("mafft --auto --keeplength --add {input.three_prime_query} {input.three_prime_ref} > {output.three_prime_aligned}")

rule assemble_full_genomes:
    input:
        five_prime = rules.mafft_noncoding_regions.output.five_prime_aligned,
        coding_region = rules.mask_alignment.output.masked_alignment,
        three_prime = rules.mafft_noncoding_regions.output.three_prime_aligned
    output:
        full_genome = "data/sequences/alignments/full_genome/assembled/fullgenome_assembled.fasta"
    run:
        shell("python scripts/assemble_full_genomes.py --five_prime {input.five_prime} --coding_region {input.coding_region} --three_prime {input.three_prime} --output {output.full_genome}")

genomic_regions = genes + ["5UTR", "3UTR", "fullgenome"]

# Collect all alignments for downstream analysis
all_alignments = rules.split_coding_region_by_gene.output + [rules.mafft_noncoding_regions.output.five_prime_aligned, rules.mafft_noncoding_regions.output.three_prime_aligned, rules.assemble_full_genomes.output.full_genome]

# Filter out sequences that contain more than max_gap_fraction gaps
rule filter_gappy_seqs:
    input:
        alignments = all_alignments
    output:
        expand("data/sequences/alignments/filtered/{region}_filtered.fasta", region=genomic_regions)
    params:
        max_gap_fraction_gene = 0.1,
        max_gap_fraction_utr = 0.5   # Allow more gaps in the UTRs, as they seem more variable (e.g., variable sequencing starts and ends in the 5'UTR and 3'UTR)
    run:
        shell("python scripts/remove_gappy_seqs_from_alignments.py --alignments {input.alignments} --max_gap_fraction_gene {params.max_gap_fraction_gene} --max_gap_fraction_utr {params.max_gap_fraction_utr} --output_dir data/sequences/alignments/filtered")

# Build trees
rule trees:
    input:
        alignments = rules.filter_gappy_seqs.output
    output:
        trees = expand("data/trees/{region}/{region}_tree_raw.nwk", region=genomic_regions)
    run:
        for region in genomic_regions:
            alignment = [aln for aln in input.alignments if region in aln][0]
            tree = f"data/trees/{region}/{region}_tree_raw.nwk"
            print(f"\n\nGenerating tree for {region} ...\n")
            shell(f"augur tree --alignment {alignment} --output {tree}")

# I don't want to build a time-resolved tree because VDPVs are independently seeded from the same source (OPV)
# Still run augur refine to add node names for downstream steps ('ancestral', 'translate')

rule refine:
    input:
        trees = rules.trees.output.trees
    output:
        trees = expand("data/trees/{region}/{region}_tree.nwk", region=genomic_regions),
        node_data = expand("data/trees/{region}/{region}_branch_lengths.json", region=genomic_regions)
    params:
        root = "mid_point"
    run:
        for region in genomic_regions:
            tree = [t for t in input.trees if region in t][0]
            output_tree = f"data/trees/{region}/{region}_tree.nwk"
            output_node_data = f"data/trees/{region}/{region}_branch_lengths.json"
            shell(
                f"augur refine \
                --tree {tree} \
                --root {params.root} \
                --output-tree {output_tree} \
                --output-node-data {output_node_data} \
                --stochastic-resolve"
                )

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        trees = rules.refine.output.trees,
        alignments = rules.filter_gappy_seqs.output
    output:
        nt_data = expand("data/trees/{region}/{region}_nt_muts.json", region=genomic_regions)
    params:
        inference = "joint"
    run:
        # Match trees and alignments by region and generate node data files
        for region in genomic_regions:
            tree = [t for t in input.trees if region in t][0]
            alignment = [aln for aln in input.alignments if region in aln][0]
            
            node_data = f"data/trees/{region}/{region}_nt_muts.json"
            print(f"\n\nReconstructing ancestral sequences and mutations for {region} ...\nTree: {tree}\nAlignment: {alignment}\n")
            shell(
                f"augur ancestral \
                --tree {tree} \
                --alignment {alignment} \
                --output-node-data {node_data} \
                --inference {params.inference} \
                --keep-ambiguous"
                )

# Translate trees for genes and full genome (exclude UTRs)
trees_to_translate = expand("data/trees/{region}/{region}_tree.nwk", region=genes) + ["data/trees/fullgenome/fullgenome_tree.nwk"]
node_data_to_translate = expand("data/trees/{region}/{region}_nt_muts.json", region=genes) + ["data/trees/fullgenome/fullgenome_nt_muts.json"]
references = expand("data/config/references/{region}_reference.gb", region=genes) + ["data/config/references/fullgenome_reference.gb"]

rule translate:
    message: "Translating nucleotide to amino acid sequences"
    input:
        tree = trees_to_translate,
        node_data = node_data_to_translate,
        reference = references
    output:
        aa_data = expand("data/trees/{region}/{region}_aa_muts.json", region=genes) + ["data/trees/fullgenome/fullgenome_aa_muts.json"]
    run:
        for region in genes + ["fullgenome"]:
            tree = [t for t in input.tree if region in t][0]
            node_data = [n for n in input.node_data if region in n][0]
            reference = [r for r in input.reference if region in r][0]
            output_node_data = f"data/trees/{region}/{region}_aa_muts.json"
            print(f"\n\nTranslating nucleotide to amino acid sequences for {region} ...\nTree: {tree}\nNode data: {node_data}\nReference: {reference}\n")
            shell(
                f"augur translate \
                --tree {tree} \
                --ancestral-sequences {node_data} \
                --reference-sequence {reference} \
                --output-node-data {output_node_data}"
                )

# Export data for visualization in auspice
rule export_all:
    message: "Exporting data files for visualization in auspice"
    input:
        trees = rules.refine.output.trees,
        metadata = rules.rivm_genotyping.output.updated_metadata,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.nt_data,
        aa_muts = rules.translate.output.aa_data,
        colors = files.colors,
        auspice_config = files.auspice_config,
        attenuation_determinants = files.sites_of_interest
    output:
        auspice_jsons = expand("data/auspice_jsons/{region}.json", region=genomic_regions)
    params:
        metadata_id_col = "Accession"
    run:
        for region in genomic_regions:
            print(f"\n\nExporting data files for {region} ...\n")
            tree = [t for t in input.trees if region in t][0]
            branch_lengths = [b for b in input.branch_lengths if region in b][0]
            nt_muts = [n for n in input.nt_muts if region in n][0]
            output_json = f"data/auspice_jsons/{region}.json"
            if "UTR" in region:
                # Export UTRs without translation
                shell(
                    f"augur export v2 \
                    --tree {tree} \
                    --metadata {input.metadata} \
                    --metadata-id-columns {params.metadata_id_col} \
                    --node-data {branch_lengths} {nt_muts} \
                    --auspice-config {input.auspice_config} \
                    --colors {input.colors} \
                    --output {output_json}"
                    )
            else:
                # Export full genome and gene trees with translation
                aa_muts = [a for a in input.aa_muts if region in a][0]
                shell(
                    f"augur export v2 \
                    --tree {tree} \
                    --metadata {input.metadata} \
                    --metadata-id-columns {params.metadata_id_col} \
                    --node-data {branch_lengths} {nt_muts} {aa_muts} \
                    --auspice-config {input.auspice_config} \
                    --colors {input.colors} \
                    --output {output_json}"
                    )

# Update strain names in auspice JSONS to include both accession number and GenBank title
# Save in new directory and adjust file names for nextstrain upload
rule strain_names:
    input:
        auspice_jsons = rules.export_all.output.auspice_jsons,
        metadata = rules.rivm_genotyping.output.updated_metadata
    output:
        updated_auspice_jsons = expand("data/nextstrain_build/enterovirus_c_{region}.json", region=genomic_regions)
    params:
        accn_col = "Accession",
        gbtitle_col = "GenBank_Title"
    run:
        for region in genomic_regions:
            auspice_json = [a for a in input.auspice_jsons if region in a][0]
            output_json = f"data/nextstrain_build/enterovirus_c_{region}.json"
            print(f"\nUpdating strain names for {region} ...")
            shell(
                f"python scripts/replace_auspice_id.py \
                --auspice_json {auspice_json} \
                --metadata {input.metadata} \
                --accn_col {params.accn_col} \
                --gbtitle_col {params.gbtitle_col} \
                --output_auspice_json {output_json}"
                )

serotypes = ["CVA1", "CVA22", "EVC113", "EVC116", "CVA11", "CVA13", "CVA17", "CVA19", "CVA20", "CVA21", "CVA24", "EVC102", "EVC104", "EVC105", "EVC109", "EVC117", "EVC118", "EVC95", "EVC96", "EVC99", "PV1", "PV2", "PV3"]

# Separate sequences by serotype for p-distance calculation
rule separate_by_serotype:
    input:
        metadata = rules.rivm_genotyping.output.updated_metadata,
        alignments = rules.filter_gappy_seqs.output
    output:
        alignments = expand("data/sequences/alignments/by_serotype/{serotype}/{serotype}_{region}_aligned.fasta", region=genomic_regions, serotype=serotypes)
    params:
        metadata_id_col = "Accession",
        metadata_type_col = "serotype_short"
    run:
        shell(
            "python scripts/split_sequences_by_serotype.py \
            --alignments {input.alignments} \
            --metadata {input.metadata} \
            --metadata_id_col {params.metadata_id_col} \
            --metadata_type_col {params.metadata_type_col} \
            --outdir data/sequences/alignments/by_serotype"
            )

# Calculate p-distances for each cluster and region
rule calculate_pdistances:
    input:
        alignments = rules.separate_by_serotype.output.alignments
    output:
        pdistance_csvs = expand("data/pdistances/by_serotype/{serotype}_pdistances.csv", serotype=serotypes)
    run:
        for serotype in serotypes:
            alignment_files = [aln for aln in input.alignments if serotype in aln]
            alignment_files = " ".join(alignment_files)
            print(f"\n\nProcessing serotype {serotype} ...")
            output = f"data/pdistances/by_serotype/{serotype}_pdistances.csv"
            shell(f"python scripts/calculate_pdistances.py --alignments {alignment_files} --output {output}")

# Make interactive p-distance plots using plotly and dash
rule plotly_pdistances:
    input:
        pdistance_csvs = rules.calculate_pdistances.output.pdistance_csvs,
        metadata = rules.rivm_genotyping.output.updated_metadata
    params:
        regions = genomic_regions,
        cutoff = 0.1
    run: 
        shell(
            "python scripts/plotly_pdistances.py \
            --pdistances {input.pdistance_csvs} \
            --metadata {input.metadata} \
            --regions {params.regions} \
            --cutoff {params.cutoff}"
            )

# Make non-interactive p-distance plots to save
rule plot_pdistances:
    input:
        pdistance_csvs = rules.calculate_pdistances.output.pdistance_csvs,
        metadata = rules.rivm_genotyping.output.updated_metadata
    params:
        serotypes = ["PV1", "PV2", "PV3", "EVC104"],
        reference_region = "VP1",
        compare_regions = ["VP2", "2A", "3D", "5UTR"],
        outdir = "plots/pdistances",
        outformat = "png"
    run:
        shell(
            "python scripts/plot_pdistances.py \
            --pdistances {input.pdistance_csvs} \
            --serotypes {params.serotypes} \
            --reference_region {params.reference_region} \
            --compare_regions {params.compare_regions} \
            --grid \
            --outdir {params.outdir} \
            --outformat {params.outformat}"
            )

# Extract candidate recombinants based on p-distance cutoffs
rule extract_candidate_recombinants:
    input:
        pdistance_csvs = rules.calculate_pdistances.output.pdistance_csvs,
    params:
        regions = ["5UTR", "2A", "2C", "3D"],
        cutoff = 0.1,
        plotdir = "plots/candidate_identification"
    output:
        candidates_df = "data/pdistances/candidate_identification/candidate_recombinants.csv"
    run:
        shell(
            "python scripts/extract_candidate_recombinants.py \
            --pdistances {input.pdistance_csvs} \
            --regions {params.regions} \
            --cutoff {params.cutoff} \
            --output {output.candidates_df} \
            --plotdir {params.plotdir}"
            )



### Build consensus sequences for different groups of sequences (all, non-polio, serotypes) for downstream analysis

# Assemble pre-masked alignment (combining UTRs and coding region)
rule pre_masked_alignment:
    input:
        five_prime = rules.mafft_noncoding_regions.output.five_prime_aligned,
        three_prime = rules.mafft_noncoding_regions.output.three_prime_aligned,
        coding_region = rules.nextalign_coding_region.output.aligned_fasta
    output:
        pre_masked_alignment = "data/sequences/alignments/full_genome/assembled/fullgenome_pre_masking_assembled.fasta"
    run:
        shell("python scripts/assemble_full_genomes.py --five_prime {input.five_prime} --coding_region {input.coding_region} --three_prime {input.three_prime} --output {output.pre_masked_alignment}")

# Extract sequence groups of interest from the alignment (non-polio sequences, serotypes)
rule extract_subalignments:
    input:
        alignment = rules.pre_masked_alignment.output.pre_masked_alignment,
        metadata = rules.rivm_genotyping.output.updated_metadata
    params:
        metadata_id_col = "Accession",
        metadata_group_col = "serotype_short"
    output:
        output_all = "data/sequences/alignments/full_genome/for_consensus/all.fasta",
        output_non_polio = "data/sequences/alignments/full_genome/for_consensus/non_polio.fasta",
        output_serotypes = expand("data/sequences/alignments/full_genome/for_consensus/{serotype}.fasta", serotype=serotypes)
    run:
        outdir = "data/sequences/alignments/full_genome/for_consensus"
        shell(
            "python scripts/extract_subalignments.py \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata_id_col {params.metadata_id_col} \
            --metadata_group_col {params.metadata_group_col} \
            --outdir {outdir}"
            )

rule consensus_sequences:
    input:
        alignments = rules.extract_subalignments.output
    output:
        consensus_sequences = "data/sequences/consensus_sequences/all_consensus_sequences.fasta"
    run:
        shell("python scripts/build_consensus_sequences.py --input {input.alignments} --output {output}")

# Slice sequences into n bp chunk and calculate p-distances against consensus sequence of each serotype
# Output dataframe structure: seq, serotype, chunk, consensus, p-distance
rule consensus_distances:
    input:
        consensus_sequences = rules.consensus_sequences.output.consensus_sequences,
        alignment = rules.pre_masked_alignment.output.pre_masked_alignment,
        metadata = rules.rivm_genotyping.output.updated_metadata
    params:
        chunk_size = 250,
        metadata_id_col = "Accession",
        metadata_group_col = "serotype_short"
    output:
        distances = "data/pdistances/candidate_identification/consensus_distances.csv"
    run:
        shell(
            "python scripts/calculate_consensus_distances.py \
            --consensus {input.consensus_sequences} \
            --alignment {input.alignment} \
            --chunk_size {params.chunk_size} \
            --metadata {input.metadata} \
            --metadata_id_col {params.metadata_id_col} \
            --metadata_group_col {params.metadata_group_col} \
            --output {output.distances}"
            )



# Identify closest neighbors of each candidate recombinant in the phylogenetic tree
# For finding recombination partners in similarity plots

rule find_neighbors:
    input:
        trees = rules.refine.output.trees,
        metadata = rules.rivm_genotyping.output.updated_metadata
    params:
        n_neighbors = 3,
        metadata_id_col = "Accession",
        metadata_type_col = "serotype_short"
    output:
        neighbors = "data/recombination_analysis/custom_simplots/neighbors_all.csv"
    run:
        shell(
            "python scripts/find_neighbors.py \
            --trees {input.trees} \
            --metadata {input.metadata} \
            --metadata_id_col {params.metadata_id_col} \
            --metadata_type_col {params.metadata_type_col} \
            --n_neighbors {params.n_neighbors} \
            --out {output.neighbors}"
            )

# Filter neighbors to include only those within a certain distance threshold
rule filter_neighbors:
    input:
        neighbors = rules.find_neighbors.output.neighbors
    params:
        distance_threshold = 0.15
    output:
        filtered_neighbors = "data/recombination_analysis/custom_simplots/neighbors_filtered_all.csv"
    run:
        shell(
            "python scripts/filter_neighbors.py \
            --neighbors {input.neighbors} \
            --distance_threshold {params.distance_threshold} \
            --out {output.filtered_neighbors}"
            )

# Generate alignments containing the candidate recombinant, its closest neighbors, the consensus sequences, and the Sabin sequences
rule neighbor_alignments:
    input:
        neighbors = rules.filter_neighbors.output.filtered_neighbors,
        alignment = rules.pre_masked_alignment.output.pre_masked_alignment,
        consensus_sequences = rules.consensus_sequences.output.consensus_sequences
    params:
        candidate_col = "seq",
        neighbor_col = "neighbor_accn",
        outdir = "data/recombination_analysis/custom_simplots/alignments",
        sabin_accessions = ["AY184219.1", "AY184220.1", "AY184221.1"]
    run:
        shell(
            "python scripts/generate_simplot_alignments.py \
            --alignment {input.alignment} \
            --neighbors {input.neighbors} \
            --neighbors_candidate_col {params.candidate_col} \
            --neighbors_neighbor_col {params.neighbor_col} \
            --consensus {input.consensus_sequences} \
            --sabin {params.sabin_accessions} \
            --outdir {params.outdir}"
            )

# Make similarity plots for the candidate recombinants
rule custom_simplots:
    input:
        alignment_dir = rules.neighbor_alignments.params.outdir,
        metadata = rules.rivm_genotyping.output.updated_metadata,
        colors = files.colors
    params:
        windowsize = 200,
        stepsize = 50,
        #query_seqs = ["AB180071.1", "MG212478.1", "KR259355.1", "MG212483.1", "MZ171089.1", "JX275065.2", "HQ738286.1", "MG212429.1", "HQ738288.1"],  # optional: make similartity plots only for subset of sequences
        outdir = "plots/simplots/candidates",
        outformat = "png"
    run:
        # If query sequences are provided, generate similarity plots only for these sequences
        query_seqs = getattr(params, "query_seqs", None)
        if query_seqs:
            query_alignments = [input.alignment_dir + "/" + accn + ".fasta" for accn in query_seqs]
        else:
            # Get all alignment files in the directory
            import os
            query_alignments = [os.path.join(input.alignment_dir, f) for f in os.listdir(input.alignment_dir) if f.endswith(".fasta")]
        for alignment in query_alignments:
            print(f"\n\nGenerating similarity plot for {alignment} ...\n")
            shell(
                "python scripts/custom_simplots.py \
                --alignment {alignment} \
                --metadata {input.metadata} \
                --colors {input.colors} \
                --windowsize {params.windowsize} \
                --stepsize {params.stepsize} \
                --sabin \
                --outdir {params.outdir} \
                --outformat {params.outformat}"
                )

# Extend SimPlot analysis to include all sequences in the alignment
rule custom_simplots_extended:
    input:
        alignment = rules.pre_masked_alignment.output.pre_masked_alignment,   # full genome alignment, pre-masked
        metadata = rules.rivm_genotyping.output.updated_metadata,
        colors = files.colors,
        neighbors = rules.filter_neighbors.output.filtered_neighbors,
        consensus_sequences = rules.consensus_sequences.output.consensus_sequences
    params:
        windowsize = 200,
        stepsize = 50,
        outdir_plot = "plots/simplots/all",
        outformat = "png",
        outdir_csv = "data/recombination_analysis/custom_simplots/results"
    run:
        shell(
            "python scripts/custom_simplots_extended.py \
            --alignment {input.alignment} \
            --neighbors {input.neighbors} \
            --consensus {input.consensus_sequences} \
            --metadata {input.metadata} \
            --colors {input.colors} \
            --windowsize {params.windowsize} \
            --stepsize {params.stepsize} \
            --outdir_plot {params.outdir_plot} \
            --outformat {params.outformat} \
            --outdir_csv {params.outdir_csv}"
            )


# After sorting the simplots, assign recombination status of each sequence (recombinant, non-recombinant, unclear) to the metadata
rule assign_recombination_status:
    input:
        metadata = rules.rivm_genotyping.output.updated_metadata
    output:
        updated_metadata = "data/metadata/evc_full_genomes_metadata_recombination_status.csv"
    run:
        shell(
            "python scripts/assign_recombination_status.py \
            --metadata_in {input.metadata} \
            --metadata_out {output.updated_metadata}"
            )


# Run VirusRecom

# Prepare input files for VirusRecom: 
# 1. Alignment file (.fasta): each sequence name requires a lineage/serotype mark
# 2. Text file of lineage/serotype marks (.txt); simple listing, one per line

rule prepare_virusrecom:
    input:
        alignment = rules.pre_masked_alignment.output.pre_masked_alignment,
        metadata = rules.rivm_genotyping.output.updated_metadata
    params:
        metadata_id_col = "Accession",
        metadata_group_col = "serotype_short"
    output:
        output_alignment = "data/recombination_analysis/virusrecom/input/virusrecom_input_alignment.fasta",
        output_reference_txt = "data/recombination_analysis/virusrecom/input/reference_names.txt"
    run:
        shell(
            "python scripts/prepare_virusrecom_input.py \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata_id_col {params.metadata_id_col} \
            --metadata_group_col {params.metadata_group_col} \
            --output_alignment {output.output_alignment} \
            --output_reference_txt {output.output_reference_txt}"
            )

# Run VirusRecom for all sequences
rule run_virusrecom:
    input:
        formatted_alignment = rules.prepare_virusrecom.output.output_alignment,
        reference_names = rules.prepare_virusrecom.output.output_reference_txt
    params:
        method = "a",                          # "p" to use polymorphic sites only ("a" for all sites)
        gaps = "y",                            # "n" to exclude sites containing gaps ("y" to include)
        window_size = 200,
        step_size = 50,
        max_region = 6000,                     # maximum allowed length of recombinant region; if "-m p": maximum number of polymorphic sites contained in a recombinant region
        percentage = 0.9,                      # cutoff threshold of proportion (cp, default: 0.9) used for searching recombination regions when mWIC/EIC >= cp
        breakpoint_scan = "y",                 # "y" to perform breakpoint scan ("n" to skip)
        bw = 100                               # window size for breakpoint scan (step size of breakpoint search is fixed to 1)
    output:
        outdir = directory("data/recombination_analysis/virusrecom/results")
    run:
        # Consider each sequence in the alignment a potential recombinant sequence
        # Loop through all sequences in the alignment
        from Bio import SeqIO
        alignment = list(SeqIO.parse(f"{input.formatted_alignment}", "fasta"))
        for seq in alignment:
            print(f"\n\nRunning VirusRecom for {seq.description} ...\n")
            # Create folder for each recombinant sequence
            outdir = f"{output.outdir}/{seq.description}"
            shell(
                "virusrecom \
                -a {input.formatted_alignment} \
                -q {seq.id} \
                -l {input.reference_names} \
                -g {params.gaps} \
                -m {params.method} \
                -w {params.window_size} \
                -s {params.step_size} \
                -mr {params.max_region} \
                -cp {params.percentage} \
                -b {params.breakpoint_scan} \
                -bw {params.bw} \
                -o {outdir}"
                )


# Identify breakpoints in the recombinant sequences based on similarity data
# Take only the sequences identified as recombinent

recombinant_ids = os.listdir("plots/simplots/all_sorted/recombinant")
recombinant_ids = [r.split("_")[1] for r in recombinant_ids]

rule identify_breakpoints:
    input:
        similarity_data = expand("data/recombination_analysis/custom_simplots/results/{seq}_similarity_results.csv", seq=recombinant_ids)
    params:
        highest_similarity_threshold = 0.98,
        high_similarity_threshold = 0.95,
        low_similarity_threshold = 0.88,
        high_difference_threshold = 0.08,
        low_difference_threshold = 0.02
    output:
        inferred_breakpoints = "data/recombination_analysis/custom_simplots/inferred_breakpoints.csv"
    run:
        shell(
            "python scripts/identify_breakpoints.py \
            --similarity_data {input.similarity_data} \
            --highest_similarity_threshold {params.highest_similarity_threshold} \
            --high_similarity_threshold {params.high_similarity_threshold} \
            --low_similarity_threshold {params.low_similarity_threshold} \
            --high_difference_threshold {params.high_difference_threshold} \
            --low_difference_threshold {params.low_difference_threshold} \
            --output {output.inferred_breakpoints}"
            )

# Replot simplots integrating the inferred breakpoints (as vertical lines)
rule replot_simplots:
    input:
        similarity_data = expand("data/recombination_analysis/custom_simplots/results/{seq}_similarity_results.csv", seq=recombinant_ids),
        breakpoints = rules.identify_breakpoints.output.inferred_breakpoints
    params:
        outdir_plot = "plots/simplots/with_breakpoints",
        outformat = "png"
    run:
        shell(
            "python scripts/replot_simplots.py \
            --similarity_data {input.similarity_data} \
            --breakpoints {input.breakpoints} \
            --outdir_plot {params.outdir_plot} \
            --outformat {params.outformat}"
            )

