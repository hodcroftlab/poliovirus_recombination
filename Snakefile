rule all:
    input:
        "data/blast_results/blast_results_polio_filtered.csv",
        "data/blast_results/blast_results_npevc_filtered.csv",
        "data/blast_results/blast_results_polio_unique_accn.csv",
        "data/blast_results/blast_results_npevc_unique_accn.csv",
        "data/blast_results/blast_results_polio_unique_submitter.csv",
        "data/blast_results/blast_results_npevc_unique_submitter.csv"

rule files:
    input:
        evc_sequences = "data/sequences/all_evc_sequences.fasta",
        metadata = "data/metadata/all_evc_metadata.csv",
        serotype_corrections = "data/metadata/serotype_corrections.json",
        sabin123_sequences = "data/sequences/sabin123_refseq.fasta",
        sabin123_metadata = "data/metadata/sabin123_metadata.csv",
        evc_annotation_gff = "data/gene_annotation/annotation.gff",
        coding_region_annotation_gff = "data/gene_annotation/annotation_coding_region.gff",
        coding_region_annotation_csv = "data/gene_annotation/annotation_coding_region.csv",
        rivm_genotyping = "data/metadata/rivm_genotyping_tool_results.csv"

files = rules.files.input

# Prepare data: filter out questionable sequences (from labs, patents, etc.), try to extract serotypes from metadata
# Safe unidentified sequences and metadata in extra files for RIVM Enterovirus Genotyping Tool

rule prepare_data:
    input:
        sequences = files.evc_sequences,
        metadata = files.metadata,
        serotype_corrections = files.serotype_corrections
    output:
        metadata_filtered = "data/metadata/all_evc_metadata_filtered.csv",
        sequences_filtered = "data/sequences/all_evc_sequences_filtered.fasta",
        metadata_unidentified = "data/metadata/unidentified_metadata.csv",
        sequences_unidentified = "data/sequences/unidentified_sequences.fasta"
    params:
        min_length = 300
    shell:
        """
        python scripts/prepare_data.py \
        --sequences {input.sequences} \
        --metadata {input.metadata} \
        --serotype_corrections {input.serotype_corrections} \
        --output_sequences {output.sequences_filtered} \
        --output_metadata {output.metadata_filtered} \
        --output_unidentified_sequences {output.sequences_unidentified} \
        --output_unidentified_metadata {output.metadata_unidentified} \
        --min_length {params.min_length}
        """

rule split_data:
    input:
        sequences = rules.prepare_data.output.sequences_filtered,
        metadata = rules.prepare_data.output.metadata_filtered,
        rivm_genotyping = files.rivm_genotyping
    output:
        polio_sequences = "data/sequences/polio_sequences.fasta",
        npevc_sequences = "data/sequences/npevc_sequences.fasta",
        polio_metadata = "data/metadata/polio_metadata.csv",
        npevc_metadata = "data/metadata/npevc_metadata.csv"
    shell:
        """
        python scripts/split_data.py \
        --sequences {input.sequences} \
        --metadata {input.metadata} \
        --genotyping {input.rivm_genotyping} \
        --output_polio_metadata {output.polio_metadata} \
        --output_npevc_metadata {output.npevc_metadata} \
        --output_polio_sequences {output.polio_sequences} \
        --output_npevc_sequences {output.npevc_sequences}
        """

rule plot_metadata_stats:
    input:
        metadata_polio = rules.split_data.output.polio_metadata,
        metadata_npevc = rules.split_data.output.npevc_metadata
    params:
        plot_dir = "plots"
    run:
        shell(
            "python scripts/plot_metadata_stats.py \
            --metadata_polio {input.metadata_polio} \
            --metadata_npevc {input.metadata_npevc} \
            --output_dir {params.plot_dir}"
        )


rule plot_seqlen_hist:
    input:
        fasta_polio = rules.split_data.output.polio_sequences,
        fasta_npevc = rules.split_data.output.npevc_sequences
    output:
        hist_polio = "plots/polio_seqlen_hist.pdf",
        hist_npevc = "plots/npevc_seqlen_hist.pdf"
    params:
        binwidth = 200,
        xstep = 400
    script:
        "scripts/plot_seqlen_hist.py"

rule filter_polio:
    input:
        fasta = rules.split_data.output.polio_sequences
    output:
        filtered_fasta = "data/sequences/polio_sequences_filtered.fasta"
    params:
        min_length = 2000
    script:
        "scripts/filter_fasta_by_length.py"

rule filter_npevc:
    input:
        fasta = rules.split_data.output.npevc_sequences
    output:
        filtered_fasta = "data/sequences/npevc_sequences_filtered.fasta"
    params:
        min_length = 300
    script:
        "scripts/filter_fasta_by_length.py"

rule plot_seqlen_hist_filtered:
    input:
        fasta_polio = rules.filter_polio.output.filtered_fasta,
        fasta_npevc = rules.filter_npevc.output.filtered_fasta
    output:
        hist_polio = "plots/polio_seqlen_hist_filtered.pdf",
        hist_npevc = "plots/npevc_seqlen_hist_filtered.pdf"
    params:
        binwidth = 200,
        xstep = 400
    script:
        "scripts/plot_seqlen_hist.py"

rule blast_polio:
    input: 
        query = rules.filter_polio.output.filtered_fasta,
        reference = rules.filter_npevc.output.filtered_fasta
    output:
        blast_out = "data/blast_results/blast_results_polio_raw.csv"
    params:
        evalue = 1e-6
    run:
        shell("makeblastdb -in {input.reference} -dbtype nucl -out data/blast_dbs/npevc_db/npevc_db")
        shell("blastn -query {input.query} -db data/blast_dbs/npevc_db/npevc_db -out {output.blast_out} -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -evalue {params.evalue} -num_threads 4")

rule blast_npevc:
    input: 
        query = rules.filter_npevc.output.filtered_fasta,
        reference = files.sabin123_sequences
    output:
        blast_out = "data/blast_results/blast_results_npevc_raw.csv"
    params:
        evalue = 1e-6
    run:
        shell("makeblastdb -in {input.reference} -dbtype nucl -out data/blast_dbs/sabin123_db/sabin123_db")
        shell("blastn -query {input.query} -db data/blast_dbs/sabin123_db/sabin123_db -out {output.blast_out} -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -evalue {params.evalue} -num_threads 4")

rule blast_filter:
    input:
        blast_results_polio = rules.blast_polio.output.blast_out,
        blast_results_npevc = rules.blast_npevc.output.blast_out,
        polio_metadata = rules.split_data.output.polio_metadata,
        npevc_metadata = rules.split_data.output.npevc_metadata,
        sabin123_metadata = files.sabin123_metadata
    output:
        filtered_blast_polio = "data/blast_results/blast_results_polio_filtered.csv",
        filtered_blast_npevc = "data/blast_results/blast_results_npevc_filtered.csv",
        unique_accn_blast_polio = "data/blast_results/blast_results_polio_unique_accn.csv",
        unique_accn_blast_npevc = "data/blast_results/blast_results_npevc_unique_accn.csv",
        unique_submitter_blast_polio = "data/blast_results/blast_results_polio_unique_submitter.csv",
        unique_submitter_blast_npevc = "data/blast_results/blast_results_npevc_unique_submitter.csv"
    params:
        min_sequence_identity = 85,
        min_match_length = 300
    script:
        "scripts/blast_filter.py"

rule filter_all_full_genomes:
    input:
        fasta = rules.prepare_data.output.sequences_filtered
    output:
        filtered_fasta = "data/sequences/all_evc_full_genomes.fasta"
    params:
        min_length = 6000
    script:
        "scripts/filter_fasta_by_length.py"

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
        aligned_fasta = rules.mafft_keep_insertions.output.aligned_fasta
    output:
        five_prime_utr = "data/sequences/alignments/noncoding_regions/5UTR.fasta",
        coding_region = "data/sequences/alignments/coding_region/coding_region.fasta",
        three_prime_utr = "data/sequences/alignments/noncoding_regions/3UTR.fasta"
    params:
        cutoff_5UTR = 1099,
        cutoff_3UTR = 7948
    run:
        shell(
            "python scripts/split_noncoding_regions.py \
            --alignment {input.aligned_fasta} \
            --cutoff_5UTR {params.cutoff_5UTR} \
            --cutoff_3UTR {params.cutoff_3UTR} \
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
        shell("python scripts/remove_gaps.py --input {input.fasta} --output {output.fasta_nogaps}")

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
# The output still looks good; the sequences are aligned and the codons are in frame
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
        rules.nextalign_coding_region.output.aligned_fasta,
        rules.remove_insertions_relative_to_reference.output.macse_trimmed,
        rules.virulign.output.nuc_alignment,
        rules.remove_insertions_relative_to_reference.output.kcalign_trimmed
    ]

# Compare the results of the different codon-aware alignment tools by computing number of incomplete codons and stop codons
rule evaluate_codon_alignments:
    input:
        alignments = get_coding_region_alignments
    output:
        evaluation = "data/alignment_evaluation/alignment_evaluation.csv"
    run:
        shell("python scripts/evaluate_codon_alignments.py --alignments {input.alignments} --output {output.evaluation}")

# Plot the number of gaps and mutations across the alignments
rule plot_gaps_mutations_pre_masking:
    input:
        alignments = get_coding_region_alignments,
        gene_annotation = files.coding_region_annotation_csv
    params:
        reference_accn = "NC_002058.3",
        output_csv = "data/alignment_evaluation/gaps_mutations_pre_masking",
        output_plot = "plots/gaps_mutations_pre_masking"
    run:
        shell(
            "python scripts/plot_gaps_mutations_across_alignment.py \
            --alignments {input.alignments} \
            --reference_accn {params.reference_accn} \
            --gene_annotation {input.gene_annotation} \
            --output_csv {params.output_csv} \
            --output_plot {params.output_plot}"
            )

# Choose which alignment to move forward with

# Reorder sequence names in the alignment: move accession number to the back of the name, to make the sequences sortable by serotype in the alignment viewer
rule reorder_sequence_names:
    input:
        alignment = rules.nextalign_coding_region.output.aligned_fasta
    output:
        alignment_reordered_names = "data/sequences/alignments/coding_region/nextalign/nextalign_coding_region_reordered_names.fasta"
    run:
        shell("python scripts/reorder_sequence_names.py --input {input.alignment} --output {output.alignment_reordered_names}")

# At this point, probably want to mask regions that seem to be badly aligned (gappy, ambigously aligned regions)

# Mask coding region alignment using MaskAlignment function implemented in DECIPHER R function to remove badly aligned regions
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

# Re-plot the number of gaps and mutations across the alignments
rule plot_gaps_mutations_post_masking:
    input:
        alignments = rules.mask_alignment.output.masked_alignment,
        gene_annotation = files.coding_region_annotation_csv
    params:
        reference_accn = "NC_002058.3",
        output_csv = "data/alignment_evaluation/gaps_mutations_post_masking",
        output_plot = "plots/gaps_mutations_post_masking"
    run:
        shell(
            "python scripts/plot_gaps_mutations_across_alignment.py \
            --alignments {input.alignments} \
            --reference_accn {params.reference_accn} \
            --gene_annotation {input.gene_annotation} \
            --output_csv {params.output_csv} \
            --output_plot {params.output_plot}"
            )

# Split the alignment by gene
rule split_coding_region_by_gene:
    input:
        alignment = rules.mask_alignment.output.masked_alignment,
        gene_annotation = files.coding_region_annotation_csv
    output:
        protein_vp4 = "data/sequences/alignments/genes/VP4_aligned.fasta",
        protein_vp2 = "data/sequences/alignments/genes/VP2_aligned.fasta",
        protein_vp3 = "data/sequences/alignments/genes/VP3_aligned.fasta",
        protein_vp1 = "data/sequences/alignments/genes/VP1_aligned.fasta",
        protein_2a = "data/sequences/alignments/genes/2A_aligned.fasta",
        protein_2b = "data/sequences/alignments/genes/2B_aligned.fasta",
        protein_2c = "data/sequences/alignments/genes/2C_aligned.fasta",
        protein_3a = "data/sequences/alignments/genes/3A_aligned.fasta",
        protein_3b = "data/sequences/alignments/genes/3B_aligned.fasta",
        protein_3c = "data/sequences/alignments/genes/3C_aligned.fasta",
        protein_3d = "data/sequences/alignments/genes/3Dpol_aligned.fasta",
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

# Collect all alignments for downstream analysis
def get_alignments(wildcards):
    return [
        rules.split_coding_region_by_gene.output.protein_vp4,
        rules.split_coding_region_by_gene.output.protein_vp2,
        rules.split_coding_region_by_gene.output.protein_vp3,
        rules.split_coding_region_by_gene.output.protein_vp1,
        rules.split_coding_region_by_gene.output.protein_2a,
        rules.split_coding_region_by_gene.output.protein_2b,
        rules.split_coding_region_by_gene.output.protein_2c,
        rules.split_coding_region_by_gene.output.protein_3a,
        rules.split_coding_region_by_gene.output.protein_3b,
        rules.split_coding_region_by_gene.output.protein_3c,
        rules.split_coding_region_by_gene.output.protein_3d,
        rules.mafft_noncoding_regions.output.five_prime_aligned,
        rules.mafft_noncoding_regions.output.three_prime_aligned
    ]

# Filter UTRs and coding regions for sequences that contain more than max_gap_fraction gaps
rule filter_gappy_seqs:
    input:
        alignments = get_alignments
    output:
        protein_vp4 = "data/sequences/alignments/filtered/VP4_filtered.fasta",
        protein_vp2 = "data/sequences/alignments/filtered/VP2_filtered.fasta",
        protein_vp3 = "data/sequences/alignments/filtered/VP3_filtered.fasta",
        protein_vp1 = "data/sequences/alignments/filtered/VP1_filtered.fasta",
        protein_2a = "data/sequences/alignments/filtered/2A_filtered.fasta",
        protein_2b = "data/sequences/alignments/filtered/2B_filtered.fasta",
        protein_2c = "data/sequences/alignments/filtered/2C_filtered.fasta",
        protein_3a = "data/sequences/alignments/filtered/3A_filtered.fasta",
        protein_3b = "data/sequences/alignments/filtered/3B_filtered.fasta",
        protein_3c = "data/sequences/alignments/filtered/3C_filtered.fasta",
        protein_3d = "data/sequences/alignments/filtered/3Dpol_filtered.fasta",
        five_prime_utr = "data/sequences/alignments/filtered/5UTR_filtered.fasta",
        three_prime_utr = "data/sequences/alignments/filtered/3UTR_filtered.fasta"
    params:
        max_gap_fraction_gene = 0.1,
        max_gap_fraction_utr = 0.5   # Allow more gaps in the UTRs, as they seem more variable (e.g., variable sequencing starts and ends in the 5'UTR and 3'UTR)
    run:
        shell("python scripts/remove_gappy_seqs_from_alignments.py --alignments {input.alignments} --max_gap_fraction_gene {params.max_gap_fraction_gene} --max_gap_fraction_utr {params.max_gap_fraction_utr} --output_dir data/sequences/alignments/filtered")

rule calculate_pdistances:
    input:
        alignments = rules.filter_gappy_seqs.output
    output:
        pdistances_csv = "data/pdistances/pdistances.csv"
    run:
        shell("python scripts/calculate_pdistances.py --alignments {input.alignments} --output {output.pdistances_csv}")



# # Build trees
# rule augur_trees:
#     input:
#         full_genome = rules.mafft_full_genomes.output.aligned_fasta,
#         five_prime_utr = rules.filter_five_prime_seqs.output.filtered_alignment,
#         protein_vp1 = rules.split_coding_region_by_gene.output.protein_vp1,
#         protein_3d = rules.split_coding_region_by_gene.output.protein_3d
#     output:
#         full_genome_tree = "data/trees/full_genome_tree.nwk",
#         five_prime_utr_tree = "data/trees/five_prime_utr_tree.nwk",
#         protein_vp1_tree = "data/trees/vp1_tree.nwk",
#         protein_3d_tree = "data/trees/3d_tree.nwk"
#     run:
#         shell("augur tree --alignment {input.full_genome} --output {output.full_genome_tree}")
#         shell("augur tree --alignment {input.five_prime_utr} --output {output.five_prime_utr_tree}")
#         shell("augur tree --alignment {input.protein_vp1} --output {output.protein_vp1_tree}")
#         shell("augur tree --alignment {input.protein_3d} --output {output.protein_3d_tree}")


# rule augur_translate:
#     input:
#         protein_vp1_tree = rules.augur_trees.output.protein_vp1_tree,
#         protein_3d_tree = rules.augur_trees.output.protein_3d_tree
#     output:
#         aa_muts_vp1 = "data/trees/aa_muts_vp1.json",
#         aa_muts_3d = "data/trees/aa_muts_3d.json"
#     run:
#         shell("augur translate --tree {input.protein_vp1_tree} --output-node-data {output.aa_muts_vp1}")
#         shell("augur translate --tree {input.protein_3d_tree} --output-node-data {output.aa_muts_3d}")









# rule filter_full_npevc_genomes:
#     input:
#         fasta = rules.split_data.output.npevc_sequences
#     output:
#         filtered_fasta = "data/sequences/npevc_full_genomes.fasta"
#     params:
#         min_length = 6000
#     script:
#         "scripts/filter_fasta_by_length.py"

# rule concat_npevc_sabin123:
#     input:
#         sabin123 = files.sabin123_sequences,
#         npevc = rules.filter_full_npevc_genomes.output.filtered_fasta
#     output:
#         concat_fasta = "data/sequences/concat_npevc_sabin123.fasta"
#     run:
#         shell("cat {input.sabin123} {input.npevc} > {output.concat_fasta}")

# rule slice_npevc_sabin123:
#     input:
#         fasta = rules.concat_npevc_sabin123.output.concat_fasta
#     output:
#         sliced_fasta = "data/sequences/sliced_concat_npevc_sabin123.fasta"
#     params:
#         n_seqs = 20
#     script:
#         "scripts/slice_fasta.py"

# rule mafft_npevc_sabin123:
#     input:
#         fasta = rules.slice_npevc_sabin123.output.sliced_fasta
#     output:
#         aligned_fasta = "data/sequences/alignments/aligned_npevc_sabin123.fasta"
#     run:
#         shell("mafft --auto {input.fasta} > {output.aligned_fasta}")

# rule replace_nonstandard_chars:
#     input:
#         fasta = rules.mafft_npevc_sabin123.output.aligned_fasta
#     output:
#         clean_fasta = "data/sequences/alignments/aligned_npevc_sabin123_clean.fasta"
#     run:
#         shell("sed -e '/^[^>]/s/[^ATGCNatgcn-]/n/g' {input.fasta} > {output.clean_fasta}")
