import datetime

FUTURE_SEASON_BY_CURRENT_SEASON = {
    "2020-10-01": "2021-10-01",
    "2021-02-01": "2022-02-01",
    "2021-10-01": "2022-10-01",
    "2022-02-01": "2023-02-01",
    "2022-10-01": "2023-10-01",
    "2023-02-01": "2024-02-01",
}

rule all:
    input:
        "results/distance_to_the_future_by_epitope_score_and_season.pdf",

rule download_metadata:
    output:
        metadata="data/metadata.tsv.xz",
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/h3n2/metadata.tsv.xz",
    conda: "env.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} {output.metadata}
        """

rule download_sequences:
    output:
        sequences="data/sequences.fasta.xz",
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/h3n2/ha/sequences.fasta.xz",
    conda: "env.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} {output.sequences}
        """

rule index_sequence:
    input:
        sequences="data/sequences.fasta.xz",
    output:
        sequence_index="data/sequence_index.tsv",
    conda: "env.yaml"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

def get_min_date_for_season(wildcards):
    """Set the minimum date to 1 year (52 weeks) before the date for the current
    season.

    """
    current_date = datetime.datetime.strptime(wildcards.season, "%Y-%m-%d")
    min_date = current_date - datetime.timedelta(weeks=52)
    return min_date.strftime("%Y-%m-%d")

rule subsample:
    input:
        metadata="data/metadata.tsv.xz",
        sequences="data/sequences.fasta.xz",
        sequence_index="data/sequence_index.tsv",
        outliers="config/outliers.txt",
    output:
        metadata="results/{season}/subsampled_metadata.tsv",
        sequences="results/{season}/subsampled_sequences.fasta",
    params:
        min_length=1700,
        exclude_ambiguous_dates_by="any",
        min_date=get_min_date_for_season,
    conda: "env.yaml"
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --exclude {input.outliers} \
            --min-length {params.min_length} \
            --non-nucleotide \
            --exclude-ambiguous-dates-by {params.exclude_ambiguous_dates_by} \
            --min-date {params.min_date:q} \
            --max-date {wildcards.season} \
            --group-by region year month \
            --subsample-max-sequences 100 \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences}
        """

rule get_nextclade_dataset:
    output:
        nextclade_directory=directory("nextclade"),
    params:
        dataset="nextstrain/flu/h3n2/ha/wisconsin-67-2005",
    conda: "env.yaml"
    shell:
        """
        nextclade dataset get \
            -n {params.dataset:q} \
            --output-dir {output.nextclade_directory}
        """

rule align:
    input:
        sequences="results/{season}/subsampled_sequences.fasta",
        nextclade_directory="nextclade"
    output:
        alignment="results/{season}/aligned.fasta",
        nextclade_annotations="results/{season}/nextclade.tsv",
        translations=[
            "results/{season}/translations/SigPep.fasta",
            "results/{season}/translations/HA1.fasta",
            "results/{season}/translations/HA2.fasta",
        ],
    params:
        cds_selection="SigPep,HA1,HA2",
        translations_template=lambda wildcards: f"results/{wildcards.season}/translations/{{cds}}.fasta",
    conda: "env.yaml"
    threads: 8
    shell:
        """
        nextclade run \
            -j {threads} \
            -D {input.nextclade_directory} \
            --cds-selection {params.cds_selection} \
            --include-reference \
            --output-fasta {output.alignment} \
            --output-translations {params.translations_template} \
            --output-tsv {output.nextclade_annotations} \
            {input.sequences}
        """

rule merge_metadata:
    input:
        metadata="results/{season}/subsampled_metadata.tsv",
        nextclade_annotations="results/{season}/nextclade.tsv",
    output:
        metadata="results/{season}/metadata.tsv",
    conda: "env.yaml"
    shell:
        """
        csvtk join \
            -t \
            --fields "strain;seqName" \
            {input.metadata} \
            {input.nextclade_annotations} > {output.metadata}
        """

rule tree:
    input:
        alignment="results/{season}/aligned.fasta",
    output:
        tree="results/{season}/tree_raw.nwk",
    params:
        tree_builder_args="-ninit 10 -n 4 -czb -nt AUTO",
    conda: "env.yaml"
    threads: 8
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args="{params.tree_builder_args}" \
            --override-default-args \
            --nthreads {threads} \
            --output {output.tree}
        """

rule root_and_prune_tree:
    input:
        tree="results/{season}/tree_raw.nwk",
    output:
        tree="results/{season}/tree_rooted.nwk",
    params:
        root="CY163680.1",
    conda: "env.yaml"
    shell:
        """
        python3 scripts/root_and_prune_tree.py \
            --tree {input.tree} \
            --root '{params.root}' \
            --output {output.tree}
        """

rule refine:
    input:
        metadata="results/{season}/metadata.tsv",
        alignment="results/{season}/aligned.fasta",
        tree="results/{season}/tree_rooted.nwk",
    output:
        tree="results/{season}/tree.nwk",
        node_data="results/{season}/branch_lengths.json",
    conda: "env.yaml"
    shell:
        """
        augur refine \
            --metadata {input.metadata} \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --keep-root \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree="results/{season}/tree.nwk",
        alignment="results/{season}/aligned.fasta",
        translations = [
            "results/{season}/translations/SigPep.fasta",
            "results/{season}/translations/HA1.fasta",
            "results/{season}/translations/HA2.fasta",
        ],
        annotation = "nextclade/genome_annotation.gff3",
    output:
        node_data = "results/{season}/muts.json",
        translations = [
            "results/{season}/translations/SigPep_withInternalNodes.fasta",
            "results/{season}/translations/HA1_withInternalNodes.fasta",
            "results/{season}/translations/HA2_withInternalNodes.fasta",
        ],
    params:
        inference = "joint",
        genes = ["SigPep", "HA1", "HA2"],
        input_translations = lambda wildcards: f"results/{wildcards.season}/translations/%GENE.fasta",
        output_translations = lambda wildcards: f"results/{wildcards.season}/translations/%GENE_withInternalNodes.fasta",
    conda: "env.yaml"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --annotation {input.annotation} \
            --genes {params.genes} \
            --translations "{params.input_translations}" \
            --output-node-data {output.node_data} \
            --output-translations "{params.output_translations}" \
            --inference {params.inference}
        """

rule distances:
    input:
        tree="results/{season}/tree.nwk",
        alignments = [
            "results/{season}/translations/SigPep_withInternalNodes.fasta",
            "results/{season}/translations/HA1_withInternalNodes.fasta",
            "results/{season}/translations/HA2_withInternalNodes.fasta",
        ],
        distance_maps = [
            "config/welsh_epitope_sites.json",
            "config/welsh_escape_by_site_and_amino_acid.json",
        ],
    output:
        distances="results/{season}/epitope_distances.json",
    params:
        genes = ["SigPep", "HA1", "HA2"],
        comparisons = ["root", "root"],
        attribute_names = ["welsh_ep", "welsh_escape"],
    conda: "env.yaml"
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes:q} \
            --compare-to {params.comparisons:q} \
            --attribute-name {params.attribute_names:q} \
            --map {input.distance_maps} \
            --output {output.distances}
        """

rule frequencies:
    input:
        tree="results/{season}/tree.nwk",
        metadata="results/{season}/metadata.tsv",
    output:
        frequencies="auspice/{season}_tip-frequencies.json",
    params:
        narrow_bandwidth = 1 / 12.0,
        wide_bandwidth = 3 / 12.0,
        proportion_wide = 0.0,
        pivot_interval = 1
    conda: "env.yaml"
    shell:
        """
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --wide-bandwidth {params.wide_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --pivot-interval {params.pivot_interval} \
            --max-date {wildcards.season} \
            --output {output}
        """

rule concatenate_aa_sequences:
    input:
        translations=[
            "results/{season}/translations/SigPep.fasta",
            "results/{season}/translations/HA1.fasta",
            "results/{season}/translations/HA2.fasta",
        ],
    output:
        aa_sequences_table="results/{season}/node_data.tsv",
    conda: "env.yaml"
    shell:
        """
        seqkit concat --quiet {input.translations} \
            | seqkit fx2tab -H -Q \
            | csvtk rename -t -f 1,2 -n strain,aa_sequence > {output.aa_sequences_table}
        """

rule convert_frequencies_to_table:
    input:
        tree="results/{season}/tree.nwk",
        frequencies="auspice/{season}_tip-frequencies.json",
    output:
        frequencies = "results/{season}/frequencies.tsv",
    conda: "env.yaml"
    shell:
        """
        python3 scripts/frequencies_to_table.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --output {output.frequencies}
        """

rule merge_node_data_and_frequencies:
    input:
        node_data="results/{season}/node_data.tsv",
        frequencies="results/{season}/frequencies.tsv",
    output:
        tip_attributes="results/{season}/tip_attributes.tsv",
    conda: "env.yaml"
    shell:
        """
        csvtk join \
            -t \
            -f strain \
            {input.node_data} \
            {input.frequencies} > {output.tip_attributes}
        """

def get_future_season_tips_by_current_season(wildcards):
    future_season = FUTURE_SEASON_BY_CURRENT_SEASON[wildcards.season]
    return f"results/{future_season}/tip_attributes.tsv"

rule distances_to_future:
    input:
        current_season_tips="results/{season}/tip_attributes.tsv",
        future_season_tips=get_future_season_tips_by_current_season,
    output:
        weighted_distances_to_future="results/{season}/weighted_distances_to_future.json",
    params:
        distance_attribute_name="weighted_distance_to_observed_future",
    conda: "env.yaml"
    shell:
        """
        python3 scripts/calculate_weighted_distances.py \
            --current-tip-attributes {input.current_season_tips} \
            --future-tip-attributes {input.future_season_tips} \
            --distance-attribute-name {params.distance_attribute_name:q} \
            --output {output.weighted_distances_to_future}
        """

rule export:
    input:
        tree="results/{season}/tree.nwk",
        metadata="results/{season}/metadata.tsv",
        branch_lengths="results/{season}/branch_lengths.json",
        muts="results/{season}/muts.json",
        epitope_distances="results/{season}/epitope_distances.json",
        weighted_distances_to_future="results/{season}/weighted_distances_to_future.json",
        auspice_config="config/auspice_config.json",
    output:
        auspice_json="auspice/{season}.json",
    conda: "env.yaml",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.muts} {input.epitope_distances} {input.weighted_distances_to_future} \
            --auspice-config {input.auspice_config} \
            --minify-json \
            --output {output.auspice_json}
        """

rule json_to_table:
    input:
        auspice_json="auspice/{season}.json",
    output:
        distances="results/{season}/distances_by_strain.tsv",
    params:
        attributes=[
            "welsh_ep",
            "welsh_escape",
            "weighted_distance_to_observed_future",
        ]
    conda: "env.yaml"
    shell:
        """
        python3 scripts/auspice_tree_to_table.py \
            --tree {input.auspice_json} \
            --attributes {params.attributes:q} \
            --output-metadata {output.distances}
        """

rule aggregate_distances:
    input:
        distances=expand("results/{season}/distances_by_strain.tsv", season=list(FUTURE_SEASON_BY_CURRENT_SEASON.keys())),
    output:
        distances="results/distances.tsv",
    conda: "env.yaml"
    shell:
        """
        tsv-append -H {input.distances} > {output.distances}
        """

rule plot_distances:
    input:
        distances="results/distances.tsv",
    output:
        distances_figure="results/distance_to_the_future_by_epitope_score_and_season.pdf",
    conda: "env.yaml"
    notebook:
        "notebooks/plot-distances.py.ipynb"
