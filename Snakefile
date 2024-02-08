FUTURE_SEASON_BY_CURRENT_SEASON = {
    "2020-10-01": "2021-10-01",
    "2021-02-01": "2022-02-01",
    "2021-10-01": "2022-10-01",
    "2022-02-01": "2023-02-01",
    "2022-10-01": "2023-10-01",
    "2023-02-01": "2024-02-01",
}

rule subsample:
    pass

rule get_nextclade_dataset:
    pass

rule align:
    input:
        sequences="results/{season}/subsampled_sequences.fasta",
    output:
        alignment="results/{season}/aligned.fasta",
        nextclade_annotations="results/{season}/nextclade.tsv",
    conda: "env.yaml"
    threads: 8
    shell:
        """
        nextclade3 run \
            -j {threads} \
            -D {input.nextclade_directory} \
            --output-fasta {output.alignment} \
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
            --fields "seqName;strain" \
            {input.metadata} \
            {input.nextclade_annotations} > {output.metadata}
        """

rule tree:
    input:
        alignment="results/{season}/aligned.fasta",
    output:
        tree="results/{season}/tree_raw.nwk",
    params:
        tree_builder_args="-ninit 10 -n 4 -czb -nt AUT",
    conda: "env.yaml"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args="{params.tree_builder_args}" \
            --output {output.tree}
        """

rule root_and_prune_tree:
    input:
        tree="results/{season}/tree_raw.nwk",
    output:
        tree="results/{season}/tree_rooted.nwk",
    params:
        root="A/Wisconsin/67/2005-egg",
    conda: "env.yaml"
    shell:
        """
        python3 scripts/root_and_prune_tree.py \
            --tree {input.tree} \
            --root {params.root} \
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

rule distances:
    input:
        tree="results/{season}/tree.nwk",
        alignments = [
            "results/{season}/translations/SigPep.fasta",
            "results/{season}/translations/HA1.fasta",
            "results/{season}/translations/HA2.fasta",
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
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --output {output.distances}
        """

rule frequencies:
    input:
        tree="results/{season}/tree.nwk",
        metadata="results/{season}/metadata.tsv",
    output:
        frequencies="results/{season}/frequencies.json",
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

rule convert_node_data_to_table:
    input:
        tree="results/{season}/tree.nwk",
        metadata="results/{season}/metadata.tsv",
        node_data="results/{season}/aa_sequences.json",
    output:
        table = "builds/{build_name}/{segment}/node-data-table.tsv",
    conda: "env.yaml"
    shell:
        """
        python3 scripts/node_data_to_table.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --jsons {input.node_data} \
            --annotations season="{wildcards.season}" \
            --output {output}
        """

rule convert_frequencies_to_table:
    input:
        tree="results/{season}/tree.nwk",
        frequencies="results/{season}/frequencies.json",
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
        python3 scripts/calculate_weighted_distance.py \
            --current-tip-attributes {input.current_season_tips} \
            --future-tip-attributes {input.future_season_tips} \
            --distance-attribute-name {params.weighted_distance_to_observed_future} \
            --output {output.weighted_distances_to_future}
        """

rule export:
    input:
        tree="results/{season}/tree.nwk",
        metadata="results/{season}/metadata.tsv",
        branch_lengths="results/{season}/branch_lengths.json",
        epitope_distances="results/{season}/epitope_distances.json",
        weighted_distances_to_future="results/{season}/weighted_distances_to_future.json",
    output:
        auspice_json="auspice/{season}.json",
        root_sequence="auspice/{season}_root-sequence.json",
    conda: "env.yaml",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.epitope_distances} {input.weighted_distances_to_future} \
            --include-root-sequence \
            --minify-json \
            --output {output.auspice_json}
        """

rule json_to_table:
    input:
        auspice_json="auspice/{season}.json",
    output:
        distances="results/{season}/distances_by_strain.tsv",
    conda: "env.yaml"
    shell:
        """
        python3 scripts/auspice_tree_to_table.py \
            --tree {input.auspice_json} \
            --attributes welsh_ep welsh_escape weighted_distance_to_future \
            --output-metadata {output.distances}
        """

rule aggregate_distances:
    input:
        distances=expand("results/{season}/distances_by_strain.tsv", season=SEASONS),
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
