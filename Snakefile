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
    pass

rule merge_metadata:
    pass

rule tree:
    pass

rule root_and_prune_tree:
    pass

rule refine:
    pass

rule distances:
    pass

rule subsample_future:
    pass

rule frequencies:
    pass

rule distances_to_future:
    pass

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
