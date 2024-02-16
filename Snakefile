import datetime

FUTURE_SEASON_BY_CURRENT_SEASON = {
    "2020-10-01": "2021-10-01",
    "2021-02-01": "2022-02-01",
    "2021-10-01": "2022-10-01",
    "2022-02-01": "2023-02-01",
    "2022-10-01": "2023-10-01",
    "2023-02-01": "2024-02-01",
}
ALL_SEASONS = sorted(
    set(list(FUTURE_SEASON_BY_CURRENT_SEASON.keys())) |
    set(list(FUTURE_SEASON_BY_CURRENT_SEASON.values()))
)

wildcard_constraints:
    season=r"\d{4}-\d{2}-\d{2}",

rule all:
    input:
        distances_by_subclade_and_escape_score="results/distance_to_the_future_by_escape_score_subclade_and_season.pdf",
        distances_by_subclade_and_lbi="results/distance_to_the_future_by_lbi_subclade_and_season.pdf",
        distances_by_historical_clade="results/distance_to_the_future_by_escape_score_historical_clade_and_season.pdf",
        escape_scores_by_historical_clade="results/escape_scores_by_historical_clade_and_season.pdf",
        auspice_jsons=expand("auspice/welsh-escape-scores_{season}.json", season=ALL_SEASONS),
        auspice_frequencies=expand("auspice/welsh-escape-scores_{season}_tip-frequencies.json", season=ALL_SEASONS),
        auspice_measurements=expand("auspice/welsh-escape-scores_{season}_measurements.json", season=ALL_SEASONS),

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
        references="config/references.txt",
    output:
        metadata="results/{season}/subsampled_metadata.tsv",
        sequences="results/{season}/subsampled_sequences.fasta",
    params:
        min_length=1700,
        exclude_ambiguous_dates_by="any",
        min_date=get_min_date_for_season,
        subsample_max_sequences=1000,
    conda: "env.yaml"
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --include {input.references} \
            --exclude {input.outliers} \
            --min-length {params.min_length} \
            --non-nucleotide \
            --exclude-ambiguous-dates-by {params.exclude_ambiguous_dates_by} \
            --min-date {params.min_date:q} \
            --max-date {wildcards.season} \
            --query "(date_submitted != 'N/A') & (date_submitted != '?') & (date_submitted < '{wildcards.season}') & (passage_category != 'egg')" \
            --group-by region year month \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences}
        """

rule get_nextclade_dataset:
    output:
        nextclade_directory=directory("nextclade"),
        reference="nextclade/reference.fasta",
        gene_map="nextclade/genome_annotation.gff3",
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
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_rate = 0.00382,
        clock_std_dev = 0.00382 / 5,
    conda: "env.yaml"
    shell:
        """
        augur refine \
            --metadata {input.metadata} \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --keep-root \
            --timetree \
            --stochastic-resolve \
            --use-fft \
            --no-covariance \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """

rule ancestral:
    input:
        tree="results/{season}/tree.nwk",
        alignment="results/{season}/aligned.fasta",
        translations = [
            "results/{season}/translations/SigPep.fasta",
            "results/{season}/translations/HA1.fasta",
            "results/{season}/translations/HA2.fasta",
        ],
        reference="nextclade/reference.fasta",
        annotation="nextclade/genome_annotation.gff3",
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
            --root-sequence {input.reference} \
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
            "config/ha1_sites.json",
            "config/welsh_epitope_sites.json",
            "config/welsh_escape_by_site_and_amino_acid.json",
        ],
    output:
        distances="results/{season}/epitope_distances.json",
    params:
        genes = ["SigPep", "HA1", "HA2"],
        comparisons = ["root", "root", "root"],
        attribute_names = ["ha1", "welsh_ep", "welsh_escape"],
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

rule scale_escape_scores_by_ha1_substitutions:
    input:
        distances="results/{season}/epitope_distances.json",
    output:
        distances="results/{season}/scaled_escape_scores.json",
    conda: "env.yaml"
    shell:
        """
        python3 scripts/scale_escape_scores_by_ha1_substititions.py \
            --distances {input.distances} \
            --output {output.distances}
        """

rule frequencies:
    input:
        tree="results/{season}/tree.nwk",
        metadata="results/{season}/metadata.tsv",
    output:
        frequencies="auspice/welsh-escape-scores_{season}_tip-frequencies.json",
    params:
        narrow_bandwidth = 1 / 12.0,
        wide_bandwidth = 3 / 12.0,
        proportion_wide = 0.0,
        pivot_interval = 1,
        min_date = get_min_date_for_season,
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
            --min-date {params.min_date} \
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
        frequencies="auspice/welsh-escape-scores_{season}_tip-frequencies.json",
    output:
        frequencies = "results/{season}/frequencies.tsv",
    conda: "env.yaml"
    shell:
        """
        python3 scripts/frequencies_to_table.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --annotations season="{wildcards.season}" \
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

rule clades:
    input:
        tree="results/{season}/tree.nwk",
        mutations="results/{season}/muts.json",
        clades="config/clades_{season}.tsv",
    output:
        clades="results/{season}/clades.json",
    params:
        membership_name="historicalclade",
        label_name="historicalclade",
    conda: "env.yaml"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.mutations} \
            --clades {input.clades} \
            --membership-name {params.membership_name:q} \
            --label-name {params.label_name:q} \
            --output {output.clades}
        """

rule lbi:
    input:
        tree="results/{season}/tree.nwk",
        branch_lengths="results/{season}/branch_lengths.json",
    output:
        lbi="results/{season}/lbi.json",
    params:
        attribute_name="lbi",
        tau=0.5,
        window=0.5,
    conda: "env.yaml"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --attribute-names {params.attribute_name} \
            --tau {params.tau} \
            --window {params.window} \
            --output {output.lbi}
        """

def get_node_data(wildcards):
    node_data = [
        f"results/{wildcards.season}/branch_lengths.json",
        f"results/{wildcards.season}/muts.json",
        f"results/{wildcards.season}/epitope_distances.json",
        f"results/{wildcards.season}/lbi.json",
    ]

    if wildcards.season in FUTURE_SEASON_BY_CURRENT_SEASON:
        node_data.extend([
            f"results/{wildcards.season}/clades.json",
            f"results/{wildcards.season}/scaled_escape_scores.json",
            f"results/{wildcards.season}/weighted_distances_to_future.json",
        ])

    return node_data

rule export:
    input:
        tree="results/{season}/tree.nwk",
        metadata="results/{season}/metadata.tsv",
        node_data=get_node_data,
        auspice_config="config/auspice_config.json",
    output:
        auspice_json="auspice/welsh-escape-scores_{season}.json",
    conda: "env.yaml",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --minify-json \
            --output {output.auspice_json}
        """

rule json_to_table:
    input:
        auspice_json="auspice/welsh-escape-scores_{season}.json",
    output:
        distances="results/{season}/distances_by_strain.tsv",
    params:
        attributes=[
            "historicalclade",
            "subclade",
            "ha1",
            "welsh_ep",
            "welsh_escape",
            "welsh_escape_per_ha1",
            "lbi",
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

rule merge_distances_with_frequencies:
    input:
        distances="results/{season}/distances_by_strain.tsv",
        frequencies="results/{season}/frequencies.tsv",
    output:
        distances="results/{season}/distances_by_strain_with_season.tsv",
    conda: "env.yaml"
    shell:
        """
        csvtk join \
            -t \
            --fields "name;strain" \
            {input.distances} \
            {input.frequencies} > {output.distances}
        """

rule aggregate_distances:
    input:
        distances=expand("results/{season}/distances_by_strain_with_season.tsv", season=list(FUTURE_SEASON_BY_CURRENT_SEASON.keys())),
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
        color_schemes="config/color_schemes.tsv",
    output:
        distances_by_subclade_and_escape_score="results/distance_to_the_future_by_escape_score_subclade_and_season.pdf",
        distances_by_subclade_and_lbi="results/distance_to_the_future_by_lbi_subclade_and_season.pdf",
        distances_by_historical_clade="results/distance_to_the_future_by_escape_score_historical_clade_and_season.pdf",
        escape_scores_by_historical_clade="results/escape_scores_by_historical_clade_and_season.pdf",
    conda: "env.yaml"
    notebook:
        "notebooks/plot-distances.py.ipynb"

rule welsh_epitope_distances_by_serum:
    input:
        tree="results/{season}/tree.nwk",
        alignments = [
            "results/{season}/translations/SigPep_withInternalNodes.fasta",
            "results/{season}/translations/HA1_withInternalNodes.fasta",
            "results/{season}/translations/HA2_withInternalNodes.fasta",
        ],
        distance_map="config/welsh_escape_by_site_and_amino_acid_by_serum/{serum}.json",
    output:
        distances="results/{season}/welsh_epitope_distances_by_serum/{serum}.json",
    params:
        genes = ["SigPep", "HA1", "HA2"],
        comparisons="root",
        attribute_names="welsh_escape_for_serum",
    conda: "env.yaml"
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_map} \
            --output {output.distances}
        """

rule convert_welsh_epitope_distances_to_table:
    input:
        tree="results/{season}/tree.nwk",
        distances="results/{season}/welsh_epitope_distances_by_serum/{serum}.json",
        distance_map="config/welsh_escape_by_site_and_amino_acid_by_serum/{serum}.json",
    output:
        distances="results/{season}/welsh_epitope_distances_by_serum/{serum}.tsv",
    params:
        distance_map_attributes=["cohort_serum", "serum", "cohort"]
    conda: "env.yaml"
    shell:
        """
        python3 scripts/convert_welsh_epitope_distances_to_table.py \
            --tree {input.tree} \
            --distances {input.distances} \
            --distance-map {input.distance_map} \
            --distance-map-attributes {params.distance_map_attributes:q} \
            --output {output.distances}
        """

def get_welsh_epitope_distances_by_serum(wildcards):
    import glob
    from pathlib import Path

    distance_maps = glob.glob("config/welsh_escape_by_site_and_amino_acid_by_serum/*.json")
    serum_samples = [
        Path(distance_map).stem
        for distance_map in distance_maps
    ]

    return [
        f"results/{wildcards.season}/welsh_epitope_distances_by_serum/{serum}.tsv"
        for serum in serum_samples
    ]

rule aggregate_welsh_epitope_distances_by_serum:
    input:
        distances=get_welsh_epitope_distances_by_serum,
    output:
        distances="results/{season}/welsh_epitope_distances_by_serum.tsv",
    conda: "env.yaml"
    shell:
        """
        csvtk concat -t {input.distances} > {output.distances}
        """

rule export_welsh_measurements:
    input:
        distances="results/{season}/welsh_epitope_distances_by_serum.tsv",
    output:
        measurements="auspice/welsh-escape-scores_{season}_measurements.json",
    conda: "env.yaml"
    params:
        default_group_by="cohort_serum",
        strain_column="strain",
        value_column="welsh_escape_for_serum",
        key="welsh_escape",
        title="Welsh et al. escape scores",
        x_axis_label="escape score",
        thresholds=[0.0],
        filters=[
            "cohort_serum",
            "cohort",
        ],
        include_columns=[
            "strain",
            "cohort_serum",
            "cohort",
            "welsh_escape_for_serum",
        ],
    shell:
        """
        augur measurements export \
            --collection {input.distances} \
            --grouping-column {params.filters} \
            --group-by {params.default_group_by} \
            --include-columns {params.include_columns:q} \
            --strain-column {params.strain_column} \
            --value-column {params.value_column} \
            --key {params.key} \
            --title {params.title:q} \
            --x-axis-label {params.x_axis_label:q} \
            --thresholds {params.thresholds} \
            --filters {params.filters} \
            --show-threshold \
            --hide-overall-mean \
            --minify-json \
            --output-json {output.measurements}
        """
