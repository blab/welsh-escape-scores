# Analysis of predictive value of escape scores from Welsh et al. 2023 experimental data

This workflow analyzes the historical value of escape scores per HA site and amino acid from [Welsh et al. 2023](https://doi.org/10.1101/2023.12.12.571235)'s antigenic escape assays.
We produced one phylogenetic tree and measurements panel per past vaccine composition meeting season starting with October 2020 and ending with February 2023.
For each season, we calculated escape scores per H3N2 HA sequence based on the amino acid mutations present in each sequence and we calculated the weighted amino acid distance of each sequence to the future population in 12 months.
To understand how predictive the escape scores would have been historically, we scaled each strain's escape scores by the number of HA1 substitutions per strain and calculated the correlation of the scaled escape scores with the distance to the future population.
When escape scores were predictive of the future, we expected higher scores to correspond to lower distances to the future.
For each tree, we assigned both the historical clade labels that were used at the time and modern "subclade" labels that reflect finer resolution of genetic diversity.

## Run the analysis

[Install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (version 7.32.4 or later) and run the following command from the top-level of this repository.

``` bash
snakemake
```

## Explore phylogenetic trees and measurements

This workflow produces Nextstrain phylogenetic trees and measurements panels (see [Lee et al. 2023](https://doi.org/10.3389/fbinf.2023.1069487)) for historical H3N2 HA sequences and the corresponding escape scores per serum for those sequences.
Each link below corresponds to a Nextstrain view for one season showing the scatterplot view of distance to the future by scaled escape score in the left panel and the measurements from escape scores per serum sample and historical clade label on the right panel.

 - [October 2020](https://nextstrain.org/groups/blab/welsh-escape-scores/2020-10-01?branches=hide&dmin=2020-05-31&l=scatter&p=grid&regression=show&scatterX=welsh_escape_nonnegative_scores_all_sera_per_ha1&scatterY=weighted_distance_to_observed_future)
 - [February 2021](https://nextstrain.org/groups/blab/welsh-escape-scores/2021-02-01?branches=hide&dmin=2020-09-07&l=scatter&p=grid&regression=show&scatterX=welsh_escape_nonnegative_scores_all_sera_per_ha1&scatterY=weighted_distance_to_observed_future)
 - [October 2021](https://nextstrain.org/groups/blab/welsh-escape-scores/2021-10-01?branches=hide&dmin=2021-03-19&l=scatter&p=grid&regression=show&scatterX=welsh_escape_nonnegative_scores_all_sera_per_ha1&scatterY=weighted_distance_to_observed_future)
 - [February 2022](https://nextstrain.org/groups/blab/welsh-escape-scores/2022-02-01?branches=hide&dmin=2021-07-07&l=scatter&p=grid&regression=show&scatterX=welsh_escape_nonnegative_scores_all_sera_per_ha1&scatterY=weighted_distance_to_observed_future)
 - [October 2022](https://nextstrain.org/groups/blab/welsh-escape-scores/2022-10-01?branches=hide&dmin=2022-05-09&l=scatter&p=grid&regression=show&scatterX=welsh_escape_nonnegative_scores_all_sera_per_ha1&scatterY=weighted_distance_to_observed_future)
 - [February 2023](https://nextstrain.org/groups/blab/welsh-escape-scores/2023-02-01?branches=hide&dmin=2022-09-12&l=scatter&p=grid&regression=show&scatterX=welsh_escape_nonnegative_scores_all_sera_per_ha1&scatterY=weighted_distance_to_observed_future)
