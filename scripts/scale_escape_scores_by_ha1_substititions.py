import argparse
from augur.utils import read_node_data, write_json
import numpy as np
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Scale escape scores by number of HA1 substitutions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--distances", required=True, help="node data JSON file of distances per strain with 'welsh_escape' and 'ha1' values")
    parser.add_argument("--output", required=True, help="node data JSON per strain with scaled escape scores")
    args = parser.parse_args()

    # Load distances.
    distances = read_node_data(args.distances)["nodes"]

    scaled_escape_scores = {}
    for node, node_values in distances.items():
        if "welsh_escape" in node_values and "ha1" in node_values and node_values["ha1"] > 0:
            scaled_escape_scores[node] = {
                "welsh_escape_per_ha1": node_values["welsh_escape"] / node_values["ha1"],
            }

    # Export distances to JSON.
    write_json({"nodes": scaled_escape_scores}, args.output)
