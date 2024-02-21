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
    parser.add_argument("--attributes-to-scale", required=True, nargs="+", default=["welsh_escape"], help="distance attributes to scale by number of HA1 substitutions")
    parser.add_argument("--output", required=True, help="node data JSON per strain with scaled escape scores")
    args = parser.parse_args()

    # Load distances.
    distances = read_node_data(args.distances)["nodes"]

    scaled_escape_scores = {}
    for node, node_values in distances.items():
        if all(attribute in node_values for attribute in args.attributes_to_scale) and "ha1" in node_values and node_values["ha1"] > 0:
            scaled_escape_scores[node] = {}
            for attribute in args.attributes_to_scale:
                scaled_escape_scores[node][f"{attribute}_per_ha1"] = node_values[attribute] / node_values["ha1"]

    # Export distances to JSON.
    write_json({"nodes": scaled_escape_scores}, args.output)
