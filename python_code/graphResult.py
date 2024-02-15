from dataclasses import dataclass


@dataclass
class GraphResult:
    """
    Output the causal graph
    """
    causal_links: list
