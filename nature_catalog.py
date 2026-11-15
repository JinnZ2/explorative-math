"""nature_catalog -- shapes a run's output can be compared against.

CC0. Co-authored. Every entry is CARRIED: stated from general knowledge, not
verified in this repo. A run's `nature_match` holds an id from here or null.
An empty match is a result, not a failure.
"""

CATALOG = {
    "phyllotaxis_golden_angle": {
        "what": "divergence angle between successive leaves/florets near 137.5 deg "
                "(sunflower heads, pine cones)",
        "status": "CARRIED",
    },
    "quasicrystal_icosahedral": {
        "what": "aperiodic long-range order in metallic alloys (Al-Mn, 1984); "
                "ratios involving phi in diffraction",
        "status": "CARRIED",
    },
    "fibonacci_chain_quasicrystal": {
        "what": "1D aperiodic Fibonacci sequences used as models of quasicrystal "
                "structure and in layered materials",
        "status": "CARRIED",
    },
}


def match(entry_id):
    if entry_id is None:
        return None
    if entry_id not in CATALOG:
        raise KeyError(f"not in catalog: {entry_id}")
    return entry_id
