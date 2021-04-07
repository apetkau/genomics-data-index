from typing import Dict, Set, Any


def get_values_from_signatures(sigs) -> Dict[str, Set[Any]]:
    ksizes = set()
    scales = set()

    for sig in sigs:
        ksizes.add(sig.minhash.ksize)
        scales.add(sig.minhash.scaled)

    return {
        'ksize': ksizes,
        'scaled': scales,
    }
