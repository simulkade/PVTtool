"""Component database loader for PVTtool.

Loads the pure-component database from puredata.json at package import time.
The JSON file is exported from puredata.mat using convert_puredata.py.
"""

import json
from pathlib import Path

_DATA_PATH = Path(__file__).parent.parent / "data" / "puredata.json"
_db_cache: list[dict] | None = None


def _load_database() -> list[dict]:
    """Load the component database from puredata.json (cached)."""
    global _db_cache
    if _db_cache is None:
        with open(_DATA_PATH) as f:
            _db_cache = json.load(f)
    return _db_cache


def get_component_data(name_or_formula: str) -> dict | None:
    """Look up a component in the database by name or formula (case-insensitive).

    Args:
        name_or_formula: Component name (e.g. 'Methane') or formula (e.g. 'CH4').

    Returns:
        Component data dict, or None if not found.
    """
    db = _load_database()
    key = name_or_formula.casefold()
    for comp in db:
        if comp["name"].casefold() == key or comp["formula"].casefold() == key:
            return comp
    return None


def get_components_data(names_or_formulas: list[str]) -> tuple[list[dict], list[int]]:
    """Look up multiple components by name or formula (case-insensitive).

    Args:
        names_or_formulas: List of component names or formulas.

    Returns:
        Tuple of (found_components, not_found_indices).
        found_components: list of component data dicts (in order).
        not_found_indices: indices into the input list that were not found.
    """
    db = _load_database()
    # Build lookup dicts
    by_name: dict[str, dict] = {}
    by_formula: dict[str, dict] = {}
    for comp in db:
        by_name[comp["name"].casefold()] = comp
        by_formula[comp["formula"].casefold()] = comp

    found: list[dict] = []
    not_found: list[int] = []

    for i, query in enumerate(names_or_formulas):
        key = query.casefold()
        comp = by_formula.get(key) or by_name.get(key)
        if comp is not None:
            found.append(comp)
        else:
            not_found.append(i)

    return found, not_found
