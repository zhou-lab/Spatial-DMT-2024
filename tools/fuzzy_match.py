"""Exact-preserving fast paths for fuzzysearch Levenshtein matching."""

import functools

from fuzzysearch import (
    LevenshteinSearchParams,
    choose_search_class,
    find_near_matches as fuzzy_find_near_matches,
)
from fuzzysearch.common import Match, group_matches


def _consolidate_matches_deterministically(matches):
    """Match fuzzysearch consolidation with stable positional tie-breaking."""
    groups = group_matches(matches)
    return sorted(
        min(
            group,
            key=lambda match: (
                match.dist,
                -(match.end - match.start),
                match.start,
                match.end,
            ),
        )
        for group in groups
    )


def find_near_matches_deterministic(subsequence, sequence, max_l_dist):
    params = LevenshteinSearchParams(max_l_dist=max_l_dist)
    search_class = choose_search_class(params)
    matches = search_class.search(subsequence, sequence, params)
    return _consolidate_matches_deterministically(matches)


@functools.lru_cache(maxsize=None)
def _exact_seed_parts(subsequence, max_l_dist):
    """Return k+1 disjoint seeds; any Levenshtein-k hit retains one exactly."""
    parts = max_l_dist + 1
    base_size, remainder = divmod(len(subsequence), parts)
    seeds = []
    start = 0
    for i in range(parts):
        width = base_size + (1 if i < remainder else 0)
        seeds.append(subsequence[start:start + width])
        start += width
    return tuple(seeds)


def find_near_matches_seeded(subsequence, sequence, max_l_dist):
    seeds = _exact_seed_parts(subsequence, max_l_dist)
    if not any(seed in sequence for seed in seeds):
        return []
    return fuzzy_find_near_matches(
        subsequence, sequence, max_l_dist=max_l_dist
    )


def find_near_matches_exact_seeded(subsequence, sequence, max_l_dist):
    """Return exact hits only when seeds prove no separate fuzzy group exists."""
    exact = []
    start = sequence.find(subsequence)
    while start >= 0:
        end = start + len(subsequence)
        exact.append(Match(start, end, 0, sequence[start:end]))
        start = sequence.find(subsequence, start + 1)
    if not exact:
        return find_near_matches_seeded(subsequence, sequence, max_l_dist)

    if len(exact) != 1:
        return fuzzy_find_near_matches(
            subsequence, sequence, max_l_dist=max_l_dist
        )

    seeds = _exact_seed_parts(subsequence, max_l_dist)
    for seed in seeds:
        seed_start = sequence.find(seed)
        while seed_start >= 0:
            seed_end = seed_start + len(seed)
            if not any(
                seed_start < match.end and seed_end > match.start
                for match in exact
            ):
                return fuzzy_find_near_matches(
                    subsequence, sequence, max_l_dist=max_l_dist
                )
            seed_start = sequence.find(seed, seed_start + 1)

    return exact


MATCHERS = {
    "fuzzy": fuzzy_find_near_matches,
    "deterministic": find_near_matches_deterministic,
    "seeded": find_near_matches_seeded,
    "exact-seeded": find_near_matches_exact_seeded,
}
