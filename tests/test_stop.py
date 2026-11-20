import math
from pathlib import Path

import pytest

from argextender.argextender import Config, flank_growth, parse_args, stop_reasons


def cfg(**kw: object) -> Config:
    base = {
        "r1": Path("a"),
        "r2": Path("b"),
        "singletons": None,
        "mapstat": Path("m"),
        "panres": Path("p"),
        "strategy": "A-scaffold",
        "out": Path("o"),
    }
    return Config(**{**base, **kw})  # type: ignore[arg-type]


def test_flank_growth() -> None:
    assert flank_growth(100, None, 1) is None  # the first assembly has no previous value
    assert flank_growth(110, 100, 1) == pytest.approx(0.10)
    assert flank_growth(133.1, 100, 3) == pytest.approx(0.10)  # B: per-round equivalent of 3 rounds
    assert flank_growth(90, 100, 1) == pytest.approx(-0.10)  # the flanks got shorter
    assert flank_growth(5, 0, 1) == math.inf
    assert flank_growth(0, 0, 1) == 0.0


def test_no_growth_rule() -> None:
    c = cfg(min_growth=2.0)
    assert stop_reasons(c, 2, 0.019, 0, False) == ["no_growth"]
    assert stop_reasons(c, 2, 0.021, 0, False) == []
    assert stop_reasons(c, 1, None, 0, False) == []
    assert stop_reasons(cfg(min_growth=0), 2, -0.5, 0, False) == []  # 0 disables it


def test_max_rounds_max_flank_and_empty_target() -> None:
    assert stop_reasons(cfg(rounds=5), 5, 0.5, 0, False) == ["max_rounds"]
    assert stop_reasons(cfg(max_flank=1000), 3, 0.5, 0.5, False) == ["max_flank"]  # 50 % = the median
    assert stop_reasons(cfg(max_flank=1000), 3, 0.5, 0.49, False) == []
    assert stop_reasons(cfg(max_flank=1000, max_flank_share=100), 3, 0.5, 0.99, False) == []
    assert stop_reasons(cfg(max_flank=0), 3, 0.5, 1.0, False) == []  # 0 disables it
    assert stop_reasons(cfg(rounds=3, min_growth=2), 3, 0.0, 0, True) == [
        "max_rounds",
        "no_growth",
        "empty_target",
    ]


def test_parse_args_assemble_every() -> None:
    base = ["--r1", "a", "--r2", "b", "--mapstat", "m", "--panres", "p", "--out", "o"]
    assert parse_args(base).strategy == "A-scaffold"  # the default
    assert parse_args(base).assemble_every == 1
    with pytest.raises(SystemExit):
        parse_args([*base, "--strategy", "A-graph"])  # not offered any more
    assert parse_args([*base, "--strategy", "B-graph"]).assemble_every == 3
    assert parse_args([*base, "--strategy", "B-graph", "--assemble-every", "2"]).assemble_every == 2
    with pytest.raises(SystemExit):
        parse_args([*base, "--strategy", "A-scaffold", "--assemble-every", "2"])
