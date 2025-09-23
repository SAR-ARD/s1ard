from __future__ import annotations
from importlib.metadata import entry_points
from typing import Type
from .base import Processor

_GROUP = "s1ard.processors"


def available_processors() -> dict[str, str]:
    eps = entry_points().select(group=_GROUP)
    return {ep.name: f"{ep.module}:{ep.attr}" for ep in eps}


def load_processor(name: str) -> Type[Processor]:
    eps = {ep.name: ep for ep in entry_points().select(group=_GROUP)}
    if name not in eps:
        installed = ", ".join(sorted(eps)) or "(none)"
        raise ModuleNotFoundError(f"Processor '{name}' not found. Installed: {installed}")
    return eps[name].load()
