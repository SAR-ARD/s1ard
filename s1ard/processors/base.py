from __future__ import annotations
from typing import Protocol, Any
from configparser import ConfigParser


class Processor(Protocol):
    """Common interface all processors should follow."""
    name: str
    
    @staticmethod
    def config_to_string(config: dict[str, Any]) -> dict[str, str]: ...
    
    @staticmethod
    def find_datasets(scene: str, outdir: str, epsg: int) -> dict[str, str] | None: ...
    
    @staticmethod
    def get_config_keys() -> list[str]: ...
    
    @staticmethod
    def get_config_section(parser: ConfigParser, **kwargs: Any) -> dict[str, str]: ...
    
    @staticmethod
    def get_metadata(scene: str, outdir: str) -> dict[str, Any]: ...
    
    @staticmethod
    def lsm_encoding() -> dict[str, int]: ...
    
    @staticmethod
    def process(scene: str, outdir: str, **kwargs: Any) -> None: ...
    
    @staticmethod
    def translate_annotation(annotation: list[str] | None, measurement: str) -> list[str]: ...
    
    @staticmethod
    def version_dict() -> dict[str, str]: ...
