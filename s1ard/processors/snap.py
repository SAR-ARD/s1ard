from __future__ import annotations
from typing import Any
from .base import Processor
from configparser import ConfigParser

from s1ard import snap


class SnapProcessor(Processor):
    """Adapter that exposes s1ard.snap functions as instance methods."""
    name = "snap"
    
    @staticmethod
    def config_to_string(config: dict[str, Any]) -> dict[str, str]:
        return snap.config_to_string(config)
    
    @staticmethod
    def find_datasets(scene: str, outdir: str, epsg: int) -> dict[str, str] | None:
        return snap.find_datasets(scene, outdir, epsg)
    
    @staticmethod
    def get_config_keys() -> list[str]:
        return snap.get_config_keys()
    
    @staticmethod
    def get_config_section(parser: ConfigParser, **kwargs: dict[str, str]) -> dict[str, str]:
        return snap.get_config_section(parser, **kwargs)
    
    @staticmethod
    def get_metadata(scene: str, outdir: str) -> dict[str, Any]:
        return snap.get_metadata(scene, outdir)
    
    @staticmethod
    def lsm_encoding() -> dict[str, int]:
        return snap.lsm_encoding()
    
    @staticmethod
    def process(scene: str, outdir: str, tmpdir: str, dem: str, measurement: str,
                spacing: int | float, dem_resampling_method: str,
                img_resampling_method: str, rlks: int | None, azlks: int | None,
                export_extra: list[str] | None, allow_res_osv: bool, clean_edges: bool,
                clean_edges_pixels: int, neighbors: list[str] | None,
                gpt_args: list[str] | None, cleanup: bool) -> None:
        return snap.process(scene=scene, outdir=outdir, tmpdir=tmpdir, dem=dem,
                            measurement=measurement, spacing=spacing,
                            dem_resampling_method=dem_resampling_method,
                            img_resampling_method=img_resampling_method,
                            rlks=rlks, azlks=azlks, export_extra=export_extra,
                            allow_res_osv=allow_res_osv, clean_edges=clean_edges,
                            clean_edges_pixels=clean_edges_pixels,
                            neighbors=neighbors, gpt_args=gpt_args, cleanup=cleanup)
    
    @staticmethod
    def translate_annotation(annotation: list[str], measurement: str) -> list[str]:
        return snap.translate_annotation(annotation, measurement)
    
    @staticmethod
    def version_dict() -> dict[str, str]:
        return snap.version_dict()
