from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Literal, Optional, Union, Any
import yaml
from pydantic import BaseModel, Field


class ImageExtensionRule(BaseModel):
    rule: Literal["explicit", "nwin_map"]
    # explicit: list of extensions (e.g., [0] or ["SCI","ERR"])
    list: Optional[List[Union[int, str]]] = None
    # nwin_map: map from NWIN -> list of extensions
    map: Optional[Dict[int, List[int]]] = None


class BufferExtractRuleStatic(BaseModel):
    rule: Literal["static"]
    left: int
    right: int


class BufferExtractRuleEdgeThreshold(BaseModel):
    rule: Literal["edge_threshold"]
    left_threshold: float
    right_threshold: float
    buffer_pixels: int
    default_buffer_pixels: int = 0


BufferExtractRule = Union[BufferExtractRuleStatic, BufferExtractRuleEdgeThreshold]


class InstrumentConfig(BaseModel):
    instrument_id: str

    observatory: Literal["ground", "space"]
    image_extension: ImageExtensionRule

    id_slice: Optional[List[int]] = None

    prefix: str
    buffer_pixels_find_trace: int = 0
    buffer_pixels_extract_trace: Optional[BufferExtractRule] = None

    diameter: Optional[float] = None
    altitude: Optional[float] = None

    gain: Optional[Union[float, Dict[str, float]]] = None
    readnoise: Optional[Union[float, Dict[str, float]]] = None

    dark_current: Optional[float] = None
    W3: Optional[float] = None
    h0: Optional[float] = None
    scintillation: Optional[float] = None

    def gain_for(self, speed: Optional[str]) -> Optional[float]:
        if self.gain is None:
            return None
        if isinstance(self.gain, (int, float)):
            return float(self.gain)
        if speed is None:
            return None
        return self.gain.get(speed)

    def readnoise_for(self, speed: Optional[str]) -> Optional[float]:
        if self.readnoise is None:
            return None
        if isinstance(self.readnoise, (int, float)):
            return float(self.readnoise)
        if speed is None:
            return None
        return self.readnoise.get(speed)


def load_instrument_config(instrument_id: str, filename: str = "instruments.yaml") -> InstrumentConfig:
    path = Path(__file__).with_name(filename)
    data = yaml.safe_load(path.read_text(encoding="utf-8"))

    if instrument_id not in data:
        raise ValueError(f"Instrument '{instrument_id}' not found in {path}")

    cfg = data[instrument_id]

    return InstrumentConfig(instrument_id=instrument_id, **cfg)