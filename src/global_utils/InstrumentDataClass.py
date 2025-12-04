from dataclasses import dataclass
import os
import yaml

@dataclass
class InstrumentConfig:
    name: str
    image_extension: int
    id_slice: tuple
    readout_speed: float
    gain: float
    overscan: dict
    prefix: str

def load_instrument_config(name, filename="instruments.yaml"):
    here = os.path.dirname(__file__)
    fullpath = os.path.join(here, filename)

    with open(fullpath, "r") as f:
        data = yaml.safe_load(f)

    if name not in data:
        raise ValueError(f"Instrument '{name}' not found in {filename}")

    cfg = data[name]

    return InstrumentConfig(
        name=name,
        image_extension=cfg["image_extension"],
        id_slice=tuple(cfg["id_slice"]),
        readout_speed=cfg.get("readout_speed"),
        gain=cfg.get("gain"),
        overscan=cfg.get("overscan"),
        prefix=cfg["prefix"]
    )
