# config.py
from pathlib import Path

class Config:
    ATOMSK_PATH = "/home/ulumuddin/Installations/atomsk-imm-temp/src/atomsk"
    DFT_WORKFLOW_ROOT = Path(__file__).resolve().parent.parent / "SIMULATION"
    POTCAR_ROOT = "/work/wz300646/copied_POTCAR/PAW_PBE_54"