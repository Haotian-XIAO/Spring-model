"""Reproducible paired experiments for the archival 0606 spring-network model."""

from .config import ModelConfig, SofteningConfig, baseline_config, smoke_config
from .experiment import generate_initial_state, run_batch, run_null_pair, run_pair

__all__ = [
    "ModelConfig",
    "SofteningConfig",
    "baseline_config",
    "smoke_config",
    "generate_initial_state",
    "run_pair",
    "run_null_pair",
    "run_batch",
]
