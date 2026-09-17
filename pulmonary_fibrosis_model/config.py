"""Configuration for the recovered 0606 model equations."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from typing import Any


@dataclass(frozen=True)
class SofteningConfig:
    """FF-like E-based intervention parameters; not biologically calibrated."""

    threshold: float = 5.0
    # Carries forward the visible 0606 softening_factor value, although it was
    # inert there and must not be interpreted as a calibrated FF magnitude.
    magnitude: float = 0.05
    duration: int = 150
    recovery_rate: float = 0.003
    recovery_tolerance: float = 1e-8
    minimum_E: float = 0.1


@dataclass(frozen=True)
class ModelConfig:
    """Centralized parameters for the archival Charline/0606 model.

    ``baseline_config`` deliberately retains the visible 0606 values.  In
    particular, the historical implementation uses k = E*A (not EA/l0),
    initializes A from spring stiffness after fibrosis initialization, and
    feeds spring stiffness to the activation function named ``calc_aK``.
    Those semantics are preserved here pending scientific review.
    """

    nx: int = 25
    ny: int = 40
    spacing: float = 1.0
    k_normal: float = 1.0
    k_fibrosis: float = 50.0
    l0_ratio: float = 0.6
    alpha: float = 0.6
    random_walk_length: int = 300

    mu: float = 0.01
    temperature: float = 1.0
    tolerance: float = 1e-7
    preconvergence_iterations: int = 2000
    agent_update_interval: int = 1500
    agent_total_iterations: int = 600
    stagnation_window: int = 20

    At: float = 10.0
    A0: float = 1.0
    epsilon_s: float = 1.0
    beta: float = 3.0
    gamma: float = 1.0
    w1: float = 1.0
    w2: float = 1.0
    memory_factor: float = 1.0
    area_update_scale: float = 0.2
    p1: float = 0.01
    D0: float = 0.3
    Dmax: float = 3.0
    E_initial: float = 1.0

    # Compact phase summaries are retained at every outer agent phase.
    output_interval: int = 1
    spatial_output_interval: int = 25
    softening: SofteningConfig = field(default_factory=SofteningConfig)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def l0(self) -> float:
        # Exact 0606 expression: l0_ratio * (2/sqrt(3)*spacing/2).
        import math

        return self.l0_ratio * self.spacing / math.sqrt(3.0)

    @property
    def activation_offset(self) -> float:
        return 0.5 * self.w1 + self.w2 * (self.A0 / self.At) ** self.beta / (
            (self.A0 / self.At) ** self.beta + 1.0
        )


def baseline_config() -> ModelConfig:
    """Return the archival 0606 parameterization, without altering its values."""

    return ModelConfig()


def smoke_config() -> ModelConfig:
    """A deliberately small, non-biological configuration for end-to-end tests."""

    return replace(
        baseline_config(),
        nx=4,
        ny=5,
        random_walk_length=12,
        preconvergence_iterations=60,
        agent_update_interval=40,
        agent_total_iterations=6,
        output_interval=1,
        spatial_output_interval=1,
        softening=SofteningConfig(
            threshold=5.0,
            magnitude=0.5,
            duration=2,
            recovery_rate=0.5,
            recovery_tolerance=1e-8,
            minimum_E=0.1,
        ),
    )
