"""Geometry, state, mechanics, remodeling, and FF softening state machine."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math

import numpy as np
from hexalattice.hexalattice import create_hex_grid

from .config import ModelConfig


@dataclass
class ModelState:
    """All mutable simulation state, with no arrays shared between pair branches."""

    points: np.ndarray
    edges: np.ndarray
    hex_cells: np.ndarray
    boundary_node_ids: np.ndarray
    boundary_edge_ids: np.ndarray
    spring_constants: np.ndarray
    A: np.ndarray
    E: np.ndarray
    strain: np.ndarray
    stress: np.ndarray
    activation: np.ndarray
    D: np.ndarray
    initial_fibrotic: np.ndarray
    softening_active: np.ndarray
    ever_softened: np.ndarray
    recovery_active: np.ndarray
    soften_age: np.ndarray
    E_pre_softening: np.ndarray
    E_target: np.ndarray
    initial_cell_areas: np.ndarray

    def copy(self) -> "ModelState":
        return ModelState(**{name: getattr(self, name).copy() for name in self.__dataclass_fields__})


@dataclass(frozen=True)
class EquilibriumInfo:
    iterations: int
    converged_by_stagnation: bool
    accepted_uphill: int
    rejected_uphill: int
    final_energy: float


def build_geometry(config: ModelConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Reproduce the 0606 notebook's hex geometry and boundary labeling."""

    centers, _ = create_hex_grid(
        nx=config.nx, ny=config.ny, do_plot=False, min_diam=config.spacing
    )
    radius = config.spacing / math.sqrt(3.0)
    angles = np.linspace(0.0, 2.0 * np.pi, 7)[:-1]
    vertices = [
        np.column_stack(
            (
                center[0] + radius * np.cos(angles + np.pi / 2.0),
                center[1] + radius * np.sin(angles + np.pi / 2.0),
            )
        )
        for center in centers
    ]
    unique_vertices, indices = np.unique(
        np.vstack(vertices).round(decimals=6), axis=0, return_inverse=True
    )
    hex_cells = indices.reshape((-1, 6)).astype(np.int32)

    edges_list: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for hexagon in hex_cells:
        for index in range(6):
            edge = tuple(sorted((int(hexagon[index]), int(hexagon[(index + 1) % 6]))))
            if edge not in seen:
                edges_list.append(edge)
                seen.add(edge)
    edges = np.asarray(edges_list, dtype=np.int32)

    boundary_id = np.zeros(len(unique_vertices), dtype=np.int32)
    boundary_l: list[int] = []
    boundary_r: list[int] = []
    boundary_t: list[int] = []
    boundary_b: list[int] = []
    for point_id in range(len(unique_vertices)):
        occurrence = np.where(hex_cells == point_id)[1]
        if (1 in occurrence and 5 not in occurrence) or (2 in occurrence and 4 not in occurrence):
            boundary_l.append(point_id)
        if (5 in occurrence and 1 not in occurrence) or (4 in occurrence and 2 not in occurrence):
            boundary_r.append(point_id)
        if (0 in occurrence and 2 not in occurrence and 4 not in occurrence) or (
            1 in occurrence and 3 not in occurrence and 5 in occurrence
        ):
            boundary_t.append(point_id)
        if (3 in occurrence and 1 not in occurrence and 5 not in occurrence) or (
            4 in occurrence and 0 not in occurrence and 2 in occurrence
        ):
            boundary_b.append(point_id)
    boundary_id[np.unique(boundary_t)] = 1
    boundary_id[np.unique(boundary_b)] = 2
    boundary_id[np.unique(boundary_l)] = 3
    boundary_id[np.unique(boundary_r)] = 4

    boundary_status = np.zeros(len(edges), dtype=np.int32)
    for edge_index, (node1, node2) in enumerate(edges):
        if boundary_id[node1] > 0 and boundary_id[node2] > 0:
            boundary_status[edge_index] = 2
        elif (boundary_id[node1] > 0) != (boundary_id[node2] > 0):
            boundary_status[edge_index] = 1
    return unique_vertices.astype(np.float64), edges, hex_cells, boundary_id, boundary_status


def compute_strain_and_stress(
    points: np.ndarray, edges: np.ndarray, l0: float, spring_constants: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    vectors = points[edges[:, 1]] - points[edges[:, 0]]
    lengths = np.linalg.norm(vectors, axis=1)
    strain = (lengths - l0) / l0
    return strain, spring_constants * strain


def compute_total_energy(
    points: np.ndarray, edges: np.ndarray, l0: float, spring_constants: np.ndarray
) -> float:
    lengths = np.linalg.norm(points[edges[:, 1]] - points[edges[:, 0]], axis=1)
    return float(0.5 * np.sum(spring_constants * (lengths - l0) ** 2))


def compute_forces(
    points: np.ndarray, edges: np.ndarray, l0: float, spring_constants: np.ndarray
) -> np.ndarray:
    vectors = points[edges[:, 1]] - points[edges[:, 0]]
    lengths = np.linalg.norm(vectors, axis=1)
    if np.any(lengths == 0.0):
        raise RuntimeError("Zero-length spring encountered during force calculation.")
    magnitudes = spring_constants * (lengths - l0)
    edge_forces = (magnitudes / lengths)[:, None] * vectors
    forces = np.zeros_like(points)
    np.add.at(forces, edges[:, 0], edge_forces)
    np.add.at(forces, edges[:, 1], -edge_forces)
    return forces


def mechanical_equilibrate(
    points: np.ndarray,
    edges: np.ndarray,
    l0: float,
    spring_constants: np.ndarray,
    free_nodes: np.ndarray,
    config: ModelConfig,
    rng: np.random.Generator,
    max_iterations: int,
) -> tuple[np.ndarray, EquilibriumInfo]:
    """The 0606 annealed force-relaxation algorithm with an explicit RNG."""

    current = points.copy()
    mu = config.mu
    temperature = config.temperature
    previous_energy = compute_total_energy(current, edges, l0, spring_constants)
    stagnation_counter = 0
    accepted_uphill = 0
    rejected_uphill = 0

    for iteration in range(1, max_iterations + 1):
        forces = compute_forces(current, edges, l0, spring_constants)
        candidate = current.copy()
        candidate[free_nodes] += mu * forces[free_nodes]
        new_energy = compute_total_energy(candidate, edges, l0, spring_constants)
        delta_energy = new_energy - previous_energy
        denominator = max(abs(new_energy), np.finfo(float).tiny)
        if abs(delta_energy) / denominator < config.tolerance:
            stagnation_counter += 1
        else:
            stagnation_counter = 0

        # Always consume one uniform per mechanical iteration.  This permits
        # paired runs to share a phase/iteration-indexed Metropolis noise
        # field even where their energy trajectories differ.
        acceptance_draw = rng.random()
        if delta_energy < 0.0:
            current = candidate
            previous_energy = new_energy
            mu *= 1.1
        else:
            acceptance_probability = math.exp(-delta_energy / max(temperature, np.finfo(float).tiny))
            if acceptance_draw < acceptance_probability:
                current = candidate
                previous_energy = new_energy
                accepted_uphill += 1
            else:
                mu *= 0.9
                rejected_uphill += 1
        temperature *= 0.99

        if stagnation_counter >= config.stagnation_window:
            return current, EquilibriumInfo(
                iteration, True, accepted_uphill, rejected_uphill, previous_energy
            )
    return current, EquilibriumInfo(
        max_iterations, False, accepted_uphill, rejected_uphill, previous_energy
    )


def polygon_areas(points: np.ndarray, hex_cells: np.ndarray) -> np.ndarray:
    """Absolute shoelace areas of the stored, deformed hexagonal cell polygons."""

    polygons = points[hex_cells]
    x = polygons[:, :, 0]
    y = polygons[:, :, 1]
    return 0.5 * np.abs(np.sum(x * np.roll(y, -1, axis=1) - y * np.roll(x, -1, axis=1), axis=1))


def make_initial_state(
    config: ModelConfig,
    biology_rng: np.random.Generator,
    shared_mechanics_rng: np.random.Generator,
) -> tuple[ModelState, EquilibriumInfo]:
    """Generate the common post-preconvergence state for one paired experiment."""

    points, edges, hex_cells, boundary_nodes, boundary_edges = build_geometry(config)
    spring_constants = np.full(len(edges), config.k_normal, dtype=np.float64)
    healthy = np.where(spring_constants == config.k_normal)[0]
    while len(healthy) / len(edges) > 1.0 - config.alpha:
        feasible = np.where(spring_constants == config.k_normal)[0]
        if len(feasible) == 0:
            break
        spring_index = int(biology_rng.choice(feasible))
        spring_constants[spring_index] = config.k_fibrosis
        for _ in range(config.random_walk_length):
            nodes = edges[spring_index]
            connected = np.where(np.any(np.isin(edges, nodes), axis=1))[0]
            spring_index = int(biology_rng.choice(connected))
            spring_constants[spring_index] = config.k_fibrosis
        healthy = np.where(spring_constants == config.k_normal)[0]

    initial_fibrotic = spring_constants == config.k_fibrosis
    free_nodes = np.where(boundary_nodes == 0)[0]
    points, pre_info = mechanical_equilibrate(
        points,
        edges,
        config.l0,
        spring_constants,
        free_nodes,
        config,
        shared_mechanics_rng,
        config.preconvergence_iterations,
    )
    strain, stress = compute_strain_and_stress(points, edges, config.l0, spring_constants)
    n_springs = len(edges)
    initial_areas = polygon_areas(points, hex_cells)
    state = ModelState(
        points=points,
        edges=edges,
        hex_cells=hex_cells,
        boundary_node_ids=boundary_nodes,
        boundary_edge_ids=boundary_edges,
        spring_constants=spring_constants.copy(),
        # Preserves 0606's A = spring_constants assignment after preconvergence.
        A=spring_constants.copy(),
        E=np.full(n_springs, config.E_initial, dtype=np.float64),
        strain=strain,
        stress=stress,
        activation=np.zeros(n_springs, dtype=np.float64),
        D=np.full(n_springs, config.D0, dtype=np.float64),
        initial_fibrotic=initial_fibrotic,
        softening_active=np.zeros(n_springs, dtype=bool),
        ever_softened=np.zeros(n_springs, dtype=bool),
        recovery_active=np.zeros(n_springs, dtype=bool),
        soften_age=np.zeros(n_springs, dtype=np.int32),
        E_pre_softening=np.full(n_springs, np.nan, dtype=np.float64),
        E_target=np.full(n_springs, np.nan, dtype=np.float64),
        initial_cell_areas=initial_areas,
    )
    return state, pre_info


def update_activation_and_density(state: ModelState, config: ModelConfig) -> None:
    """Preserve the 0606 activation and deterministic density equations."""

    state.strain, state.stress = compute_strain_and_stress(
        state.points, state.edges, config.l0, state.spring_constants
    )
    a_k = (state.spring_constants / config.At) ** config.beta / (
        (state.spring_constants / config.At) ** config.beta + 1.0
    )
    a_epsilon = (state.strain / config.epsilon_s) ** config.gamma / (
        (state.strain / config.epsilon_s) ** config.gamma + 1.0
    )
    state.activation = state.activation * (1.0 - config.memory_factor) + config.memory_factor * (
        config.w1 * a_epsilon + config.w2 * a_k - config.activation_offset
    )

    # Area remodeling remains separate from the FF intervention.
    state.A = state.A + config.area_update_scale * state.activation * state.D

    X = config.Dmax / config.D0
    p2 = config.p1 * len(state.points) / (len(state.edges) * config.D0)
    p3 = (p2 * (1.0 - 1.0 / X)) / (config.w2 - 0.5)
    state.D = state.D + (config.p1 * len(state.points) / len(state.edges)) + (
        p3 * state.activation - p2
    ) * state.D
    state.D = np.clip(state.D, config.D0, config.Dmax)


def apply_softening(state: ModelState, config: ModelConfig, enable_softening: bool) -> None:
    """Apply one event at most once per eligible spring.

    At trigger, E falls immediately from its current value to
    ``max(minimum_E, pre_E * (1 - magnitude))``.  It is held for exactly
    ``duration`` outer agent phases, then returns geometrically toward the
    recorded pre-event E.  Recovery terminates using an explicit tolerance,
    rather than floating-point equality.  A is never modified here.
    """

    soft = config.softening
    if enable_softening:
        # Age pre-existing events before adding new ones, so a duration of N
        # holds the target E through N complete outer phases.
        active = state.softening_active
        if np.any(active):
            state.soften_age[active] += 1
            finished = active & (state.soften_age >= soft.duration)
            state.softening_active[finished] = False
            state.recovery_active[finished] = True

        recovering = state.recovery_active
        if np.any(recovering):
            state.E[recovering] += soft.recovery_rate * (
                state.E_pre_softening[recovering] - state.E[recovering]
            )
            restored = recovering & (
                np.abs(state.E_pre_softening - state.E) <= soft.recovery_tolerance
            )
            state.E[restored] = state.E_pre_softening[restored]
            state.recovery_active[restored] = False

        eligible = (state.spring_constants > soft.threshold) & ~state.ever_softened
        if np.any(eligible):
            state.ever_softened[eligible] = True
            state.softening_active[eligible] = True
            state.recovery_active[eligible] = False
            state.soften_age[eligible] = 0
            state.E_pre_softening[eligible] = state.E[eligible]
            state.E_target[eligible] = np.maximum(
                soft.minimum_E, state.E[eligible] * (1.0 - soft.magnitude)
            )
            state.E[eligible] = state.E_target[eligible]

    # Exact archival equation, now allowing E to causally alter k.
    state.spring_constants = state.E * state.A


def state_fingerprint(state: ModelState) -> str:
    """Stable content hash used to prove common pair initialization."""

    digest = hashlib.sha256()
    for name in (
        "points",
        "edges",
        "hex_cells",
        "boundary_node_ids",
        "boundary_edge_ids",
        "spring_constants",
        "A",
        "E",
        "strain",
        "stress",
        "activation",
        "D",
        "initial_fibrotic",
        "softening_active",
        "ever_softened",
        "recovery_active",
        "soften_age",
        "E_pre_softening",
        "E_target",
        "initial_cell_areas",
    ):
        array = np.ascontiguousarray(getattr(state, name))
        digest.update(name.encode("utf-8"))
        digest.update(str(array.dtype).encode("utf-8"))
        digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
        digest.update(array.tobytes())
    return digest.hexdigest()
