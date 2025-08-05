from typing import List
import numpy as np
from qiskit import QuantumCircuit

from consts import IRREDUCIBLE_POLYS, QUBIT_NUM


def to_base_p(x: int, p: int, n: int) -> np.ndarray:
    digits = []
    for _ in range(n):
        digits.append(x % p)
        x //= p
    return np.array(digits, dtype=int)


def binary_to_vector(x: int) -> np.ndarray:
    return np.fromiter(f'{x:b}', dtype=int)[::-1]


def calculate_galois_matrices() -> List[np.ndarray]:
    irreducible_poly = IRREDUCIBLE_POLYS[QUBIT_NUM]
    x_vals = []
    for index in range(QUBIT_NUM):
        x = np.zeros(QUBIT_NUM, dtype=int)
        x[index] = 1
        x_vals.append(x)
    for extended_index in range(QUBIT_NUM - 1):
        prev_x = x_vals[-1]
        x = np.concatenate(([0], prev_x[:-1]))
        x += prev_x[-1] * binary_to_vector(irreducible_poly)[:-1]
        x_vals.append(x)

    galois_matrices = []
    for index in range(QUBIT_NUM):
        galois_matrix = np.zeros((QUBIT_NUM, QUBIT_NUM), dtype=int)
        for s in range(QUBIT_NUM):
            for t in range(QUBIT_NUM):
                galois_matrix[s][t] = x_vals[s + t][index]
        galois_matrices.append(galois_matrix)
    return galois_matrices


def gf_mul(a: int, b: int, galois_matrices: List[np.ndarray], lsb_bits: int = QUBIT_NUM) -> int:
    lsb_bits = min(QUBIT_NUM, lsb_bits)
    out = 0
    a_v = to_base_p(a, 2, QUBIT_NUM)
    b_v = to_base_p(b, 2, QUBIT_NUM)
    for index in range(lsb_bits):
        galois_matrix = galois_matrices[index]
        index_val = (a_v.T @ galois_matrix @ b_v) % 2
        out += index_val.item() * (1 << index)
    return out


def calculate_a(j: int, galois_matrices: List[np.ndarray]) -> np.ndarray:
    a_arr = np.zeros(QUBIT_NUM, dtype=int)
    for index in range(QUBIT_NUM):
        inner_term = gf_mul(1 << index, 1 << index, galois_matrices, lsb_bits=QUBIT_NUM)
        exponent = gf_mul(j, inner_term, galois_matrices, lsb_bits=2)
        a_arr[index] = exponent if exponent % 2 == 0 else 4 - exponent  # take the conjugate into account
    return a_arr


def calculate_b(j: int, galois_matrices: List[np.ndarray]) -> np.ndarray:
    b_arr = np.zeros(2 * QUBIT_NUM - 1, dtype=int)
    x = 1
    for index in range(2 * QUBIT_NUM - 1):
        if index > 0:
            x = gf_mul(x, 2, galois_matrices)
        b_arr[index] = gf_mul(j, x, galois_matrices, lsb_bits=1)
    return b_arr


def mub_circuit(n: int, j: int) -> QuantumCircuit:
    """
    Build the MUB circuit U(j) for n qubits.
    The circuit does H^⊗n, then S^a on each qubit, then CZ gates.
    """
    qc = QuantumCircuit(n, name=f"MUB generator for j={j}")
    # H-part: Hadamard on all qubits
    for q in range(n):
        qc.h(q)
    # Compute S and CZ parameters for this j
    galois_matrices = calculate_galois_matrices()
    a = calculate_a(j, galois_matrices)
    b = calculate_b(j, galois_matrices)
    print(f"a: {a}\n b: {b}")
    # S-part: apply S^a_t on qubit t
    for t in range(n):
        # a[t] == 0: do nothing
        if a[t] == 1:
            qc.s(t)
        elif a[t] == 2:
            qc.z(t)
        elif a[t] == 3:
            qc.sdg(t)
    # CZ-part: apply CZ(s,t) where b[(s,t)] == 1
    for s in range(QUBIT_NUM):
        for t in range(s + 1, QUBIT_NUM):
            if b[s + t] == 1:
                qc.cz(s, t)
    return qc
