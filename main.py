import itertools
from typing import List
import numpy as np
import math
from matplotlib import pyplot as plt
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator, Aer
from qiskit.compiler import transpile
from qiskit.visualization import plot_histogram, circuit_drawer
from qiskit.quantum_info import Statevector
from PIL import Image
import io

# Known irreducible polynomials (binary) for small n:
IRREDUCIBLE_POLYS = {1: 0b11,
                     2: 0b111,  # x^2 + x + 1
                     3: 0b1011,  # x^3 + x + 1
                     4: 0b10011,  # x^4 + x + 1
                     5: 0b100101,  # x^5 + x^2 + 1
                     6: 0b1000011,  # x^6 + x + 1
                     7: 0b10000011,  # etc.
                     8: 0b100011011,
                     9: 0b1000000101,
                     10: 0b10000001001}

QUBIT_NUM = 6


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


def validate_orthogonality(base: List[Statevector]) -> bool:
    for i in range(len(base)):
        for j in range(i + 1, len(base)):
            dot_prod = base[i].data.conj().dot(base[j].data)
            if not np.isclose(dot_prod, 0):
                return False
    return True


def validate_mubs(base1: List[Statevector], base2: List[Statevector]) -> bool:
    assert len(base1) == len(base2)
    dim = len(base1)
    for k in range(dim):
        for l in range(dim):
            dot_prod = base1[k].data.conj().dot(base2[l].data)
            if not np.isclose(abs(dot_prod ** 2), 1 / dim):
                return False
    return True


def draw_circuits_grid(circuits: List[QuantumCircuit]) -> None:
    count = len(circuits)
    cols = math.ceil(math.sqrt(count))
    rows = math.ceil(count / cols)

    fig, axs = plt.subplots(rows, cols, figsize=(cols * 4, rows * 3))
    axs = axs.flatten() if isinstance(axs, (list, np.ndarray)) else [axs]

    for ax in axs[count:]:
        ax.axis("off")

    for i, circuit in enumerate(circuits):
        # Render circuit to matplotlib figure
        fig_circuit = circuit_drawer(circuit, output='mpl')

        # Save to buffer
        buf = io.BytesIO()
        fig_circuit.savefig(buf, format='png')
        buf.seek(0)
        img = Image.open(buf)

        # Display as image
        axs[i].imshow(img)
        axs[i].axis("off")
        axs[i].set_title(circuit.name or f"Circuit {i}")

        plt.close(fig_circuit)

    plt.tight_layout()
    plt.show()


def generate_mubs(run_simultion: bool = False) -> None:
    standard_basis_inputs = [''.join(bits) for bits in itertools.product('01', repeat=QUBIT_NUM)]
    standard_basis = [Statevector.from_label(i) for i in standard_basis_inputs]
    standard_basis_circuit = QuantumCircuit(QUBIT_NUM, name="MUB generator for standard basis")
    mubs = [standard_basis]
    circuits = [standard_basis_circuit]
    for j in range(2 ** QUBIT_NUM):
        qc = mub_circuit(QUBIT_NUM, j)
        circuits.append(qc)

        # Find MUB
        new_base = []
        for circuit_input in standard_basis:
            output_sv = circuit_input.evolve(qc)
            new_base.append(output_sv)
        # Validate result
        assert validate_orthogonality(new_base)
        for base in mubs:
            assert validate_mubs(base, new_base)

        mubs.append(new_base)

        if run_simultion:
            backend = Aer.get_backend('statevector_simulator')
            compiled_circuit = transpile(qc, backend)
            job = backend.run(compiled_circuit, shots=1024)
            result = job.result()
            counts = result.get_counts()
            # Plot the results
            plot_histogram(counts)

    draw_circuits_grid(circuits)


if __name__ == '__main__':
    generate_mubs()
