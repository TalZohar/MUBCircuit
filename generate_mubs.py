import io
import itertools
import math
from typing import List

import numpy as np
from PIL import Image
from matplotlib import pyplot as plt
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Statevector
from qiskit.visualization import circuit_drawer, plot_histogram
from qiskit_aer import Aer

from mub_circuit import QUBIT_NUM, mub_circuit


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
