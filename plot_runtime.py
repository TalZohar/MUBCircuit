import time
import random

import mub_circuit
from typing import Callable, List

from matplotlib import pyplot as plt


def plot_runtime(func: Callable[[int], None], ns: List[int]) -> None:
    """Measure runtime of `func(n)` for each n in `ns` and plot it against n³."""
    runtimes = []

    for n in ns:
        start = time.perf_counter()
        func(n)
        end = time.perf_counter()
        runtimes.append(end - start)

    # Plot measured runtimes vs. ideal cubic curve (normalized)
    plt.plot(ns, runtimes, marker='o', label="Measured runtime")

    plt.xlabel("n")
    plt.ylabel("Time (seconds)")
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_circuit_construction_runtime() -> None:
    def random_circuit_generator(n: int) -> None:
        mub_circuit.QUBIT_NUM = n
        mub_circuit.mub_circuit(n, random.randrange(n))  # provide random j

    plot_runtime(random_circuit_generator, list(range(1, 15)))

if __name__ == '__main__':
    plot_circuit_construction_runtime()