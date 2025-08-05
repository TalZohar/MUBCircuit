Run `python generate_mubs.py` to create and plot the MUB circuits.
* These circuits are of `QUBIT_NUM` size, defined in `consts.py`.
* The program produces `2**QUBIT_NUM + 1` circuits, where each circuit produces an orthonormal basis (each basis vector corresponds to a different standard input).
* The produced bases of the circuits are mutually unbiased

Run `python plot_runtime.py` to see the runtime of calculating a random circuit as a function of `QUBIT_NUM`.
This should display a runtime graph of `O(n**3)`, like the efficient implementation described in the paper.