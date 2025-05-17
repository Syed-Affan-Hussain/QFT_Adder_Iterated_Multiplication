# QFT Quantum Addition Algorithm - Corrected Implementation
import numpy as np
from numpy import pi
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, Aer, execute
import timeit

# _____________________________________________________________________________________________________________________________________________________________________________________________________________________
def qft(circuit, qubits):
    """Applies the Quantum Fourier Transform to a register of qubits."""
    n = len(qubits)
    for j in range(n):
        # Apply Hadamard gate
        circuit.h(qubits[j])
        # Apply controlled phase rotations
        for k in range(j + 1, n):
            # Angle is pi / 2^(k-j)
            angle = pi / (2**(k - j))
            circuit.cp(angle, qubits[k], qubits[j]) # Controlled phase rotation CPhase(theta) = diag(1, 1, 1, exp(i*theta))

    # Swap qubits to reverse the order (standard QFT output)
    for j in range(n // 2):
        circuit.swap(qubits[j], qubits[n - 1 - j])

# __________________________________________________________________________________________________________________________________________________________________________________________________________________
def inverse_qft(circuit, qubits):
    """Applies the Inverse Quantum Fourier Transform to a register of qubits."""
    n = len(qubits)
    # Swap qubits to reverse the order (undo the swap from QFT)
    for j in range(n // 2):
        circuit.swap(qubits[j], qubits[n - 1 - j])

    for j in reversed(range(n)):
        # Apply controlled phase rotations (negative angles for inverse)
        for k in reversed(range(j + 1, n)):
            # Angle is -pi / 2^(k-j)
            angle = -pi / (2**(k - j))
            circuit.cp(angle, qubits[k], qubits[j])
        # Apply Hadamard gate
        circuit.h(qubits[j])

# __________________________________________________________________________________________________________________________________________________________________________________________________________________
def controlled_phase_adder(circuit, reg_a, reg_b, n):
    """
    Applies controlled phase rotations to add reg_b to reg_a in the phase basis.
    Assumes reg_a has n+1 qubits and reg_b has n qubits for adding n-bit numbers
    and storing the n+1 bit sum in reg_a.
    The controlled rotations are applied between reg_b[k] (control) and reg_a[j] (target).
    Angle for CPhase(theta) between b_k and a_j is theta = pi / 2^(j-k).
    """
    # Iterate through qubits of reg_a (target)
    for j in range(n + 1):
        # Iterate through qubits of reg_b (control)
        for k in range(n):
            # Ensure k < j for the standard phase addition circuit structure
            # and that the control qubit index k is valid for reg_b
            if j - k > 0 and k < n:
                 # Angle for CPhase(theta) between b_k and a_j is theta = pi / 2^(j-k).
                angle = pi / (2**(j - k))
                # Apply controlled phase rotation: controlled by reg_b[k], target reg_a[j]
                circuit.cp(angle, reg_b[k], reg_a[j])

# __________________________________________________________________________________________________________________________________________________________________________________________________________________
def qft_adder_circuit(num_qubits_a, num_qubits_b):
    """
    Creates the quantum circuit for the QFT adder.
    Adds an n-bit number (in reg_b) to an n-bit number (in reg_a),
    storing the (n+1)-bit sum in reg_a.
    Requires num_qubits_a = num_qubits_b + 1.
    """
    # reg_a will hold the sum, needs num_qubits_b + 1 qubits
    reg_a = QuantumRegister(num_qubits_a, "reg_a")
    # reg_b holds the number to be added
    reg_b = QuantumRegister(num_qubits_b, "reg_b")
    # Classical register to measure the sum
    classic_reg = ClassicalRegister(num_qubits_a, "classic_reg")

    # Create the quantum circuit
    circuit = QuantumCircuit(reg_a, reg_b, classic_reg, name="qft_adder")

    # --- QFT Adder Algorithm Steps ---

    # 1. Apply QFT to reg_a
    # Note: The standard QFT adder applies QFT to the register that will hold the sum (reg_a).
    # The initial state of reg_a should represent the first number.
    # The initial state of reg_b should represent the second number.
    # The QFT is applied AFTER initializing reg_a with the first number.
    # We will handle initial state preparation outside this function.

    # 2. Apply controlled phase additions using reg_b to modify reg_a in the phase basis
    # This effectively adds the value of reg_b to the phase representation of reg_a
    # The controlled_phase_adder function assumes reg_a is in the phase basis (after QFT)
    # and reg_b is in the computational basis.
    # We need to apply QFT to reg_a *before* calling controlled_phase_adder.
    # Let's structure the main execution flow to handle this.

    # The circuit construction itself just defines the gates.
    # The actual QFT and inverse QFT calls will be in the main execution block
    # after initializing the qubits.

    return circuit, reg_a, reg_b, classic_reg

# __________________________________________________________________________________________________________________________________________________________________________________________________________________
# Helper functions for binary/decimal conversion
def bin_to_dec(b):
    """Converts a binary string to a decimal integer."""
    return int(b, 2)

def dec_to_bin(u, num_bits):
    """
    Converts a decimal integer to a binary string of a specified length.
    Pads with leading zeros if necessary.
    """
    if u < 0:
        raise ValueError("Input must be a non-negative integer")
    binary_string = bin(u)[2:] # Convert to binary string, remove '0b' prefix
    if len(binary_string) > num_bits:
         raise ValueError(f"Integer {u} requires more than {num_bits} bits.")
    # Pad with leading zeros to match the desired number of bits
    padded_binary_string = binary_string.zfill(num_bits)
    return padded_binary_string

# _____________________________________________________________________________________________________________________________________________________________________________________________________________________
# Main execution block
if __name__ == "__main__":
    # Get integer inputs from the user
    try:
        int1 = int(input("Enter First non-negative Integer to add: "))
        int2 = int(input("Enter Second non-negative Integer to add: "))
        if int1 < 0 or int2 < 0:
             print("Please enter non-negative integers.")
             exit()
    except ValueError:
        print("Invalid input. Please enter integers.")
        exit()

    # Determine the number of qubits needed
    # The sum of two n-bit numbers can be up to n+1 bits.
    # We need enough qubits in reg_a to store the sum.
    # Let n be the number of bits required for the larger input number.
    # reg_b will have n qubits, reg_a will have n+1 qubits.
    max_int = max(int1, int2)
    # Number of bits for the larger number
    num_bits_b = len(bin(max_int)[2:]) if max_int > 0 else 1 # At least 1 bit

    # reg_a needs one more qubit for the potential carry
    num_qubits_a = num_bits_b + 1
    num_qubits_b = num_bits_b # reg_b size

    # Convert integers to binary strings with appropriate padding
    # reg_b will hold the n-bit number (padded to num_bits_b)
    # reg_a will hold the first n-bit number (padded to num_bits_b),
    # but is a (n+1)-bit register. We initialize the first n bits of reg_a.
    bs1 = dec_to_bin(int1, num_bits_b)
    bs2 = dec_to_bin(int2, num_bits_b)

    print(f"\nAdding {int1} ({bs1}) and {int2} ({bs2})")
    print(f"Using {num_qubits_a} qubits for reg_a and {num_qubits_b} for reg_b.")

    # Start timer
    start = timeit.default_timer()

    # Create the circuit and registers
    circuit, reg_a, reg_b, classic_reg = qft_adder_circuit(num_qubits_a, num_qubits_b)

    # --- Prepare initial state ---
    # Initialize reg_a with the binary representation of the first integer (int1)
    # The LSB of the binary string corresponds to the first qubit (index 0) in Qiskit registers.
    # We need to reverse the binary string for correct mapping to qubits.
    reversed_bs1 = bs1[::-1]
    for i in range(num_bits_b): # Initialize the first num_bits_b qubits of reg_a
        if reversed_bs1[i] == '1':
            circuit.x(reg_a[i])

    # Initialize reg_b with the binary representation of the second integer (int2)
    reversed_bs2 = bs2[::-1]
    for i in range(num_bits_b): # Initialize all qubits of reg_b
        if reversed_bs2[i] == '1':
            circuit.x(reg_b[i])

    circuit.barrier() # Separator for visualization

    # --- Apply the QFT Adder Algorithm ---

    # 1. Apply QFT to reg_a
    qft(circuit, reg_a)
    circuit.barrier()

    # 2. Apply controlled phase additions using reg_b to reg_a
    # This adds the value of reg_b to reg_a in the phase basis.
    controlled_phase_adder(circuit, reg_a, reg_b, num_bits_b) # Pass num_bits_b as n for controlled_phase_adder
    circuit.barrier()

    # 3. Apply Inverse QFT to reg_a
    inverse_qft(circuit, reg_a)
    circuit.barrier()

    # --- Measurement ---
    # Measure reg_a (which now holds the sum) into the classical register
    circuit.measure(reg_a, classic_reg)

    # --- Simulation ---
    # Use the qasm_simulator backend
    backend = Aer.get_backend('qasm_simulator')

    # Execute the circuit
    # We only need 1 shot to get the result as the adder is deterministic
    job = execute(circuit, backend, shots=1)
    result = job.result()
    counts = result.get_counts(circuit)

    # --- Interpret Results ---
    # The measurement result (key in the counts dictionary) is a binary string.
    # The bits in the measurement string are ordered from MSB to LSB,
    # corresponding to the classical register bits.
    # Qiskit measures qubits into classical bits in the order of the classical register indices.
    # The classical register classic_reg[i] corresponds to reg_a[i].
    # The QFT output is typically reversed, but our `qft` and `inverse_qft` include swaps.
    # The classical register stores the result in MSB...LSB order from left to right.
    # So the measurement string '101' means classic_reg[0]=1, classic_reg[1]=0, classic_reg[2]=1.
    # If reg_a[0] is LSB and reg_a[n] is MSB, and classic_reg[i] measures reg_a[i],
    # then the measurement string is MSB...LSB.
    # We need to take the first key from the counts dictionary, which is the measured binary string.
    measured_binary_sum_msb_lsb = list(counts.keys())[0]

    # Convert the binary string (MSB...LSB) to decimal
    decimal_sum = bin_to_dec(measured_binary_sum_msb_lsb)

    print("\nCircuit Diagram:")
    print(circuit.draw())

    print(f"\nMeasured binary sum (MSB...LSB): {measured_binary_sum_msb_lsb}")
    print(f"Decimal sum: {decimal_sum}")

    # Check the result classically
    classical_sum = int1 + int2
    print(f"Classical sum: {classical_sum}")

    if decimal_sum == classical_sum:
        print("Quantum addition result matches classical result.")
    else:
        print("Mismatch between quantum and classical addition results.")


    # Stop timer
    stop = timeit.default_timer()
    execution_time = stop - start

    print(f"\nProgram Executed in {execution_time:.6f} seconds")

