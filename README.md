# ket_operation_simpy_simulation
A symbolic framework for simulating pure and mixed quantum states, gates, measurements, and entanglement in Python with sympy.

---

## Features

- Symbolic manipulation of **pure states** and **mixed states**
- Built-in **quantum gates**: Pauli (X, Y, Z), Hadamard, CNOT, rotation gates
- **Measurement** and projection (including mixed-state measurements)
- **Tensor product** and **partial trace**
- **Fidelity** and **Schmidt decomposition**
- Exact **density matrix** representation from mixed states
- Protocol demos: **quantum teleportation**, **entanglement distillation**, and more

---

## Data Structure

This simulator uses clean symbolic Python lists to represent quantum states:

- `pure_state`:  
  A list of `(coefficient, basis)` tuples, where `basis` is a tuple like `(0, 1, 1)`  
  Example:  
  ```python
  [(1/sqrt(2), (0, 0)), (1/sqrt(2), (1, 1))]
  ```

- `mixed_state`:  
  A list of `(probability, pure_state)` entries.  
  Example:  
  ```python
  [(1/2, [(1, (0,))]), (1/2, [(1, (1,))])]
  ```

All operations (gates, measurement, fidelity, etc.) work on these structures.

---

## Requirements

- Python ≥ 3.8
- [`sympy`](https://www.sympy.org/)

Install dependencies with:

```bash
pip install sympy
```

---

## How to Use

Everything is in the main notebook:

> ✅ `Ket_operations.ipynb` — contains all symbolic definitions and demos (teleportation, fidelity distillation, Schmidt decomposition, etc.)

You can open and run it locally or in Colab.

---

## Acknowledgment

> Part of this project was created with the assistance of [ChatGPT](https://chat.openai.com/).
