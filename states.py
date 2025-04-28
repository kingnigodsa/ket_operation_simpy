from sympy import Matrix, simplify, KroneckerProduct

class PureState:
    """
    A quantum pure state |ψ⟩ represented as a SymPy column vector.
    """
    def __init__(self, amplitudes: Matrix):
        """
        amplitudes: a column Matrix of symbolic or numeric entries.
        """
        if amplitudes.cols != 1:
            raise ValueError("Amplitude must be a column vector.")
        self.vec = simplify(amplitudes)

    @classmethod
    def from_tuples(cls, terms, dim=None):
        """
        Build a PureState from a list of (coefficient, basis_index) pairs.
        - terms: list of (coeff, int) where 0 <= int < dim.
        - dim: dimension of the Hilbert space. If None, inferred as max index + 1.
        """
        if dim is None:
            dim = max(idx for _, idx in terms) + 1
        amp = [0]*dim
        for coeff, idx in terms:
            amp[idx] = coeff
        return cls(Matrix(amp).reshape(dim, 1))

    def apply_gate(self, gate: Matrix) -> "PureState":
        """
        Return a new PureState = gate * |ψ⟩.
        Gate must be a square Matrix of same dimension.
        """
        if gate.shape != (self.vec.rows, self.vec.rows):
            raise ValueError("Gate matrix must match state dimension.")
        return PureState(gate * self.vec)

    def tensor(self, other: "PureState") -> "PureState":
        """
        Return the tensor (Kronecker) product |ψ⟩ ⊗ |φ⟩.
        """
        return PureState(KroneckerProduct(self.vec, other.vec))

    def to_mixed(self, prob=1) -> "MixedState":
        """
        Promote this pure state to a MixedState with given weight.
        """
        from states import MixedState
        return MixedState([(prob, self)])

    def inner(self, other: "PureState"):
        """
        Compute ⟨ψ|φ⟩.
        """
        return simplify(self.vec.conjugate().T * other.vec)

    def normalize(self) -> "PureState":
        """
        Return a normalized copy of this state.
        """
        norm = simplify(simplify((self.vec.conjugate().T * self.vec)[0]))
        return PureState(self.vec / norm**0.5)


class MixedState:
    """
    A mixed quantum state ρ = ∑ p_i |ψᵢ⟩⟨ψᵢ|.
    """
    def __init__(self, ensemble):
        """
        ensemble: list of (probability, PureState) pairs.
        """
        total = sum(p for p, _ in ensemble)
        if simplify(total - 1) != 0:
            raise ValueError("Probabilities must sum to 1.")
        self.ensemble = ensemble

    def density_matrix(self) -> Matrix:
        """
        Build and return the density matrix ρ = ∑ p |ψ⟩⟨ψ|.
        """
        ρ = Matrix.zeros(self.ensemble[0][1].vec.rows)
        for p, ψ in self.ensemble:
            ρ += p * (ψ.vec * ψ.vec.conjugate().T)
        return simplify(ρ)

    def measure(self, projectors: list[Matrix]):
        """
        Perform a projective measurement.
        - projectors: list of d×d Hermitian idempotent Matrices that sum to I.
        Returns (outcome_index, post_state: MixedState).
        """
        ρ = self.density_matrix()
        probs = [simplify((proj * ρ).trace()) for proj in projectors]
        # Symbolic probabilities may sum to 1; user must substitute numerics before sampling.
        # Here we just pick the first non-zero branch for demonstration:
        for i, p in enumerate(probs):
            if p != 0:
                ρ_post = simplify(projectors[i] * ρ * projectors[i])
                ρ_post = ρ_post / simplify(ρ_post.trace())
                # Decompose back into pure ensemble via eigen-decomposition
                evals, evecs = ρ_post.diagonalize()
                new_ensemble = []
                for ev, evec in zip(evals, evecs):
                    if ev != 0:
                        vec = evec.col_insert(0, Matrix([]))  # extract column
                        new_ensemble.append((ev, PureState(vec)))
                return i, MixedState(new_ensemble)
        raise RuntimeError("All outcome probabilities are zero.")

    def partial_trace(self, subsys_dims, keep):
        """
        Compute the partial trace over unwanted subsystems.
        - subsys_dims: list of dimensions (e.g. [2,2] for two qubits).
        - keep: list of subsystem indices to retain (0-based).
        Returns a new MixedState on the reduced space.
        """
        from utils import partial_trace_ensemble  # implement separately
        new_ensemble = partial_trace_ensemble(self.ensemble, subsys_dims, keep)
        return MixedState(new_ensemble)
