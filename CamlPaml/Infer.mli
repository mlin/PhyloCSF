(** Probabilistic inference on phylogenetic trees given the substitution (P) matrices on each branch. These are low-level routines that should usually be used through the higher-level {{:PhyloModel}[PhyloModel]} abstractions. *)

type intermediate

type leaf = [`Certain of int | `Distribution of float array | `Marginalize]

type workspace

val new_workspace : T.t -> float array -> workspace

(** [prepare tree p_matrices root_prior leaves] infers ancestral states, given:
- [tree] the phylogenetic tree
- [p_matrices.(i)] is the substitution matrix for the branch leading TO node [i] FROM its parent.
- [root_prior] the prior probability distribution over characters at the root.
- [leaves.(i)] is the observed probability distribution over extant characters at leaf [i]. Since the extant character is usually known with certainty, the distribution usually has 1 for the appropriate entry and 0 otherwise.
@return an abstract value from which various information about ancestral states can be extracted (see below)
*)
val prepare : ?workspace:workspace -> T.t -> P.matrix array -> float array -> leaf array -> intermediate

(** calculate the probability of the leaves under the substitution model *)
val likelihood : intermediate -> float

(** compute the posterior probability distribution over characters at the specified node *)
val node_posterior : intermediate -> int -> float array

(** [branch_posteriors intermediate k] computes the posterior probability of each possible substitution on the specified branch [k]. That is, entry [(i,j)] of the returned matrix is [P(parent(k) = i && k = j | Leaves)] *)
val branch_posteriors : intermediate -> int -> float array array

(** add branch posteriors to the appropriate entries in the given matrix. (useful for summing over many sites) *)
val add_branch_posteriors : ?weight:float -> intermediate -> int -> float array array -> unit
