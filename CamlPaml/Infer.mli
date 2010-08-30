(** Probabilistic inference on phylogenetic trees given the substitution (P) matrices on each branch. These are low-level routines that should usually be used through the higher-level {{:PhyloModel}[PhyloModel]} abstractions. *)

type intermediate

(** Specifying a leaf (extant) character. Usually the extant character is known with certainty; in this case, use [`Certain] with the index of the character, e.g. [`Certain (Code.Codon61.index ('A','T','G'))]. Alternatively, you can specify an arbitrary probability distribution over extant characters. Lastly, you can specify to marginalize a leaf out of the likelihood calculations entirely. *)
type leaf = [`Certain of int | `Distribution of float array	| `Marginalize]

type workspace

(** The calculations use a workspace of [((2 * T.size tree - T.leaves tree) * k)] [float]s where [k] is the alphabet size. As a performance optimization, you can create a workspace with [new_workspace tree k] and use it across multiple calls to [prepare]; otherwise it will allocate by itself. *)
val new_workspace : T.t -> int -> workspace

(** [prepare tree p_matrices root_prior leaves] infers ancestral states, given:
- [tree] the phylogenetic tree
- [p_matrices.(i)] is the substitution matrix for the branch leading TO node [i] FROM its parent.
- [root_prior] the prior probability distribution over characters at the root.
- [leaves] is an array of [leaf]s (see above), the appropriate number for the tree
- [workspace] is an appropriately sized workspace; it will be allocated if not given.
@return an abstract value from which various information about ancestral states can be extracted (see below). If using a shared workspace, be sure to get all the results you need before the next call to [prepare].
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
