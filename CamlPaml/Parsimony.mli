(** inference of ancestral states based on the maximum parsimony heuristic *)

(** output signature of Parsimony functor *)
module type S = sig
	(** arbitrary user-defined character type. perhaps a nucleotide, or an amino acid*)
	type ty
	
	(** Infer ancestral states given the tree and the characters observed at the leaves. When several parsimonious assignments to an ancestral node are possible, the algorithm chooses the assignment matching a plurality of the leaves in the subtree under that node. If there is still a tie, then it is chosen according to the order of the character type.
	@return [n,states] where [n] is the minimum number of substitutions required to explain the observed leaves, and [states] is the inferred character at each node in the tree.
	*)
	val infer : T.t -> ty array -> int*(ty array)

module Make : functor (Ty : Set.OrderedType) -> S with type ty = Ty.t
