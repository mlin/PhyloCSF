(** rate matrices for continuous-time Markov models (Q matrices) *)

type matrix = float array array

(** validate that the given array is a square matrix in which all off-diagonal entries are nonnegative and rows sum to 0.
@raise Invalid_argument if not
*)
val validate : ?tol:float -> matrix -> unit

(** diagonalized representation of a rate matrix, from which various useful information can be extracted *)
module Diag : sig
	type t
	
	(** Diagonalize the rate matrix. The rate matrix needs not be reversible (the internal representation uses complex arithmetic). *)
	val of_Q : ?tol:float -> matrix -> t
	val to_Q : t -> matrix
	
	(** return [n] for an [n-by-n] matrix *)
	val dim : t -> int
	
	(** determine if the rate matrix is reversible (<=> the eigenvalues are real) *)
	val reversible : t -> bool
	
	(** for reversible matrices, compute the equilibrium frequencies. The behavior is undefined if the rate matrix is not reversible. *)
	val equilibrium : t -> float array
	
	(** multiply all the rates by a positive scale factor *)
	val scale : t -> float -> t
	
	(** compute the substitution matrix for running the Markov process for time [t] ([=exp(Qt)]). The implementation uses "weak memoization" to cache the results for specific values of [t], memory/GC allowing. *)
	val to_Pt : t -> float -> P.matrix
	
	(** compute the derivative of each entry of [exp(Qt)] with respect to [t]. *)
	val dPt_dt : q:t -> t:float -> float array array
	
	(** compute the derivative of each entry of [exp(Qt)] with respect to some parameter [x]
	@param dQ_dx the derivative of each entry of [Q] with respect to [x]
	*)
	val dPt_dQ_dx : q:t -> t:float -> dQ_dx:(float array array) -> float array array
	
(** Given a {{:P}P matrix} [P], compute a rate matrix [Q] such that [exp(Q)] is approximately equal to [P].

This procedure first computes the matrix logarithm [log(P)], which may or may not be a valid rate matrix. The matrix log is then adjusted to ensure that a valid rate matrix [Q] results, using a method suggested by Israel, Rosenthal & Wei (2001). This is why [exp(Q)] is only {e approximately} equal to [P]. Note that, in any case, the resulting [Q] matrix is {e not} necessarily reversible.
*)
val of_P : ?tol:float -> P.matrix  -> matrix
