(** substitution probability matrices (P matrices) *)

type matrix = Gsl_matrix.matrix

(** validate that given array is a positive square matrix in which the rows sum to 1
@raise Invalid_argument if not*)
val validate : ?tol:float -> matrix -> unit
