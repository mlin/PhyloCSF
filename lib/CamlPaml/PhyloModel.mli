(** unified representation for statistical phylogenetic models, in which the substitution process is a continuous-time Markov process on the phylogeny *)

(** {1 Fully parameterized models}

A 'fully parameterized' model has all of its parameters (each entry of the rate matrix, all branch lengths, etc.) specified as [float]s.
*)

type t

(** [make tree rate_matrices] creates a new model given the tree and rate matrices for each branch (except the branch leading to the root, i.e. [Array.length rate_matrices = T.size tree - 1]).

For the common case of a homogeneous substitution process, there is just one {{:Q} Q matrix} shared throughout the tree. Each entry of [rate_matrices] can therefore just reference the same [Q.Diag.t] value, or, as a special case, you can pass an array of length 1. More generally, even when giving a full array, you would usually want fewer actual independent rate matrices.

@param prior The prior distribution over characters at the root (common ancestor). Defaults to the equilibrium frequencies of [rate_matrices.(Array.length rate_matrices - 1)], the rate matrix on the branch leading to the right child of the root (raises [Invalid_argument] if it's not reversible) *)
val make : ?prior:(float array) -> T.t -> Q.Diag.t array -> t

val tree : t -> T.t
val q : t -> int -> Q.Diag.t (** retrieve the {{:Q}Q matrix} for a specific branch *)
val p : t -> int -> P.matrix (** retrieve the {{:P}P matrix} for a specific branch *)
val prior : t -> float array

(** Simulate evolution according to the model. First, a root character is chosen according to the prior. Then, character substitutions are performed down the tree from the root according to the substitution probabilities determined by the rate matrices and branch lengths.
@param root (optional) a specific character to place at the root; chosen from the prior if not specified.
@param a (optional) a preallocated array of length at least [T.size (tree model)], which will be overwritten and returned, to save memory allocation.
@return the assignment of characters to each node in the tree, in the order of the tree nodes (leaves first)
*)
val simulate : t -> (?root:int -> ?a:(int array) -> unit -> int array)

(** Prepare likelihood calculations for the given configuration of extant characters (leaves).

@param workspace (optional) a preallocated workspace for {{:PhyloLik}PhyloLik.prepare}

@return  {{:PhyloLik}PhyloLik.intermediate} values, from which additional information can be extracted.
*)
val prepare_lik : ?workspace:PhyloLik.workspace -> t -> PhyloLik.leaf array -> PhyloLik.intermediate


(** {1 Symbolic parameterizations}

In symbolic parameterizations (p14ns), parameters are specified using symbolic variables.
*)

module P14n : sig
	type q_p14n = Expr.t array array (** a symbolic expression for each entry of the rate matrix *)

	(** p14n of a model *)
	type model_p14n = {
		q_p14ns : q_p14n array;			(** the rate matrix parameterization for each branch of the tree. As with [make], the array can be length 1 to specify a homogeneous process (the same rate matrix on all branches). *)
		q_scale_p14ns : Expr.t array;	(** a positive scale factor, by which each entry of each rate matrix is divided. Again, can be length 1 for homogeneous processes. *)
		q_domains : Fit.domain array;		(** the domains of all the variables used in the rate matrix p14ns. For example, if [q_domains.(2) = Fit.Pos], then the domain of each occurrence of [(Expr.Var 2)] in the rate matrix p14ns is positive reals. *)

		tree_shape : T.t;				(** the tree topology (branch length settings are ignored) *)
		tree_p14n : Expr.t array;		(** a symbolic expression for each branch length ([tree_exprs.(i)] specifies the length of the branch leading to node [i] from its parent) *)
		tree_domains : Fit.domain array		(** the domains of all the variables used in the branch length p14ns. Note, variables are NOT shared between rate matrix and tree p14ns. For example [(Expr.Var 3)] in [q_p14n] does not refer to the same variable as [(Expr.Var 3)] in [tree_p14n].*)

	}
	
	(** convenience function, sets each [Q.(i).(i)] to minus the sum of all [Q.(i).(j)] for [j <> i]. The resulting expression will have [O(n)] terms, so if you know a more compact way to compute this entry, it would be better to specify it explicitly.*)
	val fill_q_diagonal : q_p14n -> unit
	
	(** in an instance of a p14n we have settings for the variables, thus determining a fully parameterized model *)
	type instance

	(** instantiate the model by giving settings for all the variables in its p14n. *)
	val instantiate : ?prior:(float array) -> model_p14n -> q_settings:(float array) -> tree_settings:(float array) -> instance

	val model : instance -> t
	val p14n : instance -> model_p14n
	val q_settings : instance -> float array
	val tree_settings : instance -> float array
	
	(** return a copy of the instance with a subset of the settings replaced (possibly reusing computed information in the original instance that is not affected by the changed parameters) *)
	val update : ?prior:(float array) -> ?q_settings:(float array) -> ?tree_settings:(float array) -> instance -> instance
	
	
	
		
