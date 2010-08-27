(** unified representation for statistical phylogenetic models, in which the substitution process is a continuous-time Markov process on the phylogeny *)

(** {1 Fully parameterized models}

A 'fully parameterized' model has all of its parameters (each entry of the rate matrix, all branch lengths, etc.) specified as [float]s.
*)

(** A "model component" specifies a complete substitution model with one fixed set of rates. Models may have multiple components, e.g. to incorporate several rate categories in a fixed-effects likelihood model. *)
type component

(** [make_component tree rate_matrices] creates a new model component given the tree and rate matrices for each branch (except the branch leading to the root, i.e. [Array.length rate_matrices = T.size tree - 1]).

For the common case of a homogeneous substitution process, there is just one {{:Q} Q matrix} shared throughout the tree. Each entry of [rate_matrices] can therefore just reference the same [Q.Diag.t] value, or, as a special case, you can pass an array of length 1. More generally, even when giving a full array, you would usually want fewer actual independent rate matrices.

@param root_prior defaults to the equilibrium frequencies of [rate_matrices.(Array.length rate_matrices - 1)], the rate matrix on the branch leading to the right child of the root (raises [Invalid_argument] if it's not reversible) *)
val make_component : ?root_prior:(float array) -> T.t -> Q.Diag.t array -> component

val tree : component -> T.t
val q : component -> int -> Q.Diag.t (** retrieve the {{:Q}Q matrix} for a specific branch *)
val p : component -> int -> P.matrix (** retrieve the {{:P}P matrix} for a specific branch *)
val root_prior : component -> float array

(** A model has one or more components, with a prior probability distribution over the components. *)
type t
val make : component array -> float array -> t

val component : t -> int -> component
val component_prior : t -> int -> float

val components : t -> component array
val components_prior : t -> float array

(** Simulate evolution according to the model. First, a model component is chosen according to the prior over components. Second, a root character is chosen according to the chosen component's prior over root characters. Lastly, character substitutions are performed down the tree from the root according to the substitution probabilities determined by the component's rate matrix and branch lengths.
@param a (optional) a preallocated array of length at least [T.size (tree model)], which will be overwritten and returned, to save memory allocation.
@return the assignment of characters to each node in the tree, in the order of the tree nodes (leaves first)
*)
val simulate : t -> (?a:(int array) -> unit -> int array)

(** Probabilistically infer ancestral characters based on the given configuration of extant characters (leaves). The leaves are an array of length [T.leaves (tree model)]. Each [leaves.(i)] is the probability distribution over extant characters for that leaf (species). Since the extant character is usually known with certainty, this distribution usually has 1 for the corresponding entry and 0 elsewhere.

@return a tuple of the probability of the leaves (likelihood of the model), posterior distribution over the model components, and {{:Infer}intermediate values} for each model component.
*)
val infer : t -> Infer.leaf array -> float*(float array)*(Infer.intermediate array)


(** {1 Symbolic parameterizations}

In symbolic parameterizations (p14ns), parameters are specified using symbolic variables.
*)

module P14n : sig
	type q_p14n = Expr.t array array (** a symbolic expression for each entry of the rate matrix *)

	(** p14n of a model component *)
	type component_p14n = {
		q_p14ns : q_p14n array;			(** the rate matrix parameterization for each branch of the tree. As with [make_component], the array can be length 1 to specify a homogeneous process (the same rate matrix on all branches). *)
		q_scale_p14ns : Expr.t array;	(** a positive scale factor, by which each entry of each rate matrix is divided. Again, can be length 1 for homogeneous processes. *)

		tree_shape : T.t;				(** the tree topology (branch length settings are ignored) *)
		tree_p14n : Expr.t array		(** a symbolic expression for each branch length ([tree_exprs.(i)] specifies the length of the branch leading to node [i] from its parent) *)
	}
	
	(** convenience function, sets each [Q.(i).(i)] to minus the sum of all [Q.(i).(j)] for [j <> i]. The resulting expression will have [O(n)] terms, so if you know a more compact way to compute this entry, it would be better to specify it explicitly.*)
	val fill_q_diagonal : q_p14n -> unit
	
	type model_p14n = {
		component_p14ns : component_p14n array;
		q_domains : Fit.domain array;		(** the domains of all the variables used in the rate matrix p14ns. For example, if [q_domains.(2) = Fit.Pos], then the domain of each occurrence of [(Expr.Var 2)] in the rate matrix p14ns is positive reals. Note, variables are shared among the rate matrix p14ns for all components. For example, all occurrences of [(Expr.Var 2)] in the [q_exprs] of different components refer to the same variable. *)
		tree_domains : Fit.domain array		(** the domains of all the variables used in the branch length p14ns. Like the rate matrix parameters, variables can be shared between the tree p14ns of different components. However, variables are NOT shared between rate matrix and tree p14ns. For example [(Expr.Var 3)] in [q_exprs] does not refer to the same variable as [(Expr.Var 3)] in [tree_exprs].*)
	}
	
	(** in an instance of a p14n we have settings for the variables, thus determining a fully parameterized model *)
	type instance

	(** instantiate the model by giving settings for all the variables in its p14n. *)
	val instantiate : ?root_priors:(float array option array) -> model_p14n -> q_settings:(float array) -> tree_settings:(float array) -> components_prior:(float array) -> instance

	val model : instance -> t
	val p14n : instance -> model_p14n
	val q_settings : instance -> float array
	val tree_settings : instance -> float array
	
	(** return a copy of the instance with a subset of the settings replaced (possibly reusing computed information in the original instance that is not affected by the changed parameters) *)
	val update : ?root_priors:(float array option array) -> ?q_settings:(float array) -> ?tree_settings:(float array) -> ?components_prior:(float array) -> instance -> instance
	
	
	
		
