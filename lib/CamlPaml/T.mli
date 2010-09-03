(** array-based representation of rooted bifurcating phylogenetic trees

For a tree with [n] leaves, the [2n-1] nodes are indexed such that iterating [0] to [2n-2] performs a botom-up traversal of the tree; indices [0] through [n-1] are the leaves and index [2n-2] is the root. By iterating in order, when you visit node [i] you are guaranteed to have visited its children already. Conversely, by iterating in reverse order, when you visit a node, you have already visited its parent. The data structure can be imported from or exported to a rooted {{:Newick}Newick} tree.*)

type t

val copy : t -> t
val size : t -> int (** total number of nodes in the tree (incl. internal nodes *)

val leaves : t -> int (** number of leaves *)
val is_leaf : t -> int -> bool	(** is node [i] a leaf? (equivalent to [i < (leaves t)]) *)
val root : t -> int		(** index of the root ([= 2*(leaves t)-2 = (size t)-1])*)

(** parent of the given node
@raise Invalid_argument if the root is given *)
val parent : t -> int -> int

(** returns the "left" and "right" child of the internal node
@raise Invalid_argument if the specified node is a leaf *)
val children : t -> int -> int*int

(** return the other child of the parent of the given node.
@raise Invalid_argument if the root is given
*)
val sibling : t -> int -> int

val branch : t -> int -> float (** branch length*)
val put_branch : t -> int -> float -> unit

val label : t -> int -> string
val put_label : t -> int -> string -> unit

(**
@param branch_tol Greatest allowable difference between branch lengths. Default 0, and can be set to infinity to ignore branch lengths.
*)
val congruent : ?tol:float -> labels:bool -> branches:bool -> t -> t -> bool
	
(** The conversion will only succeed for rooted, bifurcating Newick trees.*)
val of_newick : ?default_branch:float -> Newick.t -> t
val to_newick : ?branches:bool -> t -> Newick.t
