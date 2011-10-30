(** Newick format trees

The Newick format is a parenthesis notation widely used for sharing phylogenetic trees, which can be used to represent rooted and unrooted trees with multifurcating nodes. The recursive data structure here is good for import/export and manipulating the tree topology, but the array-based representation in {{:T}[T]} is much more convenient for probabilistic inference purposes (though restricted to bifurcating trees).

A lexer/parser is of course included, typically used as follows:

[let newick_tree = CamlPaml.NewickParser.parse CamlPaml.NewickLexer.token (Lexing.from_string str)]

@see <http://evolution.genetics.washington.edu/phylip/newicktree.html> information about Newick format
*)

(** optional branch length *)
type branch = float option

(** tree data structure *)
type t = Node of (t list)*string*branch

(** compute the parenthesis representation of the tree (a semicolon is NOT appended) *)
val to_string : t -> string

(** number of nodes in the tree (including internal nodes) *)
val size : t -> int

(** number of leaves in the tree (nodes with no subtrees) *)
val leaves : t -> int

(** remove all branch lengths from the tree *)
val remove_branch_lengths : t -> t

(** [subtree keep tree] removes any nodes that have a label [x] such that [keep x = false]. If an internal node is left with no subtrees as a result of removal of stuff below it, that node is also removed. If an internal node is left with exactly one subtree, then the node is removed, the subtree is connected to the parent of the node, and the branch length of the node is added to the subtree if applicable.

This function cannot be used to remove unnamed nodes except as side effects as described above.
*)
val subtree : (string -> bool) -> t -> t option

(** [reorder cmp t] reorders leaves in the tree to be consistent with the ordering on labels defined by [cmp]. The topology of the tree is not changed, just its arbitrary left-to-right ordering. If the tree is bifurcating, all the leaves will be sorted left-to-right by their labels. In a multifurcating tree, more complicated things can happen (internal nodes are represented by their smallest leaf). The labels, if any, of internal nodes are not used.
*)
val reorder : (string -> string -> int) -> t -> t

(** Compute the total branch length of the tree.
	@param count_root (default false) if true, include the branch length of the topmost node of the tree. If [count_root=false], this branch length is allowed to be unspecified.
	@raise Invalid_argument if the length of any branch is unspecified *)
val total_length : ?count_root:bool -> t -> float
