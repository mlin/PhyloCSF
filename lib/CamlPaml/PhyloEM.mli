(** supporting functions for estimating the parameters of phylogenetic models by
expectation-maximization (EM) *)

(** {1 E-step} *)

(** sufficient statistics for the E-step, composed of estimates of how many times each
character-character substitution occurred on each branch of the tree (summed over all sites) *)
type sufficient_statistics = float array array array

val new_sufficient_statistics : PhyloModel.t -> sufficient_statistics

(** compute the E-step statistics for one site, given the leaves (as would be given to
{{:PhyloModel}PhyloModel.prepare_lik}). Update the given [sufficient_statistics] and return the
probability of the leaves under the model. This should be called for each site to collect the
statistics for the whole alignment. *)
val collect_sufficient_statistics : ?workspace:PhyloLik.workspace -> PhyloModel.t -> PhyloLik.leaf array -> sufficient_statistics -> float

(** any entry in the sufficient statistics less than [tol] is set to zero. useful if e.g. planning
to compress the sufficient statistics for transport between compute nodes *)
val clean_sufficient_statistics : ?tol:float -> sufficient_statistics -> unit

(** {1 M-step} *)

(** expected log-likelihood of the model, where the expectation is over the soft assignment of the
ancestral characters represented by the sufficient statistics of the E-step. (objective function of
the M-step, sometimes called the Q function) *)
val ell : PhyloModel.t -> sufficient_statistics -> float

(** compute the portion of the ell arising from the prior over the ancestral sequence *)
val prior_ell : PhyloModel.t -> sufficient_statistics -> float

(** compute the portion of the ell arising from the substitution process *)
val branches_ell : PhyloModel.t -> sufficient_statistics -> float

(** compute the portion of the ell arising from the substitutions on a specific branch*)
val branch_ell : PhyloModel.t -> sufficient_statistics -> int -> float

(** compute the gradient of the ell with respect to the variables in the model's rate matrix p14n,
based on their symbolic expressions *)
val d_ell_dQ_dx : PhyloModel.P14n.instance -> sufficient_statistics -> float array

(** compute the partial derivative with respect to one specific rate matrix variable only *)
val d_ell_dQ_dxi : PhyloModel.P14n.instance -> sufficient_statistics -> int -> float

(** compute the gradient with respect to the variables in the model's p14n for branch lengths, based
on their symbolic expressions*)
val d_ell_dtree : PhyloModel.P14n.instance -> sufficient_statistics -> float array

(** compute the derivative with respect to one specific branch length variable only *)
val d_ell_dbranch : PhyloModel.P14n.instance -> int -> sufficient_statistics -> float array
