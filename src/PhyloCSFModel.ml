open Batteries_uni
open Printf
open CamlPaml
open CamlPaml.Expr

module PM = CamlPaml.PhyloModel
module Codon = CamlPaml.Code.Codon64

(* Construct rate matrix parameterization from symmetric exchangeability matrix s*)
let pivar i = Var i
let make_q_p14n s =
	let sq =
		Array.init Codon.dim
			fun i ->
				let si = s.(i)
				Array.init Codon.dim
					fun j ->
						if i = j then
							Val 0.
						else
							Mul (Val si.(j),pivar j) (* q[i,j] = s[i,j] * pi[j] *)
	PM.P14n.fill_q_diagonal sq
	sq
let qscale sq =
	let factor = ref (Val 0.)
	for i = 0 to Codon.dim - 1 do
		let sqii = sq.(i).(i)
		if sqii <> Val 0. then
			factor := Sub (!factor,Mul(pivar i,sq.(i).(i)))
	simplify !factor

(* Construct branch length parameterizations (all branches scaled by a variable scale factor) *)
let make_tree_p14n tree_shape = Array.init (T.size tree_shape - 1) (fun br -> Mul (Var 0,Val (T.branch tree_shape br)))

(* Complete parameterization for either of the codon models (coding or non-coding *)
let make_p14n s tree_shape =
	let tree_exprs = make_tree_p14n tree_shape
	let q_exprs = make_q_p14n s
	
	{ PM.P14n.q_p14ns = [| q_exprs |]; q_scale_p14ns = [| qscale q_exprs |]; q_domains = [||];
		tree_shape = (T.copy tree_shape); tree_p14n = tree_exprs; tree_domains = [|Fit.Pos|] }

(* dot product... *)
let dot v1 v2 =
	if (Array.length v1 <> Array.length v2) then invalid_arg "dot"
	let tot = ref 0.
	for i = 0 to Array.length v1 - 1 do
		tot := !tot +. v1.(i) *. v2.(i)
	!tot

(* Instantiate a codon model *)
let new_instance ?(tree_scale=1.) s pi tree_shape =
	let p14n = make_p14n s tree_shape
	let q_settings = Array.copy pi
	let tree_settings = [| tree_scale |]
	PM.P14n.instantiate p14n ~prior:pi ~q_settings:q_settings ~tree_settings:tree_settings

type lpr_leaves = {
	lpr_leaves : float;
	elpr_anc : float;
	inst : PM.P14n.instance
}

(* Calculate log-probability of the alignment
   Simultaneously, calculate expected log-probability of the ancestral sequence
*)
let lpr_leaves inst leaves t =
	let ts = PM.P14n.tree_settings inst
	ts.(0) <- t
	let inst = PM.P14n.update ~tree_settings:ts inst
	let workspace = PhyloLik.new_workspace (PM.tree (PM.P14n.model inst)) Codon.dim
	let lpr = ref 0.
	
	let anc_lprior = Array.map log (PM.prior (PM.P14n.model inst))
	let elpr_anc = ref 0.
	leaves |> Array.iter
		fun lvs ->
			let info = PM.prepare_lik ~workspace:workspace (PM.P14n.model inst) lvs
			lpr := !lpr +. log (PhyloLik.likelihood info)
			let pr_root = PhyloLik.node_posterior info (T.root (PM.tree (PM.P14n.model inst)))
			elpr_anc := !elpr_anc +. dot pr_root anc_lprior
	{ lpr_leaves = !lpr; elpr_anc = !elpr_anc; inst }
	
let maximize_lpr ?(init=1.) ?(lo=1e-2) ?(hi=10.) ?(accuracy=0.01) f g =
	let good_init = Fit.find_init ~maxtries:250 ~logspace:true ~init:init ~lo:lo ~hi:hi ~f:(fun x -> g (f x)) ()
	if lo < good_init && good_init < hi then
		let maximizer = Fit.make_maximizer ~f:(fun x -> g (f x)) ~lo:lo ~hi:hi ~init:good_init

		let go = ref true
		while !go do
			maximizer#iterate ()
			let x = maximizer#maximum ()
			let lb,ub = maximizer#interval ()
			go := ((ub -. lb) /. x) > accuracy

		let x = maximizer#maximum ()
		x, f x
	else
		good_init, f good_init

(*  Complete PhyloCSF model, with coding and non-coding ECMs *)
type model = {
	coding_model : PM.P14n.instance;
	noncoding_model : PM.P14n.instance;
}

let make s1 pi1 s2 pi2 tree_shape =
	let coding_model = new_instance s1 pi1 tree_shape
	let noncoding_model = new_instance s2 pi2 tree_shape
	{ coding_model = coding_model; noncoding_model = noncoding_model }

let sf = sprintf "%.2f"
let db x = (10. *. x /. log 10.)
let dbs x = sprintf "%.2f" (db x)

type score = {
  score : float;
  anc_comp_score : float;
  diagnostics: (string*string) list;
}

let llr_FixedLik t model leaves =
	let { lpr_leaves = lpr_c; elpr_anc = elpr_anc_c } = lpr_leaves model.coding_model leaves t
	let { lpr_leaves = lpr_n; elpr_anc = elpr_anc_n} = lpr_leaves model.noncoding_model leaves t
	let diag = ["rho",(sf t); "L(C)",(dbs lpr_c);"L(NC)",(dbs lpr_n)]
	{ score = db (lpr_c -. lpr_n);
	  anc_comp_score = db (elpr_anc_c -. elpr_anc_n);
	  diagnostics = diag }
	
let llr_MaxLik ?init model leaves =
	let rho_c, { lpr_leaves = lpr_c; elpr_anc = elpr_anc_c } = maximize_lpr ?init (lpr_leaves model.coding_model leaves) (fun { lpr_leaves } -> lpr_leaves)
	let rho_n, { lpr_leaves = lpr_n; elpr_anc = elpr_anc_n } = maximize_lpr ?init (lpr_leaves model.noncoding_model leaves) (fun { lpr_leaves } -> lpr_leaves)
	let diag = (match init with Some w0 -> ["rho_0", (sf w0)] | None -> []) @ [ "rho_C",(sf rho_c); "rho_N", (sf rho_n); "L(C)", (dbs lpr_c); "L(NC)", (dbs lpr_n)]
	{ score = db (lpr_c -. lpr_n);
	  anc_comp_score = db (elpr_anc_c -. elpr_anc_n);
	  diagnostics = diag }

type strategy =
	| FixedLik
	| MaxLik

let score strategy model leaves =
	let compute_score = match strategy with
		| FixedLik -> llr_FixedLik 1.
		| MaxLik -> llr_MaxLik ~init:1.
	compute_score model leaves

	
