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

(* Instantiate a codon model *)
let new_instance ?(tree_scale=1.) s pi tree_shape =
	let p14n = make_p14n s tree_shape
	let q_settings = Array.copy pi
	let tree_settings = [| tree_scale |]
	PM.P14n.instantiate p14n ~prior:pi ~q_settings:q_settings ~tree_settings:tree_settings

(* Calculate log-probability of the alignment *)
let lpr_leaves inst leaves t =
	let ts = PM.P14n.tree_settings inst
	ts.(0) <- t
	let inst = PM.P14n.update ~tree_settings:ts inst
	let workspace = PhyloLik.new_workspace (PM.tree (PM.P14n.model inst)) Codon.dim
	let lik = ref 0.
	leaves |> Array.iter
		fun lvs ->
			let lvslik = PhyloLik.likelihood (PM.prepare_lik ~workspace:workspace (PM.P14n.model inst) lvs)
			lik := !lik +. log lvslik
	!lik, inst

let maximize_lpr ?(init=1.) ?(lo=1e-2) ?(hi=10.) ?(accuracy=0.01) f =
	let good_init = Fit.find_init ~maxtries:250 ~logspace:true ~init:init ~lo:lo ~hi:hi ~f:(fun x -> fst (f x)) ()
	let maximizer = Fit.make_maximizer ~f:(fun x -> fst (f x)) ~lo:lo ~hi:hi ~init:good_init

	let go = ref true
	while !go do
		maximizer#iterate ()
		let x = maximizer#maximum ()
		let lb,ub = maximizer#interval ()
		go := ((ub -. lb) /. x) > accuracy

	let x = maximizer#maximum ()

	let fx, inst = f x
	x, inst, fx

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
let db x = sprintf "%.2f" (10. *. x /. log 10.)

let llr_FixedLik t model leaves =
	let lpr_c, inst_c = lpr_leaves model.coding_model leaves t
	let lpr_n, inst_n = lpr_leaves model.noncoding_model leaves t
	let diag = ["rho",(sf t); "L(C)",(db lpr_c);"L(NC)",(db lpr_n)]
	lpr_c -. lpr_n, diag

let llr_MaxLik ?init model leaves =
	let rho_c, inst_c, lpr_c = maximize_lpr ?init (lpr_leaves model.coding_model leaves)
	let rho_n, inst_n, lpr_n = maximize_lpr ?init (lpr_leaves model.noncoding_model leaves)
	let diag = (match init with Some w0 -> ["rho_0", (sf w0)] | None -> []) @ [ "rho_C",(sf rho_c); "rho_N", (sf rho_n); "L(C)", (db lpr_c); "L(NC)", (db lpr_n)]
	lpr_c -. lpr_n, diag

type strategy =
	| FixedLik
	| MaxLik

let score strategy model leaves =
	let compute_score = match strategy with
		| FixedLik -> llr_FixedLik 1.
		| MaxLik -> llr_MaxLik ~init:1.
	let score, diag = compute_score model leaves
	(10. *. score /. (log 10.)), diag



