(** omega (dN/dS) test, similar to using PAML to test for omega<1, but with a few of our own touches *)

open Batteries_uni
open Printf

open CamlPaml
open Expr

module PM = PhyloModel
module Codon = Code.Codon64

(*
rate matrix variables:
0 = kappa, transition/transversion factor
1 = omega, dN/dS ratio
2 = sigma, stop codon frequency factor
3,4,5 freq of A, G. and C relative to freq of T in codon position 1
6,7,8 " in codon position 2
9,10,11 " in codon position 3
*)
let kappa = Var 0
let omega = Var 1
let sigma = Var 2
let pi_exprs =
	let sc (n1,n2,n3) =
		let i1 = Code.DNA.index n1
		let i2 = Code.DNA.index n2
		let i3 = Code.DNA.index n3

		let f1 = Div ((if i1=3 then Val 1.0 else (Var (3+i1))),
			          (Add(Var 3,(Add(Var 4,(Add(Var  5,Val 1.0)))))))
		let f2 = Div ((if i2=3 then Val 1.0 else (Var (6+i2))), 
			          (Add(Var 6,(Add(Var 7,(Add(Var  8,Val 1.0)))))))
		let f3 = Div ((if i3=3 then Val 1.0 else (Var (9+i3))),
			          (Add(Var 9,(Add(Var 10,(Add(Var 11,Val 1.0)))))))
		Mul (f1,Mul(f2,f3))
	let denom = Sub(Val 1.0,
	                Mul(Sub(Val 1.0, sigma),
					    Add(sc ('T','A','A'),Add(sc ('T','A','G'),sc ('T','G','A')))))
	Array.init Codon.dim (fun i -> Div (sc (Codon.of_index i),denom))
	
let q_p14n =
	let sq =
		Array.init Codon.dim
			fun i ->
				let aai = Codon.translate_index i
				let ((i1,i2,i3) as iii) = Codon.of_index i
				Array.init Codon.dim
					fun j ->
						let open List
						let ((j1,j2,j3) as jjj) = Codon.of_index j
						let diffs = filter (fun (a,b) -> a <> b) [(i1,j1); (i2,j2); (i3,j3)]

						if length diffs = 1 then
							let aaj = Codon.translate_index j
							let transition = match hd diffs with
								| 'A','G'
								| 'G','A'
								| 'C','T'
								| 'T','C' -> true
								| _ -> false
								
							let kappa_part = if transition then kappa else Val 1.
							let omega_part =
								if not (Codon.is_stop iii || Codon.is_stop jjj) && aai <> aaj then
									omega
								else
									Val 1. 
							
							simplify (Mul (pi_exprs.(j),Mul(kappa_part,omega_part)))
						else
							Val 0.							
	PM.P14n.fill_q_diagonal sq
	sq

let q_scale =
	let factor = ref (Val 0.)
	for i = 0 to Codon.dim - 1 do
		factor := Sub (!factor,Mul(pi_exprs.(i),q_p14n.(i).(i)))
	simplify !factor

(* Construct branch length parameterizations (all branches scaled by a variable scale factor) *)
let make_tree_p14n tree_shape =
	Array.init (T.size tree_shape - 1) (fun br -> Mul (Var 0,Val (T.branch tree_shape br)))

(* Complete parameterization for the model *)
let make_p14n tree_shape =
	{ PM.P14n.q_p14ns = [| q_p14n |];
	  q_scale_p14ns = [| q_scale |];
	  q_domains = Array.make 12 Fit.NonNeg;
	  tree_shape = T.copy tree_shape;
	  tree_p14n = make_tree_p14n tree_shape;
	  tree_domains = [| Fit.Pos |]; }

let new_instance ?(kappa=1.0) ?(omega=1.0) ?(sigma=1.0) ?(tree_scale=1.0) tree_shape =
	let p14n = make_p14n tree_shape
	let q_settings = Array.concat [ [| kappa; omega; sigma |]; (Array.make 9 1.0) ]
	let tree_settings = [| tree_scale |]
	PM.P14n.instantiate p14n ~q_settings:q_settings ~tree_settings:tree_settings

(* update the rate matrix f3x4 parameters based on the alignment leaves*)
let update_f3x4 inst leaves =
	let counts = Array.init 3 (fun _ -> Array.make 4 1)
	leaves |> Array.iter
		fun lv ->
			lv |> Array.iter
				function
					| `Certain which_codon ->
						let n1,n2,n3 = Codon.of_index which_codon
						let inc a b = counts.(a).(b) <- counts.(a).(b) + 1
						inc 0 (Code.DNA.index n1)
						inc 1 (Code.DNA.index n2)
						inc 2 (Code.DNA.index n3)
					| _ -> ()
	let qs = PM.P14n.q_settings inst
	
	(*
	TODO: Pond et al. Correcting the Bias of Empirical Frequency Parameter Estimators in Codon
	      Models. PLoS One 5:e11230 (2010).
	*)
	
	qs.(3) <- float counts.(0).(0) /. float counts.(0).(3)
	qs.(4) <- float counts.(0).(1) /. float counts.(0).(3)
	qs.(5) <- float counts.(0).(2) /. float counts.(0).(3)

	qs.(6) <- float counts.(1).(0) /. float counts.(1).(3)
	qs.(7) <- float counts.(1).(1) /. float counts.(1).(3)
	qs.(8) <- float counts.(1).(2) /. float counts.(1).(3)
	
	qs.(9) <- float counts.(2).(0) /. float counts.(2).(3)
	qs.(10) <- float counts.(2).(1) /. float counts.(2).(3)
	qs.(11) <- float counts.(2).(2) /. float counts.(2).(3)
	
	PM.P14n.update ~q_settings:qs inst
	
(* Calculate log-probability of the alignment *)
let lpr_leaves inst leaves =
	let workspace = PhyloLik.new_workspace (PM.tree (PM.P14n.model inst)) Codon.dim
	let lik = ref 0.
	leaves |> Array.iter
		fun lvs ->
			let lvslik = PhyloLik.likelihood (PM.prepare_lik ~workspace:workspace (PM.P14n.model inst) lvs)
			lik := !lik +. log lvslik
	!lik

(* a little bit of regularization: half-cauchy prior distributions for rho and kappa *)
let pi = acos (-1.0)
let cauchy_cdf scale x = atan (x /. scale) /. pi +. 0.5
let half_cauchy_lpdf ?(mode=0.0) ~scale x =
	if x < 0.0 || scale <= 0.0 || mode < 0.0 then invalid_arg "half_cauchy_lpdf"
	let numer = 1.0 /. (pi *. scale *. (1.0 +. ((x -. mode) /. scale) ** 2.0))
	let denom = 1.0 -. cauchy_cdf scale (0.0 -. mode)
	assert (denom > 0.0)
	log numer -. log denom

let lpr_rho = half_cauchy_lpdf ~mode:1.0 ~scale:0.5 (* rho = tree_scale (as in manuscript) *)
let lpr_kappa k = log (Gsl_randist.gamma_pdf ~a:7.0 ~b:0.25 (k -. 1.0 +. epsilon_float))

(* find MAP estimates of kappa & rho *)
let kr_map leaves inst =
	let f_rho inst rho =
		let ts = PM.P14n.tree_settings inst
		ts.(0) <- rho
		let inst_rho = PM.P14n.update ~tree_settings:ts inst
		(lpr_rho rho +. lpr_leaves inst_rho leaves), inst_rho
	let f_kappa inst kappa =
		let qs = PM.P14n.q_settings inst
		qs.(0) <- kappa
		let inst_kappa = PM.P14n.update ~q_settings:qs inst
		(lpr_kappa kappa +. lpr_leaves inst_kappa leaves), inst_kappa
	let round inst =
		let _, inst_rho, _ =
			PhyloCSFModel.maximize_lpr
				~init:(PM.P14n.tree_settings inst).(0)
				~lo:0.001
				~hi:10.
				~accuracy:0.01
				f_rho inst
		let _, inst_kappa, lpr =
			PhyloCSFModel.maximize_lpr
				~init:(PM.P14n.q_settings inst).(0)
				~lo:1.0
				~hi:10.
				~accuracy:0.01
				f_kappa inst_rho
		inst_kappa, lpr
	(* 3 rounds, cyclic coordinate ascent *)
	round (fst (round (fst (round inst))))

let sf = sprintf "%.2f"
let db x = sprintf "%.2f" (10. *. x /. log 10.)
	
let score ?(omega_H1=0.2) ?(sigma_H1=0.01) tree_shape leaves =
	(* H0: omega=1, sigma=1, MLE(kappa), MLE(rho) *)
	let inst0, lpr_H0 = kr_map leaves (update_f3x4 (new_instance ~kappa:2.5 tree_shape) leaves)
	let qs0 = PM.P14n.q_settings inst0
	let ts0 = PM.P14n.tree_settings inst0
	
	(* H1: omega=0.2, sigma=0.01, MLE(kappa), MLE(rho) *)
	let inst1, lpr_H1 =	
		let qs = PM.P14n.q_settings inst0
		qs.(1) <- omega_H1
		qs.(2) <- sigma_H1
		kr_map leaves (PM.P14n.update ~q_settings:qs inst0)
	let qs1 = PM.P14n.q_settings inst1
	let ts1 = PM.P14n.tree_settings inst1
	
	let diag = ["L(H0)", (db lpr_H0);
	            "rho_H0", (sf ts0.(0)); "kappa_H0", (sf qs0.(0));
	            "omega_H0", (sf qs0.(1)); "sigma_H0", (sf qs0.(2));
	            "L(H1)", (db lpr_H1);
	            "rho_H1", (sf ts1.(0)); "kappa_H1", (sf qs1.(0));
	            "omega_H1", (sf qs1.(1)); "sigma_H1", (sf qs1.(2)) ]
	
	(10. *. (lpr_H1 -. lpr_H0) /. log 10.), diag
