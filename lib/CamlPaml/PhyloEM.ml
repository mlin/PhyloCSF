open List
open Printf
let (|>) x f = f x

module PM = PhyloModel

type sufficient_statistics = float array array array

let new_sufficient_statistics m =
	let n = Q.Diag.dim (PM.q m 0)
	let t0 = PM.tree m
	Array.init (T.size t0 - 1) (fun _ -> Array.make_matrix n n 0.)

let collect_sufficient_statistics ?workspace m leaves ss =
	let inter = PM.prepare_lik ?workspace m leaves
	let lik = PhyloLik.likelihood inter
	let t = PM.tree m
	for i = 0 to T.root t - 1 do
		PhyloLik.add_branch_posteriors inter i ss.(i)
	lik
	
let clean_sufficient_statistics ?(tol=1e-6) ss =
	ss |> Array.iteri	
		fun br ssbr ->
			ssbr |> Array.iteri
				fun a row ->
					row |> Array.iteri 
						fun b ssab ->
							if ssab < tol then
								row.(b) <- 0.
	
let branch_ell m ss br =
	let n = Q.Diag.dim (PM.q m 0)
	let tot = ref 0.
	
	let ssbr = ss.(br)
	let pbr = PM.p m br

	for k = 0 to n-1 do
		let pbrk = Gsl_matrix.row pbr k
		let ssbrk = ssbr.(k)
		for l = 0 to n-1 do
			let ssbrkl = ssbrk.(l)
			if ssbrkl > 0. then
				tot := !tot +. ssbrkl *. log pbrk.{l}

	!tot
		
let branches_ell m ss =
	let t = PM.tree m
	let nbr = T.size t - 1
	let tot = ref 0.
	for br = 0 to nbr-1 do
		tot := !tot +. branch_ell m ss br
	!tot

let prior_ell m ss =
	let n = Q.Diag.dim (PM.q m 0)
	let tot = ref 0.
	
	let tree = PM.tree m

	let rp = PM.prior m
	let rlc,_ = T.children tree (T.root tree)
	for k = 0 to n-1 do
		let ss_root_k = Array.fold_left (+.) 0. ss.(rlc).(k)
		if ss_root_k > 0. then
			(* root sequence prior *)
			tot := !tot +. ss_root_k *. (log rp.(k))
	!tot


let ell m ss = prior_ell m ss +. branches_ell m ss

(* for consistency it would be nice to do this dbranch. only downside is computation of dQ_dxi would be repeated *)
let d_ell_dQ_dxi inst ss i =
	let m = PM.P14n.model inst
	let p14n = PM.P14n.p14n inst
	let q_settings = PM.P14n.q_settings inst
	
	let tot = ref 0.
	
	let n = Q.Diag.dim (PM.q m 0)
	let tree = PM.tree m
	let nbr = T.size tree - 1

	let previous = ref []

	for br = 0 to nbr-1 do
		let any = ref false
		let dQ_dxi =
			try
				let (_,_,last_dQ_dxi,last_any) =
					find
						function
							| (last_q_p14n,last_scale_p14n,last_dQ_dxi,last_any) when last_q_p14n == p14n.PM.P14n.q_p14ns.(br) && last_scale_p14n == p14n.PM.P14n.q_scale_p14ns.(br) -> true
							| _ -> false
						!previous
				any := last_any
				last_dQ_dxi
			with
				| Not_found  -> 
					let scale = Expr.eval p14n.PM.P14n.q_scale_p14ns.(br) q_settings
					let scale2 = scale *. scale
					let dscale_dxi = Expr.eval (Expr.deriv p14n.PM.P14n.q_scale_p14ns.(br) i) q_settings
					let dQ_dxi =
						p14n.PM.P14n.q_p14ns.(br) |> Array.map
							Array.map
								fun expr ->
									let qij = Expr.eval expr q_settings
									let dqij_dxi = Expr.eval (Expr.deriv expr i) q_settings
									(* Doing the quotient rule numerically. We could do it all symbolically
										with Expr.deriv, but then we'd be recomputing scale & dscale for
										each entry, and those often depend on all the entries in the matrix! *)
									let rslt = (dqij_dxi *. scale -. qij *. dscale_dxi) /. scale2
									if rslt <> 0. then any := true
									rslt
					previous := (p14n.PM.P14n.q_p14ns.(br),p14n.PM.P14n.q_scale_p14ns.(br),dQ_dxi,!any) :: !previous
					dQ_dxi
		if !any then
			let ssbr = ss.(br)
			let t = T.branch tree br
			let dPt_dxi = Q.Diag.dPt_dQ_dx ~q:(PM.q m br) ~t:t ~dQ_dx:dQ_dxi
			let p = PM.p m br

			for a = 0 to n-1 do
				let ssbra = ssbr.(a)
				let dPt_dxi_a = dPt_dxi.(a)
				let pa = Gsl_matrix.row p a
				for b = 0 to n-1 do
					let ssbrab = ssbra.(b)
					if ssbrab > 0. then
						tot := !tot +. ssbrab *. dPt_dxi_a.(b) /. pa.{b}
	!tot

let d_ell_dQ_dx inst ss = Array.init (Array.length (PM.P14n.p14n inst).PM.P14n.q_domains) (d_ell_dQ_dxi inst ss)

(* numerical estimate - for debugging *)
let nd_ell_dQ_dx inst ss =
	let q_settings = PM.P14n.q_settings inst
	let np = Array.length (PM.P14n.p14n inst).PM.P14n.q_domains
	Array.init np
		fun i ->
			let f x =
				let settings = Array.copy q_settings
				settings.(i) <- x
				let xinst =	PM.P14n.update ~q_settings:q_settings inst
				ell (PM.P14n.model xinst) ss
			let { Gsl_fun.res = res } = Gsl_diff.central f q_settings.(i)
			res

let d_ell_dbranch inst br ss =
	let m = PM.P14n.model inst
	let p14n = PM.P14n.p14n inst
	let tree_settings = PM.P14n.tree_settings inst
	let np = Array.length (PM.P14n.p14n inst).PM.P14n.tree_domains

	(* precompute d(ELL)/dt *)
	let tree_dELL_dt =
		let n = Q.Diag.dim (PM.q m 0)
		let tree = PM.tree m
		let q = PM.q m br

		let dPt_dt = Q.Diag.dPt_dt ~q:q ~t:(T.branch tree br)
		let p = PM.p m br
		let tot = ref 0.
		let ssbr = ss.(br)
		for a = 0 to n-1 do
			let ssbra = ssbr.(a)
			let dPt_dt_a = dPt_dt.(a)
			let pa = Gsl_matrix.row p a
			for b = 0 to n-1 do
				let ssbrab = ssbra.(b)
				if ssbrab > 0. then
					tot := !tot +. ssbra.(b) *. dPt_dt_a.(b) /. pa.{b}
		!tot

	Array.init np
		fun i ->
			(* dt/dx on this branch*)
			let dt_dxi = Expr.eval (Expr.deriv p14n.PM.P14n.tree_p14n.(br) i) tree_settings
			(* d(ELL)/dx = d(ELL)/dt dt/dx *)
			tree_dELL_dt *. dt_dxi

let d_ell_dtree inst ss =
	let nbr = T.size (PM.P14n.p14n inst).PM.P14n.tree_shape - 1
	let np = Array.length (PM.P14n.p14n inst).PM.P14n.tree_domains
	
	let rslt = Array.make np 0.
	for br = 0 to nbr-1 do
		let rsltbr = d_ell_dbranch inst br ss
		for i = 0 to np-1 do
			rslt.(i) <- rslt.(i) +. rsltbr.(i)
	rslt
