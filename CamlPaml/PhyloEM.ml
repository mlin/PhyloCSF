open List
open Printf
let (|>) x f = f x

module PM = PhyloModel

type sufficient_statistics = (float array)*(float array array array array)

let new_sufficient_statistics m =
	let cs = PM.components m
	let n = Q.Diag.dim (PM.q cs.(0) 0)
	let t0 = PM.tree cs.(0)	
	(Array.make (Array. length cs) 0.), (Array.init (Array.length cs) (fun _ -> Array.init (T.size t0 - 1) (fun _ -> Array.make_matrix n n 0.)))

let collect_sufficient_statistics m leaves (wts,ss) =
	let lik, post, inters = PM.infer m leaves
	let cs = PM.components m
	ss |> Array.iteri
		fun i ssi ->
			let t = PM.tree cs.(i)
			for j = 0 to T.root t - 1 do
				Infer.add_branch_posteriors ~weight:post.(i) inters.(i) j ssi.(j)
			wts.(i) <- wts.(i) +. post.(i)
	lik
	
let clean_sufficient_statistics ?(tol=1e-6) (wts,ss) =
	ss |> Array.iteri
		fun i ssi ->
			ssi |> Array.iteri	
				fun br ssibr ->
					ssibr |> Array.iteri
						fun a row ->
							row |> Array.iteri 
								fun b ssab ->
									if ssab < tol then
										row.(b) <- 0.
	wts |> Array.iteri
		fun i wtsi ->
			if wtsi < tol then
				wts.(i) <- 0.
	
let components_posterior (wts,_) =
	let z = Array.fold_left (+.) 0. wts
	wts |> Array.map (fun x -> x /. z)

let branch_ell ?component m (_,ss) br =
	let cs = PM.components m
	let n = Q.Diag.dim (PM.q cs.(0) 0)
	let tot = ref 0.
	
	let lo, hi = match component with
		| Some i -> i,i
		| None -> 0,(Array.length cs - 1)

	for i = lo to hi do
		let ssibr = ss.(i).(br)
		let pibr = PM.p cs.(i) br

		for k = 0 to n-1 do
			let pibrk = Gsl_matrix.row pibr k
			let ssibrk = ssibr.(k)
			for l = 0 to n-1 do
				let ssibrkl = ssibrk.(l)
				if ssibrkl > 0. then
					tot := !tot +. ssibrkl *. log pibrk.{l}

	!tot
		
let branches_ell ?component m ss =
	let t0 = PM.tree (PM.component m 0)
	let nbr = T.size t0 - 1
	let tot = ref 0.
	for br = 0 to nbr-1 do
		tot := !tot +. branch_ell ?component m ss br
	!tot

let prior_ell ?component m (_,ss) =
	let cs = PM.components m
	let n = Q.Diag.dim (PM.q cs.(0) 0)
	let prior = PM.components_prior m
	let tot = ref 0.
	
	let lo, hi = match component with
		| Some i -> i,i
		| None -> 0,(Array.length cs - 1)

	for i = lo to hi do
		let ssi = ss.(i)
		let tree = PM.tree cs.(i)

		(* terms for priors *)
		let rp = PM.root_prior cs.(i)
		let rlc,_ = T.children tree (T.root tree)
		let ssi_tot = ref 0.
		for k = 0 to n-1 do
			let ssi_root_k = Array.fold_left (+.) 0. ssi.(rlc).(k)
			let rpk = rp.(k)
			if ssi_root_k > 0. then
				(* root sequence prior *)
				tot := !tot +. ssi_root_k *. (log rpk)
			ssi_tot := !ssi_tot +. ssi_root_k	
		if !ssi_tot > 0. then
			(* component prior *)
			tot := !tot +. !ssi_tot *. (log prior.(i))

	!tot

let ell ?component m ss = prior_ell ?component m ss +. branches_ell ?component m ss

(* for consistency it would be nice to do this dbranch. only downside is computation of dQ_dxi would be repeated *)
let d_ell_dQ_dxi inst (_,ss) i =
	let m = PM.P14n.model inst
	let p14n = PM.P14n.p14n inst
	let q_settings = PM.P14n.q_settings inst
	
	let tot = ref 0.
	PM.components m |> Array.iteri
		fun j comp ->
			let n = Q.Diag.dim (PM.q comp 0)
			let tree = PM.tree comp
			let nbr = T.size tree - 1
			let p14nj = p14n.PM.P14n.component_p14ns.(j)
						
			let previous = ref []
			
			for br = 0 to nbr-1 do
				let any = ref false
				let dQ_dxi =
					try
						let (_,_,last_dQ_dxi,last_any) =
							find
								function
									| (last_q_p14n,last_scale_p14n,last_dQ_dxi,last_any) when last_q_p14n == p14nj.PM.P14n.q_p14ns.(br) && last_scale_p14n == p14nj.PM.P14n.q_scale_p14ns.(br) -> true
									| _ -> false
								!previous
						any := last_any
						last_dQ_dxi
					with
						| Not_found  -> 
							let scale = Expr.eval p14nj.PM.P14n.q_scale_p14ns.(br) q_settings
							let scale2 = scale *. scale
							let dscale_dxi = Expr.eval (Expr.deriv p14nj.PM.P14n.q_scale_p14ns.(br) i) q_settings
							let dQ_dxi =
								p14nj.PM.P14n.q_p14ns.(br) |> Array.map
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
							previous := (p14nj.PM.P14n.q_p14ns.(br),p14nj.PM.P14n.q_scale_p14ns.(br),dQ_dxi,!any) :: !previous
							dQ_dxi
				if !any then
					let ssjbr = ss.(j).(br)
					let t = T.branch tree br
					let dPt_dxi = Q.Diag.dPt_dQ_dx ~q:(PM.q comp br) ~t:t ~dQ_dx:dQ_dxi
					let p = PM.p comp br

					for a = 0 to n-1 do
						let ssjbra = ssjbr.(a)
						let dPt_dxi_a = dPt_dxi.(a)
						let pa = Gsl_matrix.row p a
						for b = 0 to n-1 do
							let ssjbrab = ssjbra.(b)
							if ssjbrab > 0. then
								tot := !tot +. ssjbrab *. dPt_dxi_a.(b) /. pa.{b}
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

let d_ell_dbranch inst br (_,ss) =
	let m = PM.P14n.model inst
	let cs = PM.components m
	let p14n = PM.P14n.p14n inst
	let tree_settings = PM.P14n.tree_settings inst
	let np = Array.length (PM.P14n.p14n inst).PM.P14n.tree_domains

	(* precompute d(ELL)/dt for each component *)
	let tree_dELL_dt =
		cs |> Array.mapi
			fun i comp ->
				let n = Q.Diag.dim (PM.q comp 0)
				let tree = PM.tree comp
				let q = PM.q comp br

				let dPt_dt = Q.Diag.dPt_dt ~q:q ~t:(T.branch tree br)
				let p = PM.p comp br
				let tot = ref 0.
				let ssbri = ss.(i).(br)
				for a = 0 to n-1 do
					let ssbria = ssbri.(a)
					let dPt_dt_a = dPt_dt.(a)
					let pa = Gsl_matrix.row p a
					for b = 0 to n-1 do
						let ssbriab = ssbria.(b)
						if ssbriab > 0. then
							tot := !tot +. ssbria.(b) *. dPt_dt_a.(b) /. pa.{b}
				!tot

	Array.init np
		fun i ->
			let tot = ref 0.
			cs |> Array.iteri
				fun j comp ->
					let p14nj = p14n.PM.P14n.component_p14ns.(j)
					
					(* dt/dx on this branch*)
					let dt_dxi = Expr.eval (Expr.deriv p14nj.PM.P14n.tree_p14n.(br) i) tree_settings
					(* d(ELL)/dx = d(ELL)/dt dt/dx *)
					tot := !tot +. tree_dELL_dt.(j) *. dt_dxi
			!tot

let d_ell_dtree inst ss =
	let nbr = T.size (PM.P14n.p14n inst).PM.P14n.component_p14ns.(0).PM.P14n.tree_shape - 1
	let np = Array.length (PM.P14n.p14n inst).PM.P14n.tree_domains
	
	let rslt = Array.make np 0.
	for br = 0 to nbr-1 do
		let rsltbr = d_ell_dbranch inst br ss
		for i = 0 to np-1 do
			rslt.(i) <- rslt.(i) +. rsltbr.(i)
	rslt
