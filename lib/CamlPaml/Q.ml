open Printf

type matrix = float array array

let validate ?(tol=1e-6) q =
	let n = Array.length q
	
	if n < 2 then invalid_arg "CamlPaml.Q.validate: trivial matrix"
	
	Array.iteri
		fun i qi ->
			if Array.length qi <> n then invalid_arg "CamlPaml.Q.validate: non-square matrix"
			let rowsum = ref 0.
			for j = 0 to n-1 do
				if i <> j && qi.(j) < 0. then invalid_arg (sprintf "CamlPaml.Q.validate: Q[%d,%d] = %.2e < 0" i j qi.(j))
				rowsum := !rowsum +. qi.(j)
			if abs_float !rowsum > tol then invalid_arg (sprintf "CamlPaml.Q.validate: sum(Q[%d,]) = %.2e != 0" i !rowsum)
		q

let check_real ?(tol=1e-6) ({ Complex.re = re; Complex.im = im } as z) = im = 0. || (abs_float im) *. 1000. < (Complex.norm z) || (abs_float re < tol && abs_float im < tol)

let real_of_complex ?(tol=1e-6) z =
	if not (check_real ~tol:tol z) then
		failwith (sprintf "CamlPaml.Q.real_of_complex %g+%gi" z.Complex.re z.Complex.im)
	z.Complex.re

let m_of_cm cm = Gsl.Matrix.of_arrays (Array.map (Array.map real_of_complex) (Gsl.Matrix_complex.to_arrays cm))

let cm_of_m m = Gsl.Matrix_complex.of_arrays (Array.map (Array.map (fun re -> { Complex.re = re; im = 0. })) (Gsl.Matrix.to_arrays m))

let v_of_cv cv = Gsl.Vector.of_array (Array.map real_of_complex (Gsl.Vector_complex.to_array cv))

let gemm ?c a b =
	let m,n = Gsl.Matrix.dims a
	let n',p = Gsl.Matrix.dims b
	if n <> n' then invalid_arg "CamlPaml.Q.gemm: incompatible dimensions"
	let c = match c with
		| Some c -> c
		| None -> Gsl.Matrix.create m p
	Gsl.Blas.gemm ~ta:Gsl.Blas.NoTrans ~tb:Gsl.Blas.NoTrans ~alpha:1. ~a:a ~b:b ~beta:0. ~c:c
	c

let zgemm ?c a b =
	let m,n = Gsl.Matrix_complex.dims a
	let n',p = Gsl.Matrix_complex.dims b
	if n <> n' then invalid_arg "CamlPaml.Q.zgemm: incompatible dimensions"
	let c = match c with
		| Some c -> c
		| None -> Gsl.Matrix_complex.create m p
	Gsl.Blas.Complex.gemm ~ta:Gsl.Blas.NoTrans ~tb:Gsl.Blas.NoTrans ~alpha:Complex.one ~a:a ~b:b ~beta:Complex.zero ~c:c
	c

let diagm ?c d a =
	let m = Gsl.Vector.length d
	let (n,p) = Gsl.Matrix.dims a
	if m <> n then invalid_arg "CamlPaml.Q.diagm: incompatible dimensions"
	let c = match c with
		| Some c -> c
		| None -> Gsl.Matrix.create m p
	
	for i = 0 to m-1 do
		let d_i = d.{i}
		for j = 0 to p-1 do
			c.{i,j} <- a.{i,j} *. d_i
			
	c

let zdiagm ?c d a =
	let m = Gsl.Vector_complex.length d
	let (n,p) = Gsl.Matrix_complex.dims a
	if m <> n then invalid_arg "CamlPaml.Q.diagm: incompatible dimensions"
	let c = match c with
		| Some c -> c
		| None -> Gsl.Matrix_complex.create m p
	
	for i = 0 to m-1 do
		let d_i = d.{i}
		for j = 0 to p-1 do
			c.{i,j} <- Complex.mul a.{i,j} d_i
			
	c

let zinvm m =
	let n,n' = Gsl.Matrix_complex.dims m
	if n <> n' then invalid_arg "CamlPaml.Q.zinvm: non-square matrix"
	let p = Gsl.Permut.make n
	let lu = Gsl.Vectmat.cmat_convert (`CM (Gsl.Matrix_complex.copy m))
	ignore (Gsl.Linalg.complex_LU_decomp lu p)
	let m' = Gsl.Vectmat.cmat_convert (`CM (Gsl.Matrix_complex.create n n))
	Gsl.Linalg.complex_LU_invert lu p m'
	match m' with
		| `CM x -> x
		| _ -> assert false

module Diag = struct
	(* diagonalized Q = S*L*S'
	   These are real for reversible models, complex for non-reversible models. *)
	type eig_r = {
		r_s : Gsl.Matrix.matrix;			(* S = right eigenvectors (in the columns) *)
		r_s' : Gsl.Matrix.matrix;			(* S' = left eigenvectors (in the rows) *) 
		r_l : Gsl.Vector.vector 			(* diag(L) = eigenvalues *)
	}
	type eig_nr = {
		nr_s : Gsl.Matrix_complex.matrix;	(* S = right eigenvectors (in the columns) *)
		nr_s' : Gsl.Matrix_complex.matrix;	(* S' = left eigenvectors (in the rows) *) 
		nr_l : Gsl.Vector_complex.vector;	(* diag(L) = eigenvalues *)
	}
	type t = {
		q : Gsl.Matrix.matrix;			(* Q *)
		eig : [`r of eig_r | `nr of eig_nr];
		pi : Gsl.Vector.vector;
		mutable have_pi : bool;
		
		mutable memoized_to_Pt : (float -> Gsl.Matrix.matrix) option;
		
		mutable tol : float;
	}

	let dim q =
		let (n,n') = Gsl.Matrix.dims q.q
		assert (n = n')
		n

	let of_Q ?(tol=1e-6) qm =
		let qm = Gsl.Matrix.of_arrays qm
		let l, s = Gsl.Eigen.nonsymmv ~protect:true (`M qm)
		let s' = zinvm s

		let rev = ref true
		let n = Gsl.Vector_complex.length l
		for i = 0 to n-1 do if not (check_real ~tol l.{i}) then rev := false

		let eig =
			if !rev then
				`r { r_s = m_of_cm s; r_s' = m_of_cm s'; r_l = v_of_cv l }
			else
				`nr { nr_s = s; nr_s' = s'; nr_l = l }

		{ q = qm; eig;
			pi = Gsl.Vector.create (fst (Gsl.Matrix.dims qm)); have_pi = false;
			memoized_to_Pt = None; tol = tol }

	let to_Q q = Gsl.Matrix.to_arrays q.q

	let reversible = function
		| { eig = `r _ } -> true
		| _ -> false

	let equilibrium q =
		if not q.have_pi then
			let eig = match q.eig with
				| `r eig -> eig
				| `nr _ -> failwith "CamlPaml.Q.equilibrium: non-reversible model"
			let n = Gsl.Vector.length eig.r_l
			let min_L = ref infinity
			let min_Lp = ref (-1)
			for i = 0 to n-1 do
				let mag_i = abs_float eig.r_l.{i}
				if mag_i < !min_L then
					min_L := mag_i
					min_Lp := i
			assert (!min_Lp >= 0)
			if (abs_float !min_L) > q.tol then
				failwith (sprintf "CamlPaml.Q.equilibrium: smallest-magnitude eigenvalue %e is unacceptably large; check rate matrix validity or increase tol" !min_L)
			let lev = Gsl.Matrix.row eig.r_s' !min_Lp
			let mass = ref 0.
			for i = 0 to n-1 do
				mass := !mass +. lev.{i}
			for i = 0 to n-1 do
				q.pi.{i} <- lev.{i} /. !mass
			q.have_pi <- true
		Gsl.Vector.to_array q.pi

	(** normalize to unity mean rate of replacement 
	let normalize ?tol q =
		let n = Gsl.Vector_complex.length q.l
		let pi = equilibrium ?tol q
		let tot = ref 0.
		for i = 0 to n-1 do
			tot := !tot +. (-1.) *. pi.(i) *. q.q.{i,i}
		let nq = Gsl.Matrix.copy q.q
		Gsl.Matrix.scale q.q (1. /. !tot)
		let nl = Gsl.Vector_complex.copy q.l
		for i = 0 to n-1 do
			nl.{i} <- Complex.div nl.{i} { Complex.re = !tot; im = 0. }
		{ q with q = nq; l = nl } *)
		
	let scale q x =
		if x <= 0. then invalid_arg "CamlPaml.Q.scale: nonpositive scale factor"
		let xq = Gsl.Matrix.copy q.q
		Gsl.Matrix.scale xq x
		let xeig = match q.eig with
			| `r { r_s; r_s'; r_l } ->
				let xl = Gsl.Vector.copy r_l
				for i = 0 to (Gsl.Vector.length xl) - 1 do
					xl.{i} <- xl.{i} *. x
				`r { r_s; r_s'; r_l = xl }
			| `nr { nr_s; nr_s'; nr_l } ->
				let xl = Gsl.Vector_complex.copy nr_l
				let cx = { Complex.re = x; im = 0. }
				for i = 0 to (Gsl.Vector_complex.length xl) - 1 do
					xl.{i} <- Complex.mul xl.{i} cx
				`nr { nr_s; nr_s'; nr_l = xl }
		{ q with q = xq; eig = xeig; memoized_to_Pt = None }
		
	let real_to_Pt q t =
		if t < 0. then invalid_arg "CamlPaml.Q.to_Pt"
		let sm = match q.eig with
			| `r { r_s; r_s'; r_l } ->
				let expLt = Gsl.Vector.copy r_l
				for i = 0 to Gsl.Vector.length expLt - 1 do
					expLt.{i} <- exp (t *. expLt.{i})
				gemm r_s (diagm expLt r_s')
			| `nr { nr_s; nr_s'; nr_l } ->
				let ct = { Complex.re = t; im = 0. }
				let expLt = Gsl.Vector_complex.copy nr_l
				for i = 0 to Gsl.Vector_complex.length expLt - 1 do
					expLt.{i} <- Complex.exp (Complex.mul ct expLt.{i})	
				m_of_cm (zgemm nr_s (zdiagm expLt nr_s'))

		let n,_ = Gsl.Matrix.dims sm
		for i = 0 to n-1 do
			let tot = ref 0.
			let smii = ref 1.

			for j = 0 to n-1 do
				tot := !tot +. sm.{i,j}

				if sm.{i,j} < 0. then
					if abs_float sm.{i,j} > q.tol then
						failwith (sprintf "CamlPaml.Q.substition_matrix: expm(%.2e*Q)[%d,%d] = %e < 0" t i j sm.{i,j})
					else
						sm.{i,j} <- 0.

				if i <> j then
					smii := !smii -. sm.{i,j}

			if abs_float (!tot -. 1.) > q.tol then
				failwith (sprintf "CamlPaml.Q.substitution matrix: sum(expm(%.2e*Q)[%d,] = %e > 1" t i !tot)
			assert (!smii <= 1. && !smii > 0.)

			sm.{i,i} <- !smii

		sm

	let rec to_Pt q t =
		match q.memoized_to_Pt with
			| Some f -> f t
			| None ->
				q.memoized_to_Pt <- Some (Tools.weakly_memoize (real_to_Pt q))
				to_Pt q t
				
	let to_Pt_gc q = q.memoized_to_Pt <- None

	let dPt_dt ~q ~t = match q.eig with
		| `r _ -> failwith "not implemented"
		| `nr { nr_s; nr_s'; nr_l } ->
			let n,_ = Gsl.Matrix.dims q.q
			let (( * ),(+),(/)) = Complex.mul, Complex.add, Complex.div
			let exp= Complex.exp
			let s, l, s' = nr_s, nr_l, nr_s'
			let ct = { Complex.re = t; Complex.im = 0. }
			Array.init n
				fun a ->
					Array.init n
						fun b ->
							let x = ref Complex.zero
							for c = 0 to n-1 do
								x := !x + (s.{a,c} * l.{c} * (exp (l.{c} * ct)) * s'.{c,b})
							real_of_complex !x
		
	(* TODO in richly parameterized models, dQ_dx is likely to be sparse. *)
	let dPt_dQ_dx ~q ~t ~dQ_dx = match q.eig with
		| `r _ -> failwith "not implemented"
		| `nr { nr_s; nr_s'; nr_l } ->
			let n,_ = Gsl.Matrix.dims q.q
			let dQt_dx = cm_of_m (Gsl.Matrix.of_arrays dQ_dx)
			
			if Gsl.Matrix_complex.dims dQt_dx <> Gsl.Matrix.dims q.q then
				invalid_arg "CamlPaml.Q.dP_dx: wrong dimension of dQ_dx"
			
			let ct = { Complex.re = t; Complex.im = 0. }
			let exp = Complex.exp
			let (( * ),(-),(/)) = Complex.mul, Complex.sub, Complex.div
			
			(* combos 1,3 or 2 (alone) -- 1,3 matches P&S? *)
			Gsl.Matrix_complex.scale dQt_dx ct (* 1 *)
			
			let f = Gsl.Matrix_complex.create n n
			for i = 0 to pred n do
				for j = 0 to pred n do
					if nr_l.{i} = nr_l.{j} then
						f.{i,j} <- exp (nr_l.{i} * ct) (* * ct *) (* 2 *)
					else
						f.{i,j} <- (exp (nr_l.{i} * ct) - exp (nr_l.{j} * ct)) / ((nr_l.{i} - nr_l.{j}) * ct (* 3 *))
			
			(* ehh not being too gentle with the allocator/GC here *)
			Gsl.Matrix_complex.mul_elements f (zgemm nr_s' (zgemm dQt_dx nr_s))
			Gsl.Matrix.to_arrays (m_of_cm (zgemm nr_s (zgemm f nr_s')))

let logm (m:Gsl.Matrix.matrix) =
	let n,n' = Gsl.Matrix.dims m
	if n <> n' then invalid_arg "CamlPaml.Q.logm: non-square matrix"
	let l, s = Gsl.Eigen.nonsymmv ~protect:true (`M m)
	let s' = zinvm s
	for i = 0 to n-1 do
		l.{i} <- Complex.log l.{i}
	
	m_of_cm (zgemm s (zdiagm l s'))

(* Correction of a square matrix with zero row-sums but potentially negative off-diagonal entries into a valid rate matrix. As suggested by Israel, Rosenthal & Wei (2001) *)
let irw2001 rawq =
	let n,n' = Gsl.Matrix.dims rawq
	assert (n = n')
	for i = 0 to n-1 do
		let g_i = ref (abs_float rawq.{i,i})
		let b_i = ref 0.
		for j = 0 to n-1 do
			if i <> j then
				if rawq.{i,j} >= 0. then
					g_i := !g_i +. rawq.{i,j}
				else
					b_i := !b_i -. rawq.{i,j}
		for j = 0 to n-1 do
			let rawqij = rawq.{i,j}
			if i <> j && rawqij < 0. then
				rawq.{i,j} <- 0.
			else if !g_i > 0. then
				rawq.{i,j} <- rawqij -. !b_i *. (abs_float rawqij) /. !g_i
	rawq

let of_P ?tol m =
	P.validate ?tol m
	Gsl.Matrix.to_arrays (irw2001 (logm m))
