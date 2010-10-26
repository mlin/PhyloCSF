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

let m_of_cm cm = Gsl_matrix.of_arrays (Array.map (Array.map real_of_complex) (Gsl_matrix_complex.to_arrays cm))

let cm_of_m m = Gsl_matrix_complex.of_arrays (Array.map (Array.map (fun re -> { Complex.re = re; im = 0. })) (Gsl_matrix.to_arrays m))

(* stupid wrappers *)

let zgemm ?c a b =
	let m,n = Gsl_matrix_complex.dims a
	let n',p = Gsl_matrix_complex.dims b
	if n <> n' then invalid_arg "CamlPaml.Q.zgemm: incompatible dimensions"
	let c = match c with
		| Some c -> c
		| None -> Gsl_matrix_complex.create m p
	Gsl_blas.Complex.gemm ~ta:Gsl_blas.NoTrans ~tb:Gsl_blas.NoTrans ~alpha:Complex.one ~a:a ~b:b ~beta:Complex.zero ~c:c
	c

let zdiagm ?c d a =
	let m = Gsl_vector_complex.length d
	let (n,p) = Gsl_matrix_complex.dims a
	if m <> n then invalid_arg "CamlPaml.Q.diagm: incompatible dimensions"
	let c = match c with
		| Some c -> c
		| None -> Gsl_matrix_complex.create m p
	
	for i = 0 to m-1 do
		let d_i = d.{i}
		for j = 0 to p-1 do
			c.{i,j} <- Complex.mul a.{i,j} d_i
			
	c

let zinvm m =
	let n,n' = Gsl_matrix_complex.dims m
	if n <> n' then invalid_arg "CamlPaml.Q.zinvm: non-square matrix"
	let p = Gsl_permut.make n
	let lu = Gsl_vectmat.cmat_convert (`CM (Gsl_matrix_complex.copy m))
	ignore (Gsl_linalg.complex_LU_decomp lu p)
	let m' = Gsl_vectmat.cmat_convert (`CM (Gsl_matrix_complex.create n n))
	Gsl_linalg.complex_LU_invert lu p m'
	match m' with
		| `CM x -> x
		| _ -> assert false

module Diag = struct
	(* diagonalized Q = S*L*S'  *)
	type t = {
		q : Gsl_matrix.matrix;			(* Q *)
		s : Gsl_matrix_complex.matrix;	(* S = right eigenvectors (in the columns) *)
		s' : Gsl_matrix_complex.matrix;	(* S' = left eigenvectors (in the rows) *) 
		l : Gsl_vector_complex.vector;	(* diag(L) = eigenvalues *)
		pi : Gsl_vector.vector;
		mutable have_pi : bool;
		
		mutable memoized_to_Pt : (float -> Gsl_matrix.matrix) option;
		
		mutable tol : float;
	}

	let dim q =
		let (n,n') = Gsl_matrix.dims q.q
		assert (n = n')
		n

	let of_Q ?(tol=1e-6) qm =
		let qm = Gsl_matrix.of_arrays qm
		let l, s = Gsl_eigen.nonsymmv ~protect:true (`M qm)
		let s' = zinvm s
		{ q = qm; s = s; s' = s'; l = l;
			pi = Gsl_vector.create (fst (Gsl_matrix.dims qm)); have_pi = false; memoized_to_Pt = None; tol = tol }

	let to_Q q = Gsl_matrix.to_arrays q.q

	let reversible q =
		let rev = ref true
		let n = Gsl_vector_complex.length q.l
		for i = 0 to n-1 do if not (check_real ~tol:q.tol q.l.{i}) then rev := false
		!rev

	(* TODO figure out what this returns for a non-reversible matrix *)
	let equilibrium q =
		if not q.have_pi then
			let min_L = ref infinity
			let min_Lp = ref (-1)
			let n = Gsl_vector_complex.length q.l
			for i = 0 to n-1 do
				let nm_lz = Complex.norm q.l.{i}
				if nm_lz < !min_L then
					min_L := nm_lz
					min_Lp := i
			assert (!min_Lp >= 0)
			if (abs_float !min_L) > q.tol then
				failwith "CamlPaml.Q.equilibrium: smallest-magnitude eigenvalue is unacceptably large; is this a valid rate matrix?"
			let lev = Gsl_matrix_complex.row q.s' !min_Lp
			let mass = ref Complex.zero
			for i = 0 to n-1 do
				mass := Complex.add !mass lev.{i}
			for i = 0 to n-1 do
				q.pi.{i} <- real_of_complex (Complex.div lev.{i} !mass)
			q.have_pi <- true
		Gsl_vector.to_array q.pi

	(** normalize to unity mean rate of replacement 
	let normalize ?tol q =
		let n = Gsl_vector_complex.length q.l
		let pi = equilibrium ?tol q
		let tot = ref 0.
		for i = 0 to n-1 do
			tot := !tot +. (-1.) *. pi.(i) *. q.q.{i,i}
		let nq = Gsl_matrix.copy q.q
		Gsl_matrix.scale q.q (1. /. !tot)
		let nl = Gsl_vector_complex.copy q.l
		for i = 0 to n-1 do
			nl.{i} <- Complex.div nl.{i} { Complex.re = !tot; im = 0. }
		{ q with q = nq; l = nl } *)
		
	let scale q x =
		if x <= 0. then invalid_arg "CamlPaml.Q.scale: nonpositive scale factor"
		let xq = Gsl_matrix.copy q.q
		Gsl_matrix.scale xq x
		let xl = Gsl_vector_complex.copy q.l
		let cx = { Complex.re = x; im = 0. }
		for i = 0 to (Gsl_vector_complex.length xl) - 1 do
			xl.{i} <- Complex.mul xl.{i} cx
		{ q with q = xq; l = xl; memoized_to_Pt = None }
		
	let real_to_Pt q t =
		if t < 0. then invalid_arg "CamlPaml.Q.to_Pt"
		let ct = { Complex.re = t; im = 0. }
		let expLt = Gsl_vector_complex.copy q.l

		for i = 0 to Gsl_vector_complex.length expLt - 1 do
			expLt.{i} <- Complex.exp (Complex.mul ct expLt.{i})	

		let sm = m_of_cm (zgemm q.s (zdiagm expLt q.s'))

		let n,_ = Gsl_matrix.dims sm
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

	let dPt_dt ~q ~t =
		let n,_ = Gsl_matrix.dims q.q
		let (( * ),(+),(/)) = Complex.mul, Complex.add, Complex.div
		let exp= Complex.exp
		let s, l, s' = q.s, q.l, q.s'
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
	let dPt_dQ_dx ~q ~t ~dQ_dx =
		let n,_ = Gsl_matrix.dims q.q
		let dQt_dx = cm_of_m (Gsl_matrix.of_arrays dQ_dx)
		
		if Gsl_matrix_complex.dims dQt_dx <> Gsl_matrix.dims q.q then
			invalid_arg "CamlPaml.Q.dP_dx: wrong dimension of dQ_dx"
		
		let ct = { Complex.re = t; Complex.im = 0. }
		let exp = Complex.exp
		let (( * ),(-),(/)) = Complex.mul, Complex.sub, Complex.div
		
		(* combos 1,3 or 2 (alone) -- 1,3 matches P&S? *)
		Gsl_matrix_complex.scale dQt_dx ct (* 1 *)
		
		let f = Gsl_matrix_complex.create n n
		for i = 0 to pred n do
			for j = 0 to pred n do
				if q.l.{i} = q.l.{j} then
					f.{i,j} <- exp (q.l.{i} * ct) (* * ct *) (* 2 *)
				else
					f.{i,j} <- (exp (q.l.{i} * ct) - exp (q.l.{j} * ct)) / ((q.l.{i} - q.l.{j}) * ct (* 3 *))
		
		(* ehh not being too gentle with the allocator/GC here *)
		Gsl_matrix_complex.mul_elements f (zgemm q.s' (zgemm dQt_dx q.s))
		Gsl_matrix.to_arrays (m_of_cm (zgemm q.s (zgemm f q.s')))

let logm (m:Gsl_matrix.matrix) =
	let n,n' = Gsl_matrix.dims m
	if n <> n' then invalid_arg "CamlPaml.Q.logm: non-square matrix"
	let l, s = Gsl_eigen.nonsymmv ~protect:true (`M m)
	let s' = zinvm s
	for i = 0 to n-1 do
		l.{i} <- Complex.log l.{i}
	
	m_of_cm (zgemm s (zdiagm l s'))

(* Correction of a square matrix with zero row-sums but potentially negative off-diagonal entries into a valid rate matrix. As suggested by Israel, Rosenthal & Wei (2001) *)
let irw2001 rawq =
	let n,n' = Gsl_matrix.dims rawq
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
	Gsl_matrix.to_arrays (irw2001 (logm m))
