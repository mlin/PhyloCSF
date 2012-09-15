(*

Subroutines for the Q matrix estimation strategy described in:

Arvestad L, Bruno WJ (1997). Estimation of Reversible Substitution Matrices from Multiple Pairs of
Sequences. J Mol Evol 45:696-703

*)

open Batteries_uni
open Printf

exception Numerical_error of string

let diag v =
	let n = Gsl_vector.length v
	let m = Gsl_matrix.create ~init:0. n n
	for i = 0 to n-1 do m.{i,i} <- v.{i}
	m

let mm a b =
	let (x,y) = Gsl_matrix.dims a
	let (y',z) = Gsl_matrix.dims b
	if (y <> y') then invalid_arg "mm: matrix dimensions mismatched"
	let c = Gsl_matrix.create ~init:0.0 x z
	Gsl_blas.(gemm ~ta:NoTrans ~tb:NoTrans ~alpha:1.0 ~a ~b ~beta:0.0 ~c)
	c

(* symmetrize the ~P matrix according to eq. 12 *)
let symmetrizeP pi p =
	let (m,n) = Gsl_matrix.dims p
	assert (m=n)
	assert (n=Gsl_vector.length pi)

	let pi_sqrt = Gsl_vector.copy pi
	for i = 0 to n-1 do pi_sqrt.{i} <- sqrt pi.{i}
	(* pi_sqrt is pi^(1/2) as in eq. 10 *)
	
	let pi_invsqrt = Gsl_vector.copy pi_sqrt
	for i = 0 to n-1 do pi_invsqrt.{i} <- 1.0 /. pi_sqrt.{i}
	(* pi_invsqrt is pi^(-1/2) as in eq. 10 *)

	(* now compute eq. 12. TODO this more efficiently than actually constructing & multiplying the diagonal matrices *)
	let x = mm (diag pi_sqrt) (mm p (diag pi_invsqrt))
	
	let xT = Gsl_matrix.create n n
	Gsl_matrix.transpose xT x
	
	Gsl_matrix.add x xT
	Gsl_matrix.scale x 0.5
	x

(* compute the eigenvectors of the symmetrized P matrix,
   and use them to estimate the eigenvectors of P and Q *)
let est_eigenvectors pi symP =
	let (m,n) = Gsl_matrix.dims symP
	assert (m=n)
	assert (n=Gsl_vector.length pi)

	let (_,u) = Gsl_eigen.symmv ~protect:true (`M symP)
	(* columns of u are the eigenvectors of symP *)

	let pi_sqrt = Gsl_vector.copy pi
	for i = 0 to n-1 do pi_sqrt.{i} <- sqrt pi.{i}
	
	let pi_invsqrt = Gsl_vector.copy pi_sqrt
	for i = 0 to n-1 do pi_invsqrt.{i} <- 1.0 /. pi_sqrt.{i}

	let rv = mm (diag pi_invsqrt) u
	Gsl_matrix.transpose_in_place rv
	(* rows of rv are the estimated right eigenvectors of P and Q (eq. 14) *)

	Gsl_matrix.transpose_in_place u
	let lv = mm u (diag pi_sqrt)
	(* rows of lv are estimated left eigenvectors of P and Q (eq. 13) *)

	lv, rv

(* estimate the eigenvalues of Q based on the observations and estimated eigenvectors *)
let est_eigenvalues sms lv rv =
	let (m,n) = Gsl_matrix.dims lv
	assert (m=n)
	assert ((m,n) = Gsl_matrix.dims rv)
	let k = Array.length sms
	assert (k>0)
	sms |> Array.iter (fun sm -> assert ((m,n) = Gsl_matrix.dims sm))

	(* compute distance estimate for each pair of sequences (sm) and nucleotide (eq. 17) *)
	let distances =
		sms |> Array.mapi
			fun which_pair sm ->
				let y = Gsl_vector.create ~init:0. n
				Array.init n
					fun i ->
						Gsl_blas.(gemv NoTrans ~alpha:1.0 ~a:sm ~x:(Gsl_matrix.row rv i) ~beta:0.0 ~y)
						let exp_d_i = (Gsl_blas.(dot (Gsl_matrix.row lv i) y))
						if (exp_d_i <= 0.) then raise (Numerical_error (sprintf "Bad distance estimate for species pair %d, residue %d: log(%e)" which_pair i exp_d_i))
						(* TODO: test cases to provoke negative exp_d_i
						http://www2.gsu.edu/~mkteer/npdmatri.html
						*)
						log exp_d_i

	(* least-squares estimates of the eigenvalues (eq. 18), arbitrarily assigning the last eigenvalue to 1.
		TODO make sure that cannot be the zero eigenvalue. *)
	let denom = distances |> Array.map (fun d -> d.(n-1) *. d.(n-1)) |> Array.fold_left (+.) 0.0
	let denomfpc = classify_float denom
	if (denomfpc = FP_nan || denomfpc = FP_infinite || denom < 1e-6) then
		raise (Numerical_error (sprintf "Bad denominator for eigenvalue estimates: %e" denom))
	let ans =
		Array.init n
			fun i ->
				if i < n-1 then
					distances |> Array.map (fun d -> d.(i) *. d.(n-1)) |> Array.fold_left (+.) 0.0 |> ( *. ) (1.0 /. denom) 
				else 1.0

	(* check all eigenvalues are nonnegative reals *)
	let ans_tot = Array.fold_left (+.) 0. ans

	match classify_float ans_tot with
		| FP_infinite | FP_nan -> raise (Numerical_error "Infinite/undefined eigenvalue estimate")
		| _ -> ()

	(* round slightly negative eigenvalue estimates up to zero *)
	ans |> Array.iteri
		fun i ans_i ->
			if ans_i < 0. then
				if (abs_float ans_i >= 1e-3 *. (float n) *. ans_tot) then
					raise (Numerical_error (sprintf "Eigenvalue estimate %d is too negative: %e" i ans_i))
				ans.(i) <- 0.

	Gsl_vector.of_array ans

(* Normalize rate matrix to unity mean rate of replacement at equilibrium *)
let normalize_Q pi q =
	let (m,n) = Gsl_matrix.dims q
	assert (m=n)
	let qdiag = Gsl_vector.of_array (Array.init m (fun i -> q.{i,i}))

	let z = 1.0 /. (Gsl_blas.dot pi qdiag)
	(* 'sign' of eigenvectors is arbitrary *)
	Gsl_matrix.scale q (if z>0. then 0. -. z else z)

	(* round slightly negative off-diagonal entries to zero *)
	for i = 0 to n-1 do
		let row_tot = ref 0.
		for j = 0 to n-1 do if i <> j then row_tot := !row_tot +. q.{i,j}

		for j = 0 to n-1 do
			if i <> j && q.{i,j} < 0. then
				if abs_float q.{i,j} >= 1e-3 *. (float n) *. !row_tot then
					raise (Numerical_error (sprintf "Rate estimate Q[%d,%d] is too negative: %e" i j q.{i,j}))
				q.{i,j} <- 0.

(* putting it all together *)
let est_Q pi sms =
	let n = Gsl_vector.length pi
	if (n < 2 || Array.length sms = 0) then invalid_arg "ArvestadBruno1997.est_Q"
	sms |> Array.iter (fun sm -> if (Gsl_matrix.dims sm <> (n,n)) then invalid_arg "ArvestadBruno1997.est_Q")

	let p = Gsl_matrix.create ~init:0. n n
	sms |> Array.iter (fun sm -> Gsl_matrix.add p sm)

	let symP = symmetrizeP pi p

	let lv, rv = est_eigenvectors pi symP
	let ev = est_eigenvalues sms lv rv

	(* reconstruct Q according to the spectral thm Q = rv'*diag(ev)*lv
	   http://www.maths.lancs.ac.uk/~gilbert/m306c/node4.html *)
	Gsl_matrix.transpose_in_place rv
	let q = mm rv (mm (diag ev) lv)
	normalize_Q pi q
	q
