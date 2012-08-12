(*

Subroutines for the Q matrix estimation strategy described in:

Arvestad L, Bruno WJ (1997). Estimation of Reversible Substitution Matrices from Multiple Pairs of
Sequences. J Mol Evol 45:696-703

*)

open Batteries_uni

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
		sms |> Array.map
			fun sm ->
				let y = Gsl_vector.create ~init:0. m
				Array.init m
					fun i ->
						Gsl_blas.(gemv NoTrans ~alpha:1.0 ~a:sm ~x:(Gsl_matrix.row rv i) ~beta:0.0 ~y)
						log (Gsl_blas.(dot (Gsl_matrix.row lv i) y))

	(* least-squares estimates of the eigenvalues (eq. 18), arbitrarily assigning the last eigenvalue to 1.
		TODO make sure that cannot be the zero eigenvalue. *)
	let denom = distances |> Array.map (fun d -> d.(m-1) *. d.(m-1)) |> Array.fold_left (+.) 0.0
	assert (denom > 10. *. epsilon_float)
	Gsl_vector.of_array
		Array.init m
			fun i ->
				if i < m-1 then
					distances |> Array.map (fun d -> d.(i) *. d.(m-1)) |> Array.fold_left (+.) 0.0 |> ( *. ) (1.0 /. denom) 
				else 1.0

