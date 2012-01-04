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

(* symmetrize the P matrix according to eq. 12 *)
let symmetrizeP pi p =
	let (m,n) = Gsl_matrix.dims p
	assert (m=n)
	assert (n=Gsl_vector.length pi)

	let pi_sqrt = Gsl_vector.copy pi
	for i = 0 to n-1 do pi_sqrt.{i} <- sqrt pi.{i}
	
	let pi_invsqrt = Gsl_vector.copy pi_sqrt
	for i = 0 to n-1 do pi_invsqrt.{i} <- 1.0 /. pi_sqrt.{i}

	let mm a b =
		let c = Gsl_matrix.create ~init:0.0 n n
		Gsl_blas.(gemm ~ta:NoTrans ~tb:NoTrans ~alpha:1.0 ~a ~b ~beta:0.0 ~c)
		c
	let x = mm (diag pi_sqrt) (mm p (diag pi_invsqrt))
	
	let xT = Gsl_matrix.create n n
	Gsl_matrix.transpose xT x
	
	Gsl_matrix.add x xT
	Gsl_matrix.scale x 0.5
	x
