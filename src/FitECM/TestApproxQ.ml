open Batteries_uni
open OUnit
open CamlPaml
open Printf

Random.init 12345
let tol = 1e-6



exception False
let is_symmetric mat =
	let open Gsl_matrix
	let (m,n) = dims mat
	if m <> n then false
	else
		try
			for i = 0 to m-1 do
				for j = i+1 to m-1 do
					if abs_float (mat.{i,j} -. mat.{j,i}) > tol then raise False
			true
		with False -> false

let random_P n () =
	let p = Gsl_matrix.create n n
	let rowsum i =
		let tot = ref 0.
		for j = 0 to n-1 do tot := !tot +. p.{i,j}
		!tot
	for i = 0 to n-1 do
		for j = 0 to n-1 do p.{i,j} <- float (Random.int 10000 + if i=j then 10000 else 1)
		let z = rowsum i
		for j = 0 to n-1 do p.{i,j} <- p.{i,j} /. z
	
	assert_bool "synthesized P matrix is valid" (try P.validate ~tol p; true with Invalid_argument _ -> false)
	assert_bool "synthesized P matrix is asymmetric" (not (is_symmetric p))

	let pi = Gsl_vector.create n
	for j = 0 to n-1 do pi.{j} <- float (1 + Random.int 10000)
	let z = pi |> Gsl_vector.to_array |> Array.enum |> fold (+.) 0.
	for j = 0 to n-1 do pi.{j} <- pi.{j} /. z

	p, pi

let test_symmetrizeP n () =
	let p, pi = random_P n ()
	let sp = ArvestadBruno1997.symmetrizeP pi p
	assert_bool "symmetrizeP outputs symmetric matrix" (is_symmetric sp)
	
	for i = 0 to n-1 do
		for j = 0  to n-1 do
			assert_bool "symmetrizeP outputs positive matrix" (sp.{i,j} > tol)

let symmetrizeP_tests = [
	"symmetrizeP 4x4" >:: (fun () -> Random.init 0; test_symmetrizeP 4 ());
	"symmetrizeP 4x4 (100x)" >:: (fun () -> Random.init 0; foreach (1 -- 100) (fun _ -> test_symmetrizeP 4 ()));
	"symmetrizeP 64x64 (100x)" >:: (fun () -> Random.init 0; foreach (1 -- 100) (fun _ -> test_symmetrizeP 64 ()))
]



let mm a b =
	let (x,y) = Gsl_matrix.dims a
	let (y',z) = Gsl_matrix.dims b
	if (y <> y') then invalid_arg "mm: matrix dimensions mismatched"
	let c = Gsl_matrix.create ~init:0.0 x z
	Gsl_blas.(gemm ~ta:NoTrans ~tb:NoTrans ~alpha:1.0 ~a ~b ~beta:0.0 ~c)
	c

let d v1 v2 =
	let maxd = ref 0.
	for i = 0 to (Array.length v1) - 1 do
		maxd := max !maxd (abs_float (log (abs_float v1.(i)) -. log (abs_float v2.(i))))
	!maxd

let ppm m =
	let (r,c) = Gsl_matrix.dims m
	for i = 0 to r-1 do
		for j = 0 to c-1 do
			printf "\t%.4f" m.{i,j}
		printf "\n"

let test_est_eigenvectors () =
	(* based on example worked out in Mathematica *)
	let pi = Gsl_vector.of_array [| 0.21; 0.31; 0.29; 0.19 |]
	let p =
		Gsl_matrix.of_arrays
			[|
				[| 2.2248  ; 0.422137; 0.212191; 0.140869 |];
				[| 0.285964; 2.37529 ; 0.20583 ; 0.132912 |];
				[| 0.153656; 0.220025; 2.40466 ; 0.221659 |];
				[| 0.155697; 0.216856; 0.338322; 2.28912  |]
			|]

	let lv =
		[|
			[| -0.21      ; -0.31      ; -0.29      ; -0.19      |];
			[| -0.191046  ; -0.308325  ;  0.300389  ;  0.198981  |];
			[|  0.00767626; -0.00545097; -0.340088  ;  0.337863  |];
			[|  0.359642  ; -0.344683  ; -0.00250571; -0.0124537 |]
		|]

	let rv =
		[|
			[| -1.       ; -1.       ; -1.        ; -1.        |];
			[| -0.909744 ; -0.994596 ;  1.03583   ;  1.04727   |];
			[|  0.0365536; -0.0175838; -1.17272   ;  1.77823   |];
			[|  1.71258  ; -1.11188  ; -0.00864038; -0.0655459 |]
		|]

	let symP = ArvestadBruno1997.symmetrizeP pi p
	let lv', rv' = ArvestadBruno1997.est_eigenvectors pi symP

	let lv' = Gsl_matrix.to_arrays lv'
	for i = 0 to 3 do
		let mind = ref infinity
		for j = 0 to 3 do
			mind := min !mind (d lv.(i) lv'.(j))
		assert (!mind < 0.001)

	Gsl_matrix.transpose_in_place rv'
	let rv' = Gsl_matrix.to_arrays rv'
	for i = 0 to 3 do
		let mind = ref infinity
		for j = 0 to 3 do
			mind := min !mind (d rv.(i) rv'.(j))
		assert (!mind < 0.001)

let all_tests = "FitECM tests" >::: ["symmetrizeP" >::: symmetrizeP_tests; "est_eigenvectors" >:: test_est_eigenvectors]

run_test_tt_main all_tests
