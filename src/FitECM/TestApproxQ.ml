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

let veq ?(tol=0.001) v1 v2 =
	let n = Gsl_vector.length v1
	if Gsl_vector.length v2 <> n then false
	else
		let ans = ref true
		for i = 0 to n-1 do
			let d = abs_float (v1.{i} -. v2.{i})
			let mag = min (abs_float v1.{i}) (abs_float v2.{i})
			assert (match classify_float d with FP_infinite | FP_nan -> false | _ -> true)
			assert (match classify_float mag with FP_infinite | FP_nan -> false | _ -> true)
			if (mag > 0. && d /. mag > tol) || (mag = 0. && d > tol) then ans := false
		!ans

let meq ?tol m1 m2 =
	let (r,c) = Gsl_matrix.dims m1
	let (r',c') = Gsl_matrix.dims m2
	if (r <> r' || c <> c') then false
	else
		let ans = ref true
		for i = 0 to r-1 do
			if not (veq ?tol (Gsl_matrix.row m1 i) (Gsl_matrix.row m2 i)) then ans := false
		!ans

let ppm m =
	let (r,c) = Gsl_matrix.dims m
	for i = 0 to r-1 do
		for j = 0 to c-1 do
			printf "\t%.4f" m.{i,j}
		printf "\n"

module Example = struct
	(* Test that the code can step-by-step reproduce a 4x4 example I previously worked out in Mathematica *)

	(* symmetric exchangeabilities *)
	let s =
		[|
			[| 0.0 ; 1.1 ; 0.51; 0.52 |];
			[| 1.1 ; 0.0 ; 0.49; 0.48 |];
			[| 0.51; 0.49; 0.0 ; 0.90 |];
			[| 0.52; 0.48; 0.90; 0.0  |]
		|]

	(* "nucleotide" frequencies *)
	let pi = [| 0.21; 0.31; 0.29; 0.19 |]

	(* reversible rate matrix q_ij = pi_j s_ij *)
	let q =
		[|
			[| -1.19777 ;  0.694982;  0.301431;  0.201361 |];
			[|  0.470794; -0.946276;  0.28961 ;  0.185872 |];
			[|  0.218277;  0.309583; -0.876371;  0.34851  |];
			[|  0.222557;  0.303265;  0.531937; -1.05776  |]
		|]

	(* pairwise substitution matrices we supposedly observe (generated as MatrixExp[0.1*q], MatrixExp[0.25*q], 
	   and MatrixExp[0.5*q]). Our code should recover q from these observations (+pi). *)
	let sms =
		[|
			[|
				[| 0.889106 ; 0.063203 ; 0.0286138; 0.019077  |];
				[| 0.042815 ; 0.911884 ; 0.0274714; 0.0177294 |];
				[| 0.0207203; 0.0294729; 0.917676 ; 0.0321309 |];
				[| 0.0210851; 0.0289269; 0.0490419; 0.900946  |]
			|];

			[|
				[| 0.752013 ; 0.137632 ; 0.0662809; 0.0440743 |];
				[| 0.0932345; 0.801305 ; 0.0641217; 0.0413385 |];
				[| 0.0479964; 0.0685439; 0.812112 ; 0.0713473 |];
				[| 0.0487137; 0.067447 ; 0.108899 ; 0.774941  |]
			|];

			[|
				[| 0.583684 ; 0.221302; 0.117296; 0.0777178 |];
				[| 0.149914 ; 0.662105; 0.114137; 0.0738438 |];
				[| 0.0849388; 0.122008; 0.674872; 0.118181 |];
				[| 0.0858987; 0.120482; 0.180381; 0.613238 |]
			|]
		|]

	(* intermediate result: sum of substitution matrices (~P) *)
	let p =
		[|
			[| 2.2248  ; 0.422137; 0.212191; 0.140869 |];
			[| 0.285964; 2.37529 ; 0.20583 ; 0.132912 |];
			[| 0.153656; 0.220025; 2.40466 ; 0.221659 |];
			[| 0.155697; 0.216856; 0.338322; 2.28912  |]
		|]

	(* intermediate result: symmetrized sum of substitution matrices (eq. 12) *)
	let symP =
		[|
			[| 2.2248  ; 0.347442; 0.180567; 0.148098 |];
			[| 0.347442; 2.37529 ; 0.212809; 0.169772 |];
			[| 0.180567; 0.212809; 2.40466 ; 0.273847 |];
			[| 0.148098; 0.169772; 0.273847; 2.28912  |]
		|]

	(* intermediate result: left and right eigenvectors of q (both row-oriented) *)
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

	(* intermediate result: eigenvalues of q *)
	let ev = [| 0.; 0.610885; 0.848496; 1.0 |]

	let test_symP () =
		let sms' = Array.map Gsl_matrix.of_arrays sms
		Gsl_matrix.add sms'.(0) sms'.(1)
		Gsl_matrix.add sms'.(0) sms'.(2)
		assert (meq (Gsl_matrix.of_arrays symP) (ArvestadBruno1997.symmetrizeP (Gsl_vector.of_array pi) sms'.(0)))

	let test_est_eigenvectors () =
		let lv', rv' = ArvestadBruno1997.est_eigenvectors (Gsl_vector.of_array pi) (Gsl_matrix.of_arrays symP)
		let lv'' = Gsl_matrix.to_arrays lv'

		let d v1 v2 =
			let maxd = ref 0.
			for i = 0 to (Array.length v1) - 1 do
				maxd := max !maxd (abs_float (log (abs_float v1.(i)) -. log (abs_float v2.(i))))
			!maxd

		for i = 0 to 3 do
			let mind = ref infinity
			for j = 0 to 3 do
				(* since the eigenvectors may be returned in arbitrary order, 
				   find the smallest distance between lv.(i) and any row of lv'' *)
				mind := min !mind (d lv.(i) lv''.(j))
			assert (!mind < 0.001)

		let rv'' = Gsl_matrix.to_arrays rv'
		for i = 0 to 3 do
			let mind = ref infinity
			for j = 0 to 3 do
				mind := min !mind (d rv.(i) rv''.(j))
			assert (!mind < 0.001)

	let test_est_eigenvalues () =
		let lv' = Gsl_matrix.of_arrays lv
		let rv' = Gsl_matrix.of_arrays rv
		let ev' = Gsl_vector.to_array (ArvestadBruno1997.est_eigenvalues (Array.map Gsl_matrix.of_arrays sms) lv' rv')
		Array.sort compare ev'

		assert (veq (Gsl_vector.of_array ev) (Gsl_vector.of_array ev'))

	let tests = ["symP" >:: test_symP; "est_eigenvectors" >:: test_est_eigenvectors; "est_eigenvalues" >:: test_est_eigenvalues]

let all_tests = ("FitECM tests" >:::
					[
						"symmetrizeP" >::: symmetrizeP_tests;
						"previously worked-out 4x4 example" >::: Example.tests
					])

run_test_tt_main all_tests
