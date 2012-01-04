open Batteries_uni
open OUnit
open CamlPaml

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

let test_symmetrizeP n () =
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

let all_tests = "symmetrizeP" >::: symmetrizeP_tests

run_test_tt_main all_tests
