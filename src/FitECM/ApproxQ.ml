open CamlPaml
open Batteries
open Printf
module JSON = Yojson.Basic

let tol = 1e-6
let smooth = 1e-6

(* JSON helpers *)
let vector_of_json = function
	| `List v ->
		Array.of_list
			v |> List.map
				function
					| `Float x -> x
					| `Int x -> float x
					| _ -> invalid_arg "expected numeric vector"
	| _ -> invalid_arg "expected vector"

let matrix_of_json = function
	| `List rows ->
		try
			let m = rows |> List.map vector_of_json |> Array.of_list
			if Array.length m > 0 then
				let nc = Array.length m.(0)
				if not (Array.for_all (fun ar -> Array.length ar = nc) m) then invalid_arg ""
			m
		with
			| Invalid_argument _ -> failwith "expected rectangular numeric matrix"
	| _ -> invalid_arg "expected matrix"

let (@) json key = match json with
	| `Assoc lst -> List.assoc key lst
	| _ -> invalid_arg "expected JSON hash"



(* parse input JSON *)
if Array.length Sys.argv < 2 then
	eprintf "Usage: ApproxECM <file.json>\n"
	exit (-1)

let fn_input = Sys.argv.(1)
let input = JSON.from_file fn_input



(* load and validate input *)
let pi = 
	try
		vector_of_json (input@"CodonFrequencies")
	with
		| Not_found
		| Invalid_argument _ -> failwith ("expected CodonFrequencies to be a numeric vector")

let n = Array.length pi
if n < 2 then failwith "CodonFrequencies vector is too short"

let pisum =
	pi |> Array.fold_left
		fun sum x ->
			if x < 0. then failwith "CodonFrequencies vector has a negative entry"
			sum +. x
		0.
if abs_float (pisum -. 1.0) > tol then failwith "CodonFrequencies vector does not sum to 1"
let pi = Gsl.Vector.of_array pi

let sms =
	try
		let lst = match input@"SubstitutionMatrices" with
			| `List lst -> lst
			| _ -> failwith "expected SubstitutionMatrices to be a list"
		Array.of_list
			List.map
				fun entry ->
					try
						let m = Gsl.Matrix.of_arrays (matrix_of_json (entry@"Matrix"))
						if Gsl.Matrix.dims m <> (n,n) then
							failwith (sprintf "Expected each substitution matrix to be %d-by-%d" n n)

						if smooth > 0. then
							for i = 0 to n-1 do
								for j = 0 to n-1 do
									m.{i,j} <- (m.{i,j} +. smooth) /. (1.0 +. (float n) *. smooth)

						CamlPaml.P.validate m
						let sps = `List [entry@"FirstSpecies"; entry@"SecondSpecies"]
						eprintf "%s\n" (JSON.to_string sps)
						m
					with
						| Not_found -> failwith "expected each member of SubstitutionMatrices to have a Matrix"
				lst
	with
		| Not_found -> failwith "expected SubstitutionMatrices"

(* run algorithm *)

eprintf "Estimating %d-by-%d rate matrix based on %d pairwise substitution matrices...\n" n n (Array.length sms)
flush stderr

let q = ArvestadBruno1997.est_Q pi sms



(* output results *)

for i = 1 to n-1 do
	for j = 0 to i-1 do
		if j > 0 then printf " "
		printf "%f" (q.{i,j} /. pi.{j})
	printf "\n"

pi |> Gsl.Vector.to_array |> Array.to_list |> List.map (sprintf "%f") |> String.concat " " |> printf "\n%s\n\n\n"

if n = 64 then
	(* FIXME: residue vector should be included in input JSON *)
	printf "AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT \
CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT \
GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT \
TTA TTC TTG TTT\n"

(*
let pos = ref 0
while !pos < n do
	let row = min 20 (n - !pos)
	Array.sub codons !pos row |> Array.to_list |> List.map (sprintf "%f") |> String.concat " " |> printf "%s\n"
	pos := !pos + row
*)