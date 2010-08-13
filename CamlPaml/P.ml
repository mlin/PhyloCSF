open Printf

type matrix = float array array

let validate ?(tol=1e-6) p =
	let n = Array.length p
	
	if n < 2 then invalid_arg "CamlPaml.P.validate: trivial matrix"
	
	Array.iteri
		fun i pi ->
			if Array.length pi <> n then invalid_arg "CamlPaml.P.validate: non-square matrix"
			let rowsum = ref 0.
			for j = 0 to n-1 do
				if pi.(j) < 0. || pi.(j) > 1. then invalid_arg (sprintf "CamlPaml.P.validate: P[%d,%d] = %f" i j pi.(j))
				rowsum := !rowsum +. pi.(j)
			if abs_float (!rowsum -. 1.) > tol then
				invalid_arg (sprintf "CamlPaml.P.validate: sum(P[%d,]) = %f != 1" i !rowsum)
		p		
