open Printf

type matrix = Gsl_matrix.matrix

let validate ?(tol=1e-6) p =
	let (n,n') = Gsl_matrix.dims p

	if n <> n' then invalid_arg "CamlPaml.P.validate: non-square matrix"
	if n < 2 then invalid_arg "CamlPaml.P.validate: trivial matrix"

	for i = 0 to n-1 do
		let pi = Gsl_matrix.row p i
		let rowsum = ref 0.
		for j = 0 to n-1 do
			if pi.{j} < 0. || pi.{j} > 1. then invalid_arg (sprintf "CamlPaml.P.validate: P[%d,%d] = %f" i j pi.{j})
			rowsum := !rowsum +. pi.{j}
		if abs_float (!rowsum -. 1.) > tol then
			invalid_arg (sprintf "CamlPaml.P.validate: sum(P[%d,]) = %f != 1" i !rowsum)
