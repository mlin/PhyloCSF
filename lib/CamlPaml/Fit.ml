open List

class maximizer f init lo hi =

	let exn = ref None
	let f' x =
		try
			0. -. f x
		with
			| any ->
				exn := Some any
				nan
	let minimizer = Gsl_min.make Gsl_min.BRENT f' ~min:init ~lo:lo ~up:hi

	object
		method maximum () = Gsl_min.minimum minimizer
		method interval () = Gsl_min.interval minimizer
		method iterate () =
			Gsl_min.iterate minimizer
			match !exn with
				| Some ex -> raise ex
				| None -> ()

let make_maximizer ~f ~init ~lo ~hi = (new maximizer f init lo hi)

exception Out_of_range of float
let find_init ?(maxtries=1000) ?(logspace=false) ~f ~init ~lo ~hi () =
	if lo >= hi || (logspace && lo <= 0.) then invalid_arg "CamlPaml.Fit.find_init"
	let width = if logspace then log hi -. log lo else hi -. lo
	let flo = f lo
	let fhi = f hi
	let x = ref init
	let fx = ref (f init)
	let i = ref 0
	Random.init 0
	while !i < maxtries && (!fx <= flo || !fx <= fhi) do
		if logspace then
			x := exp (log lo +. Random.float width)
		else
			x := lo +. Random.float width
		fx := f !x
		incr i
	if !i = maxtries then
		if flo > fhi then
			x := lo
		else
			x := hi
	!x

class multi_maximizer f df init =

	let exn = ref None
	let n = Array.length init
	
	let f' ~x =
		if !exn = None then
			try
				0. -. f (Gsl_vector.to_array x)
			with
				| any ->
					exn := Some any
					nan
		else
			nan
	let df' ~x ~g =
		if !exn = None then
			try
				let dfv = Gsl_vector.of_array (df (Gsl_vector.to_array x))
				Gsl_vector.scale dfv (-1.) 
				Gsl_vector.set_zero g
				Gsl_vector.add g dfv
			with
				| any -> exn := Some any
	let gsl_fdf = {
		Gsl_fun.multim_f = f';
		Gsl_fun.multim_df = df';
		Gsl_fun.multim_fdf = (fun ~x ~g -> let rslt = f' ~x in (if !exn = None then (df' ~x ~g; rslt) else rslt));
	}
	
	
	let minimizer = Gsl_multimin.Deriv.make Gsl_multimin.Deriv.VECTOR_BFGS2 n gsl_fdf ~x:(Gsl_vector.of_array init) ~step:1. ~tol:1e-6
	
	object
		method maximum () = 
			let x = Gsl_vector.create n
			let g = Gsl_vector.create n
			let opt = 0. -. Gsl_multimin.Deriv.minimum ~x:x ~g:g minimizer
			Gsl_vector.scale g (-1.)

			opt, (Gsl_vector.to_array x), (Gsl_vector.to_array g)
		method iterate () =
			Gsl_multimin.Deriv.iterate minimizer
			match !exn with Some ex -> raise ex | None -> ()
			
let make_multi_maximizer ~f ~df ~init = (new multi_maximizer f df init)

type domain = Real | Pos | Neg | NonPos | NonNeg | Probability

let check_domain domain x =
	match domain with
		| Real when x=x -> true
		| Pos when x > 0. -> true
		| Neg when x < 0. -> true
		| NonPos when x <= 0. -> true
		| NonNeg when x >= 0. -> true
		| Probability when x >= 0. && x <= 1. -> true
		| _ -> false


exception False
let enforce_domain ?(rejection_value=neg_infinity) f domains =
	fun x ->
		if (Array.length x) <> Array.length domains then invalid_arg "CamlPaml.Fit.enforce_domains: length mismatch"
		try
			Array.iteri
				fun i domain -> if not (check_domain domain x.(i)) then raise False
				domains
			f x
		with
			| False -> rejection_value




