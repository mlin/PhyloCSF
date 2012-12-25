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
	let minimizer = Gsl.Min.make Gsl.Min.BRENT f' ~min:init ~lo:lo ~up:hi

	object
		method maximum () = Gsl.Min.minimum minimizer
		method interval () = Gsl.Min.interval minimizer
		method iterate () =
			Gsl.Min.iterate minimizer
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
				0. -. f (Gsl.Vector.to_array x)
			with
				| any ->
					exn := Some any
					nan
		else
			nan
	let df' ~x ~g =
		if !exn = None then
			try
				let dfv = Gsl.Vector.of_array (df (Gsl.Vector.to_array x))
				Gsl.Vector.scale dfv (-1.) 
				Gsl.Vector.set_zero g
				Gsl.Vector.add g dfv
			with
				| any -> exn := Some any
	let gsl_fdf = {
		Gsl.Fun.multim_f = f';
		Gsl.Fun.multim_df = df';
		Gsl.Fun.multim_fdf = (fun ~x ~g -> let rslt = f' ~x in (if !exn = None then (df' ~x ~g; rslt) else rslt));
	}
	
	
	let minimizer = Gsl.Multimin.Deriv.make Gsl.Multimin.Deriv.VECTOR_BFGS2 n gsl_fdf ~x:(Gsl.Vector.of_array init) ~step:1. ~tol:1e-6
	
	object
		method maximum () = 
			let x = Gsl.Vector.create n
			let g = Gsl.Vector.create n
			let opt = 0. -. Gsl.Multimin.Deriv.minimum ~x:x ~g:g minimizer
			Gsl.Vector.scale g (-1.)

			opt, (Gsl.Vector.to_array x), (Gsl.Vector.to_array g)
		method iterate () =
			Gsl.Multimin.Deriv.iterate minimizer
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




