open List
open Printf

module Int = struct
	type t = int
	let compare = compare
module IntSet = Set.Make(Int)


type intermediate = {
	tree : T.t;
	pms : P.matrix array;
	alpha : float array array;
	mutable z : float;
	mutable have_alpha : bool;
	beta : float array array;
	mutable have_beta : bool
}

let dot (v1:float array) (v2:float array) =
	let l = Array.length v1
	if l <> Array.length v2 then invalid_arg "dot"
	(*SuperFast.dot v1 v2 l*)
	let z = ref 0.
	for i = 0 to l-1 do z := !z +. v1.(i) *. v2.(i)
	!z

let empty = [||]
let prepare tree pms prior leaves =
	let k = Array.length prior
	let n = T.size tree
	let nl = T.leaves tree
	
	if nl <> Array.length leaves then invalid_arg "CamlPaml.Infer.prepare: length(leaves) != leaves(t)"
	if k <> Array.length leaves.(0) then invalid_arg "CamlPaml.infer.prepare: length(root_prior) != length(leaves.(0)) (alphabet size)"
	if Array.length pms < n-1 then invalid_arg "CamlPaml.Infer.prepare: not enough P matrices"

	let alpha =
		Array.init n
			fun i ->
				if i < nl then
					let ll = Array.length leaves.(i)
					if ll = 0 then
						Array.make k 1. (* marginalize *)
					else if ll = k then
						leaves.(i)
					else
						invalid_arg (sprintf "CamlPaml.Infer.prepare: incorrect leaf dimensionality %d (expected %d)" ll k)
				else
					empty
					
	let beta = Array.init n (fun i -> if i = n-1 then prior else empty)					

	{ tree = tree; pms = pms; alpha = alpha; z = nan; have_alpha = false; beta = beta; have_beta = false }

(* Inside algo (aka Felsenstein algo.) *)
let ensure_alpha x =
	if not x.have_alpha then
		let n = T.size x.tree
		let nl = T.leaves x.tree
		let k = Array.length x.alpha.(0)
		for i = nl to n-1 do
			let (lc,rc) = T.children x.tree i
			assert (lc >= 0 && rc >= 0)
			assert (lc < i && rc < i)
			let alc = x.alpha.(lc)
			let arc = x.alpha.(rc)
			let ls = x.pms.(lc)
			let rs = x.pms.(rc)
			if Array.length rs <> k || Array.length ls <> k then
				invalid_arg "CamlPaml.Infer: incorrect substitution matrix dimension"
			let ai =
				Array.init k
					fun a ->
						let lsa = ls.(a)
						let rsa = rs.(a)
						if Array.length rsa <> k || Array.length lsa <> k then
							invalid_arg "CamlPaml.Infer: incorrect substitution matrix dimension"						
						(dot lsa alc) *. (dot rsa arc)
			x.alpha.(i) <- ai
		x.z <- dot x.alpha.(n-1) x.beta.(n-1)
		x.have_alpha <- true

(* Outside algo *)
let ensure_beta x =
	ensure_alpha x
	if not x.have_beta then
		let k = Array.length x.alpha.(0)
		let inter = Array.make k nan
		for i = (T.root x.tree)-1 downto 0 do
			let p = T.parent x.tree i
			assert (p > i)
			let s = T.sibling x.tree i
			let ps = x.pms.(i)
			let ss = x.pms.(s)
			let bp = x.beta.(p)
			let xas = x.alpha.(s)

			for a = 0 to k-1 do
				inter.(a) <- bp.(a) *. (dot ss.(a) xas)

			let bi =
				Array.init k
					fun b ->
						let tot = ref 0.
						for a = 0 to k-1 do
							tot := !tot +. inter.(a) *. ps.(a).(b)
						!tot
			x.beta.(i) <- bi		
		x.have_beta <- true
		
let likelihood x =
	ensure_alpha x
	x.z

let node_posterior inferred i =
	ensure_beta inferred
	let k = Array.length inferred.alpha.(0)
	if inferred.z = 0. then
		Array.make k 0.
	else
		let nleaves = T.leaves inferred.tree

		if i < nleaves then
			Array.copy inferred.alpha.(i)
		else
			Array.init k (fun x -> inferred.alpha.(i).(x) *. inferred.beta.(i).(x) /. inferred.z)

let add_branch_posteriors ?(weight=1.0) inferred branch ecounts =
	ensure_beta inferred
	let k = Array.length inferred.alpha.(0)
	
	if Array.length ecounts <> k then invalid_arg "CamlPaml.Infer.add_branch_posteriors: wrong matrix dimension"
	if branch < 0 || branch >= (T.root inferred.tree) then invalid_arg "CamlPaml.Infer.add_branch_posteriors: invalid branch"
	
	let z = inferred.z
	if z > 0. then
		let tree = inferred.tree
		let sm = inferred.pms.(branch)
		let p = T.parent tree branch
		let sib = T.sibling tree branch

		let bp = inferred.beta.(p)
		let ab = inferred.alpha.(branch)
		assert (Array.length ab = k)

		let xas = inferred.alpha.(sib)
		assert (Array.length xas = k)
		let sms = inferred.pms.(sib)

		for a = 0 to k-1 do
			let bpa = bp.(a)
			if bpa > 0. then
				let sma = sm.(a)
				assert (Array.length sma = k)
				let smsa = sms.(a)
				assert (Array.length smsa = k)

				let bpa_sibtot = bpa *. (dot smsa xas)

				let ecountsa = ecounts.(a)
				if Array.length ecountsa <> k then invalid_arg "CamlPaml.Infer.add_branch_posteriors: wrong matrix dimension"

				for b = 0 to k-1 do
					let pr = bpa_sibtot *. sma.(b) *. ab.(b) /. z
					ecountsa.(b) <- ecountsa.(b) +. weight *. pr

let branch_posteriors inferred branch =
	let k = Array.length inferred.alpha.(0)
	let a = Array.make_matrix k k 0.
	add_branch_posteriors inferred branch a
	a

