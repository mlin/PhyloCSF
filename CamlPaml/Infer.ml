open List
open Printf

module Int = struct
	type t = int
	let compare = compare
module IntSet = Set.Make(Int)

type leaf = [`Certain of int | `Distribution of float array | `Marginalize]

let raw_leaf (k,i) = Gsl_vector.of_array (Array.init k (fun j -> if i = j then 1. else 0.))
let raw_marg k = Gsl_vector.of_array (Array.make k 1.)
let hashcons_raw_leaf = Tools.weakly_memoize raw_leaf
let hashcons_raw_marg = Tools.weakly_memoize raw_marg

let get_leaf k leaf = match leaf with
	| `Certain i -> hashcons_raw_leaf (k,i)
	| `Distribution ar -> Gsl_vector.of_array ar
	| `Marginalize -> hashcons_raw_marg k

type intermediate = {
	tree : T.t;
	pms : Gsl_matrix.matrix array;
	leaves : leaf array;
	alpha : Gsl_matrix.matrix;
	mutable z : float;
	mutable have_alpha : bool;
	beta : Gsl_matrix.matrix;
	mutable have_beta : bool
}

type workspace = Gsl_matrix.matrix
let new_workspace tree dim =
	let rows = 2 * (T.size tree) - (T.leaves tree)
	Gsl_matrix.create rows dim

let empty = [||]
let prepare ?workspace tree pms prior leaves =
	let k = Array.length prior
	let n = T.size tree
	let nl = T.leaves tree
	
	if nl <> Array.length leaves then invalid_arg "CamlPaml.Infer.prepare: length(leaves) != leaves(t)"
	if Array.length pms < n-1 then invalid_arg "CamlPaml.Infer.prepare: not enough P matrices"
	
	let workspace = match workspace with Some x -> x | None -> new_workspace tree k
	let rows,cols = Gsl_matrix.dims workspace
	if rows < (2*n-nl) || cols <> k then invalid_arg "CamlPaml.Infer.prepare: inappropriate workspace dimensions"
	let alpha = Bigarray.Array2.sub_left workspace 0 (n-nl)
	let beta = Bigarray.Array2.sub_left workspace (n-nl) n
	
	for a = 0 to k-1 do beta.{n-1,a} <- prior.(a)
	
	{ tree = tree; pms = pms; leaves = leaves; alpha = alpha; z = nan; have_alpha = false; beta = beta; have_beta = false }

(* Inside algo (aka Felsenstein algo.) *)
let alpha_get x br =
	let k = snd (Gsl_matrix.dims x.alpha)
	let nl = T.leaves x.tree
	if br < nl then get_leaf k x.leaves.(br) else Gsl_matrix.row x.alpha (br-nl)
	
let ensure_alpha x =
	if not x.have_alpha then
		let n = T.size x.tree
		let nl = T.leaves x.tree
		let k = snd (Gsl_matrix.dims x.alpha)
		for i = nl to n-1 do
			let (lc,rc) = T.children x.tree i
			assert (lc >= 0 && rc >= 0)
			assert (lc < i && rc < i)
			let ls = x.pms.(lc)
			let rs = x.pms.(rc)
			let alc = alpha_get x lc
			let arc = alpha_get x rc

			for a = 0 to k-1 do
				let lsa = Gsl_matrix.row ls a
				let rsa = Gsl_matrix.row rs a
				x.alpha.{i-nl,a} <- (Gsl_blas.dot lsa alc) *. (Gsl_blas.dot rsa arc)

		x.z <- Gsl_blas.dot (alpha_get x (n-1)) (Gsl_matrix.row x.beta (n-1))
		x.have_alpha <- true

(* Outside algo *)
let ensure_beta x =
	ensure_alpha x
	if not x.have_beta then
		let k = snd (Gsl_matrix.dims x.beta)
		let inter = Gsl_vector.create k
		let ps_colb = Gsl_vector.create k
		for i = (T.root x.tree)-1 downto 0 do
			let p = T.parent x.tree i
			assert (p > i)
			let s = T.sibling x.tree i
			let ps = x.pms.(i)
			let ss = x.pms.(s)
			let bp = Gsl_matrix.row x.beta p
			let xas = alpha_get x s

			for a = 0 to k-1 do
				inter.{a} <- bp.{a} *. (Gsl_blas.dot (Gsl_matrix.row ss a) xas)

			for b = 0 to k-1 do
				for a = 0 to k-1 do ps_colb.{a} <- ps.{a,b}
				x.beta.{i,b} <- Gsl_blas.dot inter ps_colb
		x.have_beta <- true
		
let likelihood x =
	ensure_alpha x
	x.z

let node_posterior inferred i =
	ensure_beta inferred
	let k = snd (Gsl_matrix.dims inferred.alpha)
	if inferred.z = 0. then
		Array.make k 0.
	else
		let nleaves = T.leaves inferred.tree

		if i < nleaves then
			Gsl_vector.to_array (alpha_get inferred i)
		else
			Array.init k (fun x -> (alpha_get inferred i).{x} *. inferred.beta.{i,x} /. inferred.z)

let add_branch_posteriors ?(weight=1.0) inferred branch ecounts =
	ensure_beta inferred
	let k = snd (Gsl_matrix.dims inferred.alpha)
	
	if Array.length ecounts <> k then invalid_arg "CamlPaml.Infer.add_branch_posteriors: wrong matrix dimension"
	if branch < 0 || branch >= (T.root inferred.tree) then invalid_arg "CamlPaml.Infer.add_branch_posteriors: invalid branch"
	
	let z = inferred.z
	if z > 0. then
		let tree = inferred.tree
		let sm = inferred.pms.(branch)
		let p = T.parent tree branch
		let sib = T.sibling tree branch

		let bp = Gsl_matrix.row inferred.beta p
		let ab = alpha_get inferred branch

		let xas = alpha_get inferred sib
		let sms = inferred.pms.(sib)

		for a = 0 to k-1 do
			let bpa = bp.{a}
			if bpa > 0. then
				let sma = Gsl_matrix.row sm a
				let smsa = Gsl_matrix.row sms a

				let bpa_sibtot = bpa *. (Gsl_blas.dot smsa xas)

				let ecountsa = ecounts.(a)
				if Array.length ecountsa <> k then invalid_arg "CamlPaml.Infer.add_branch_posteriors: wrong matrix dimension"

				for b = 0 to k-1 do
					let pr = bpa_sibtot *. sma.{b} *. ab.{b} /. z
					ecountsa.(b) <- ecountsa.(b) +. weight *. pr

let branch_posteriors inferred branch =
	let k = snd (Gsl_matrix.dims inferred.alpha)
	let a = Array.make_matrix k k 0.
	add_branch_posteriors inferred branch a
	a

