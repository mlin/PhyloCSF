open List

type t = {
	labels : string array;
	parents : int array;
	children : (int*int) array;
	siblings : int array;
	branches : float array
}

let infer_siblings parents children =
	let n = Array.length parents
	Array.init n
		fun i ->
			if i < (n-1) then
				let p = parents.(i)
				assert (p >= 0 && p < n)
				let (l,r) = children.(p)
				if l == i then
					r
				else if r == i then
					l
				else
					assert false
			else
				(-1)

let copy t = { t with labels = Array.copy t.labels; branches = Array.copy t.branches }
let size t = Array.length t.parents
let root t = (size t)-1
let leaves t = ((size t) + 1)/2
let is_leaf t i = i >= 0 && i < leaves t
let parent t i = if i < 0 || i >= root t then invalid_arg "CamlPaml.T.parent" else t.parents.(i)
let children t i = if i < leaves t || i > root t then invalid_arg "CamlPaml.T.children" else t.children.(i)
let sibling t i = if i < 0 || i >= root t then invalid_arg "CamlPaml.T.sibling" else t.siblings.(i)
let branch t i = t.branches.(i)
let put_branch t i b = t.branches.(i) <- b
let label t i = t.labels.(i)
let put_label t i tag = t.labels.(i) <- tag

exception False
let congruent ?(tol=1e-6) ~labels ~branches t1 t2 =
	if (size t1 <> size t2) then
		false
	else
		try
			for i = 0 to size t1 - 1 do
				if i < root t1 && parent t1 i <> parent t2 i then raise False
				assert (i = root t1 || sibling t1 i = sibling t2 i)
				assert (i < leaves t1 || children t1 i = children t2 i)
				if labels && label t1 i <> label t2 i then raise False
				if branches && abs_float (branch t1 i -. branch t2 i) > tol then raise False
			true
		with
			| False -> false

let of_newick ?(default_branch=nan) nt =
	let bitch () = invalid_arg "CamlPaml.T.of_newick: input is not a rooted, bifurcating tree"
	let n = Newick.size nt
	if n < 3 || (n mod 2 = 0) then bitch ()

	let labels = Array.make n ""
	let parents = Array.make n (-1)
	let children = Array.make n (-1,-1)
	let branches = Array.make n default_branch


	let rec find_leaves = function
		| (Newick.Node ([],_,_)) as leaf -> [leaf]
		| Newick.Node (l :: r :: [],lbl,bl) -> find_leaves l @ find_leaves r
		| _ -> bitch ()
	let leaves = Array.of_list (find_leaves nt)
	let nl = Array.length leaves
	assert (nl = (n+1)/2)

	let which_leaf leaf =
		let j = ref 0
		try
			while not (leaves.(!j) == leaf) do
				incr j
		with
			| Invalid_argument _ -> assert false
		!j

	let fresh_i =
		let i = ref (nl-1)
		fun () ->
			incr i
			!i

	let rec fill = function
		| (Newick.Node ([],lbl,bl)) as leaf ->
			let i = which_leaf leaf
			labels.(i) <- lbl
			branches.(i) <- match bl with Some x -> x | None -> default_branch
			i
		| Newick.Node (l :: r :: [],lbl,bl) ->
			let lc = fill l
			let rc = fill r
			let i = fresh_i ()
			labels.(i) <- lbl
			branches.(i) <- match bl with Some x -> x | None -> default_branch
			parents.(lc) <- i
			parents.(rc) <- i
			children.(i) <- (lc,rc)
			i
		| _ -> bitch ()

	let root = fill nt
	assert (root = n-1)

	{ labels = labels; parents = parents; branches = branches; children = children; siblings = infer_siblings parents children }

let to_newick ?(branches=false) t =
	let branch i =
		if branches then
			match classify_float t.branches.(i) with
				| FP_infinite | FP_nan -> None
				| _ -> Some t.branches.(i)
		else
			None
	let rec cons i =
		if is_leaf t i then
			Newick.Node ([],t.labels.(i),branch i)
		else
			let l,r = t.children.(i)
			Newick.Node ([cons l; cons r],t.labels.(i),branch i)
	cons (root t)
