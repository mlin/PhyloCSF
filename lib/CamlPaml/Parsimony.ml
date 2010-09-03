open List
let (|>) x f = f x

module type S = sig
	type ty
	val infer : T.t -> ty array -> int*(ty array)

module Make = functor (Ty : Set.OrderedType) -> struct
	type ty = Ty.t
	module TySet = Set.Make(Ty)

	let infer t lvs =
		let n = T.size t
		let nl = T.leaves t
		if (Array.length lvs) <> nl then
			invalid_arg "Parsimony.infer: wrong number of leaves in input"
		let sets = Array.make n TySet.empty
		let counts = Array.make n []
		let subs = ref 0
		for i = 0 to n-1 do
			if i < nl then
				sets.(i) <- TySet.singleton lvs.(i)
				counts.(i) <- [(lvs.(i),1)]
			else
				let lc, rc = T.children t i
				assert ((lc <> (-1)) && (rc <> (-1)))
				let inter = TySet.inter sets.(lc) sets.(rc)
				if (TySet.cardinal inter) = 0 then
					sets.(i) <- TySet.union sets.(lc) sets.(rc)
					incr subs
				else
					sets.(i) <- inter
				let cts =
					TySet.elements sets.(i) |> map
						fun x -> (x,(try assoc x counts.(lc) with Not_found -> 0) +
										(try assoc x counts.(rc) with Not_found -> 0))
				counts.(i) <- cts
		let inferred = Array.make n lvs.(0)
		for i = n-1 downto 0 do
			if i < nl then
				inferred.(i) <- lvs.(i)
			else
				if (i <> T.root t) && TySet.mem inferred.(T.parent t i) sets.(i) then
					inferred.(i) <- inferred.(T.parent t i)
				else
					(* assert ((i = T.root t) || (inferred.(T.parent t i) <> "")) *)
					(* TODO prefer synonymous subs? *)
					inferred.(i) <- fst (hd (sort (fun (_,ct1) (_,ct2) -> compare ct2 ct1) counts.(i)))
		!subs, inferred
