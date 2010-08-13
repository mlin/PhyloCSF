open List
open Printf
let (|>) x f = f x

type component = {
	tree : T.t;
	qms : Q.Diag.t array;
	pms : P.matrix array;
	root_prior : float array
}

type t = {
	components : component array;
	prior : float array
}

let make_component ?root_prior t qms =
	let qms = if Array.length qms = 1 then Array.make (T.size t - 1) qms.(0) else Array.copy qms
	if Array.length qms <> T.size t - 1 then invalid_arg "CamlPaml.PhyloModel.make_component"
	for i = 0 to T.root t - 1 do
		if T.branch t i < 0. then invalid_arg "CamlPaml.PhyloModel.make_component"
	let pms = Array.init (T.size t - 1) (fun br -> Q.Diag.to_Pt qms.(br) (T.branch t br))
	let root_prior = match root_prior with	
		| Some pr -> pr
		| None ->
			assert (snd (T.children t (T.root t)) = Array.length qms - 1)
			let q = qms.(Array.length qms - 1)
			if not (Q.Diag.reversible q) then invalid_arg "CamlPaml.PhyloModel.make_component"
			Q.Diag.equilibrium q
	{ tree = T.copy t; qms = qms; pms = pms; root_prior = root_prior }
	
let tree { tree = t } = T.copy t
let q { qms = qms } br = qms.(br)
let p { pms = pms } br = pms.(br) (* Array.map Array.copy pms.(br) *)
let root_prior { root_prior = rp } = Array.copy rp

let make components prior =
	if components = [||] || Array.length components <> Array.length prior then invalid_arg "CamlPaml.PhyloModel.make"
	let t0 = components.(0).tree
	let n = Q.Diag.dim components.(0).qms.(0)
	components |> Array.iter
		fun { qms = qms; tree = t } ->
			if T.size t0 <> T.size t then invalid_arg "CamlPaml.PhyloModel.make: incongruent trees"
			if exists (fun q -> Q.Diag.dim q <> n) (Array.to_list qms) then invalid_arg "CamlPaml.PhyloModel.make: inconsistent rate matrix dimensionality"
	{ components = Array.copy components; prior = Array.copy prior }

let component { components = components } i = components.(i)
let component_prior { prior = prior } i = prior.(i)

let components { components = components } = Array.copy components
let components_prior { prior = prior } = Array.copy prior

let infer m leaves =
	let rslts = m.components |> Array.map (fun c -> Infer.prepare c.tree c.pms c.root_prior leaves)
	let liks = rslts |> Array.map Infer.likelihood |> Array.mapi (fun i lik -> lik *. m.prior.(i))
	let tot_lik = Array.fold_left (+.) 0. liks
	let posterior = liks |> Array.map (fun p -> p /. tot_lik)
	tot_lik, posterior, rslts

let random_chooser ?checksum weights =
	let n = Array.length weights
	let cum = Array.copy weights
	for i = 1 to n-1 do
		cum.(i) <- cum.(i-1) +. cum.(i)
	let z = cum.(n-1)
	match checksum with
		| Some (sum,tol) when abs_float (z -. sum) > tol -> invalid_arg (sprintf "random_chooser checksum |%f - %f| > %f" z sum tol)
		| _ -> ()
	(* binary search for i s.t. cum.(i-1) < x <= cum.(i) *)
	let rec bs x lo hi =
		assert (lo <= hi)
		let mid = (lo+hi)/2
		if x > cum.(mid) then
			bs x (mid+1) hi
		else if (mid = 0 || x > cum.(mid-1)) then
			mid
		else
			bs x lo (mid-1)
	fun () -> bs (Random.float z) 0 n

let checksum = 1., 1e-6

let simulate m =
	let branch_choosers = m.components |> Array.map (fun {pms = pms } -> pms |> Array.map (Array.map (random_chooser ~checksum:checksum)))
	let root_choosers = m.components |> Array.map (fun { root_prior = p } -> random_chooser ~checksum:checksum p)
	let component_chooser = random_chooser ~checksum:checksum m.prior
	fun ?a () ->
		let component = component_chooser ()
		let t = m.components.(component).tree
		let a = match a with
			| Some a -> a
			| None -> Array.make (T.size t) (-1)

		a.(T.root t) <- root_choosers.(component) ()
		for i = T.root t - 1 downto 0 do
			let p = a.(T.parent t i)
			a.(i) <- branch_choosers.(component).(i).(p) ()
		a	

module P14n = struct
	type q_p14n = Expr.t array array

	type component_p14n = {
		q_p14ns : q_p14n array;
		q_scale_p14ns : Expr.t array;

		tree_shape : T.t;
		tree_p14n : Expr.t array
	}
	
	type model_p14n = {
		component_p14ns : component_p14n array;
		q_domains : Fit.domain array;
		tree_domains : Fit.domain array
	}

	type instance = {
		model : t;
		p14n : model_p14n;
		q_settings : float array;
		tree_settings : float array
	}
	
	let fill_q_diagonal q =
		let n = Array.length q
		for i = 0 to n-1 do
			if Array.length q.(i) <> n then invalid_arg "CamlPaml.PhyloModel.P14n.fill_q_diagonal"
			let tot = ref (Expr.Val 0.)
			for j = 0 to n-1 do
				if i <> j then
					tot := Expr.Add (q.(i).(j), !tot)
			q.(i).(i) <- Expr.simplify (Expr.Sub (Expr.Val 0., !tot))

	let instantiate_tree shape exprs settings =
		let tree = T.copy shape
		for br = 0 to T.root tree - 1 do
			T.put_branch tree br (Expr.eval exprs.(br) settings)
		tree
		
	let instantiate_q exprs scale_expr settings =
		let scale = Expr.eval scale_expr settings
		
		if scale <= 0. then
			failwith "CamlPaml.P14n.instantiate_q: Q scale evaluated to a non-positive value"
		
		let qm = exprs |> Array.map (Array.map (fun expr -> Expr.eval expr settings /. scale))
											
		let q = Q.Diag.of_Q qm
		
		q
		
	let instantiate_qs p14ns scale_p14ns settings =
		let previous = ref [] (* memoized results from previous branches...I'm assuming the number of independent rate matrix parameterizations to be sublinear in the size of the tree, otherwise this memoization is quadratic... *)
		Array.init (Array.length p14ns)
			fun br ->
				try
					let (_,_,q) =
						find
							function
								| (q_p14n,scale_p14n,q) when q_p14n == p14ns.(br) && scale_p14n == scale_p14ns.(br) -> true
								| _ -> false
							!previous
					q
				with
					| Not_found ->
						let q = instantiate_q p14ns.(br) scale_p14ns.(br) settings
						previous := (p14ns.(br),scale_p14ns.(br),q) :: !previous
						q

	let instantiate_component ?root_prior p14n q_settings tree_settings =
		let qms = instantiate_qs p14n.q_p14ns p14n.q_scale_p14ns q_settings
		let tree = instantiate_tree p14n.tree_shape p14n.tree_p14n tree_settings
		(* TODO cache to_P, although it's not clear it would help *)
		make_component ?root_prior tree qms

	let instantiate ?root_priors p14n ~q_settings ~tree_settings ~components_prior =
		let components =
			p14n.component_p14ns |> Array.mapi
				fun i cp14n ->
					let root_prior = match root_priors with
						| Some ar -> (match ar.(i) with Some rp -> Some rp | None -> None)
						| None -> None
					instantiate_component ?root_prior cp14n q_settings tree_settings
		let model = make components components_prior
		{ p14n = { p14n with component_p14ns = Array.copy p14n.component_p14ns };
			q_settings = Array.copy q_settings; tree_settings = Array.copy tree_settings;
			model = model }
			
	let update ?root_priors ?q_settings ?tree_settings ?components_prior inst =
		let components =
			inst.p14n.component_p14ns |> Array.mapi
				fun i cp14n ->
					let pms_changed = ref false
					let newq = match q_settings with
						| Some q_settings ->
							pms_changed := true
							instantiate_qs cp14n.q_p14ns cp14n.q_scale_p14ns q_settings
						| None -> inst.model.components.(i).qms
					let newtree = match tree_settings with
						| Some tree_settings ->
							pms_changed := true
							instantiate_tree cp14n.tree_shape cp14n.tree_p14n tree_settings
						| None -> inst.model.components.(i).tree
					let root_prior = match root_priors with
						| Some ar -> (match ar.(i) with Some rp -> Array.copy rp | None -> inst.model.components.(i).root_prior)
						| None -> inst.model.components.(i).root_prior
					if !pms_changed then
						make_component ~root_prior:root_prior newtree newq
					else
						{ inst.model.components.(i) with root_prior = root_prior }
		let model = { components = components;
						prior = (match components_prior with Some pr -> pr | None -> inst.model.prior) }
		{ inst with model = model; q_settings = (match q_settings with Some set -> Array.copy set | None -> inst.q_settings);
			tree_settings = (match tree_settings with Some set -> Array.copy set | None -> inst.tree_settings) }
						
					
	
	let model { model = model } = model
	let p14n { p14n = p14n } = { p14n with component_p14ns = Array.copy p14n.component_p14ns }
	let q_settings { q_settings = q_settings } = Array.copy q_settings
	let tree_settings { tree_settings = tree_settings } = Array.copy tree_settings
	

