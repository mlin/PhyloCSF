open List
open Printf
let (|>) x f = f x

type t = {
	tree : T.t;
	qms : Q.Diag.t array;
	pms : P.matrix array;
	prior : float array
}

let make ?prior t qms =
	let qms = if Array.length qms = 1 then Array.make (T.size t - 1) qms.(0) else Array.copy qms
	if Array.length qms <> T.size t - 1 then invalid_arg "CamlPaml.PhyloModel.make"
	for i = 0 to T.root t - 1 do
		if Q.Diag.dim qms.(i) <> Q.Diag.dim qms.(0) || T.branch t i < 0. then invalid_arg "CamlPaml.PhyloModel.make"
	let pms = Array.init (T.size t - 1) (fun br -> Q.Diag.to_Pt qms.(br) (T.branch t br))
	let prior = match prior with	
		| Some pr -> Array.copy pr
		| None ->
			let q = qms.(snd (T.children t (T.root t)))
			if not (Q.Diag.reversible q) then invalid_arg "CamlPaml.PhyloModel.make"
			Q.Diag.equilibrium q
	{ tree = T.copy t; qms = qms; pms = pms; prior = prior }
	
let tree { tree = t } = T.copy t
let q { qms = qms } br = qms.(br)
let p { pms = pms } br = pms.(br) (* Array.map Array.copy pms.(br) *)
let prior { prior = rp } = Array.copy rp

let prepare_lik ?workspace m leaves = PhyloLik.prepare ?workspace m.tree m.pms m.prior leaves

let checksum = 1., 1e-6

let simulate m =
	let branch_choosers = m.pms |> Array.map (fun pm -> (Array.map (Tools.random_chooser ~checksum:checksum) (Gsl_matrix.to_arrays pm)))
	let root_chooser = Tools.random_chooser ~checksum:checksum m.prior
	fun ?root ?a () ->
		let t = m.tree
		let a = match a with
			| Some a -> a
			| None -> Array.make (T.size t) (-1)

		a.(T.root t) <- (match root with
		                   | None -> root_chooser ()
		                   | Some ch when ch >= 0 && ch < Q.Diag.dim m.qms.(0) -> ch
		                   | _ -> invalid_arg "CamlPaml.PhyloModel.simulate (root)")
		for i = T.root t - 1 downto 0 do
			let p = a.(T.parent t i)
			a.(i) <- branch_choosers.(i).(p) ()
		a	

module P14n = struct
	type q_p14n = Expr.t array array

	type model_p14n = {
		q_p14ns : q_p14n array;
		q_scale_p14ns : Expr.t array;
		q_domains : Fit.domain array;

		tree_shape : T.t;
		tree_p14n : Expr.t array;
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

	let instantiate_tree shape exprs domains settings =
		Array.iteri (fun i domain -> if not (Fit.check_domain domain settings.(i)) then invalid_arg ("CamlPaml.P14n.instantiate_tree: domain violation on variable " ^ (string_of_int i))) domains
		
		let tree = T.copy shape
		for br = 0 to T.root tree - 1 do
			T.put_branch tree br (Expr.eval exprs.(br) settings)
		tree
		
	let instantiate_q exprs scale_expr domains settings =
		Array.iteri (fun i domain -> if not (Fit.check_domain domain settings.(i)) then invalid_arg ("CamlPaml.P14n.instantiate_q: domain violation on variable " ^ (string_of_int i))) domains

		let scale = Expr.eval scale_expr settings
		
		if scale <= 0. then
			failwith "CamlPaml.P14n.instantiate_q: Q scale evaluated to a non-positive value"
		
		let qm = exprs |> Array.map (Array.map (fun expr -> Expr.eval expr settings /. scale))
											
		let q = Q.Diag.of_Q qm
		
		q
		
	let instantiate_qs p14ns scale_p14ns domains settings =
		if Array.length p14ns <> Array.length scale_p14ns then invalid_arg ("CamlPaml.P14n.instantiate: different numbers of rate matrix and scale p14ns")
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
						let q = instantiate_q p14ns.(br) scale_p14ns.(br) domains settings
						previous := (p14ns.(br),scale_p14ns.(br),q) :: !previous
						q

	let instantiate ?prior p14n ~q_settings ~tree_settings =
		let qms = instantiate_qs p14n.q_p14ns p14n.q_scale_p14ns p14n.q_domains q_settings
		let tree = instantiate_tree p14n.tree_shape p14n.tree_p14n p14n.tree_domains tree_settings
		{ model = make ?prior tree qms; p14n = p14n; q_settings = Array.copy q_settings; tree_settings = Array.copy tree_settings }
			
	let update ?prior ?q_settings ?tree_settings inst =
		let pms_changed = ref false
		let newq = match q_settings with
			| Some q_settings ->
				pms_changed := true
				instantiate_qs inst.p14n.q_p14ns inst.p14n.q_scale_p14ns inst.p14n.q_domains q_settings
			| None -> inst.model.qms
		let newtree = match tree_settings with
			| Some tree_settings ->
				pms_changed := true
				instantiate_tree inst.p14n.tree_shape inst.p14n.tree_p14n inst.p14n.tree_domains tree_settings
			| None -> inst.model.tree
		let prior = match prior with
			| Some ar -> Array.copy ar
			| None -> inst.model.prior
		let model = if !pms_changed then make ~prior newtree newq else { inst.model with prior = prior }
		{ inst with model = model;
					q_settings = (match q_settings with Some set -> Array.copy set | None -> inst.q_settings);
					tree_settings = (match tree_settings with Some set -> Array.copy set | None -> inst.tree_settings) }
					
	let model { model = model } = model
	let p14n { p14n = p14n } = p14n
	let q_settings { q_settings = q_settings } = Array.copy q_settings
	let tree_settings { tree_settings = tree_settings } = Array.copy tree_settings
	

