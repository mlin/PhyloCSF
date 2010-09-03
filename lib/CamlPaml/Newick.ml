open List
open Printf

type branch = float option
type t = Node of (t list)*string*branch

let lbl_bl lbl = function
	| Some bl -> sprintf "%s:%f" lbl bl
	| None -> lbl
		
let rec to_string = function
	| Node ([],lbl,bl) -> lbl_bl lbl bl
	| Node (st,lbl,bl) ->
		sprintf "(%s)%s"
			String.concat "," (map to_string st)
			lbl_bl lbl bl
	
let rec size (Node (st,_,_)) = 1 + fold_left (+) 0 (map size st)

let rec leaves = function
	| Node ([],_,_) -> 1
	| Node (st,_,_) -> fold_left (+) 0 (map leaves st)

let rec remove_branch_lengths (Node (st,lbl,_)) = Node (map remove_branch_lengths st,lbl,None)

let rec list_keep_some f = function
	| [] -> []
	| x :: rest -> match f x with Some x -> x :: list_keep_some f rest | None -> list_keep_some f rest

let maybe_add bl1 bl2 = match (bl1,bl2) with
	| (Some x,Some y) -> Some (x +. y)
	| _ -> None

let rec subtree keep = function
	| (Node ([],lbl,_) as leaf) when lbl = "" || keep lbl -> Some leaf
	| Node (st,lbl,bl) when st <> [] && (lbl = "" || keep lbl) ->
		match list_keep_some (subtree keep) st with
			| [] -> None
			| (Node (sst,snm,sbl)) :: [] -> Some (Node (sst,snm,maybe_add bl sbl))
			| st -> Some (Node (st,lbl,bl))
	| _ -> None

let reorder compare t =
	let rec f = function
		| Node ([],lbl,_) as leaf -> (lbl,leaf)
		| Node (subtrees,lbl,bl) ->
			let f_subtrees = sort (fun (lbl1,_) (lbl2,_) -> compare lbl1 lbl2) (map f subtrees)
			let minlbl = fst (hd f_subtrees)
			(minlbl, Node (map snd f_subtrees,lbl,bl))
	snd (f t)

let rec total_length' = function
	| Node (st,_,Some bl) -> bl +. (fold_left (+.) 0. (map total_length' st))
	| _ -> invalid_arg "CamlPaml.Newick.total_length: unspecified branch length"
	
let total_length ?(count_root=false) t =
	if count_root then
		total_length' t
	else
		match t with
			| Node (st,lbl,_) -> total_length' (Node (st,lbl,Some 0.))
