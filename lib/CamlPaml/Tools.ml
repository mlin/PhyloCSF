open Printf
open List

let weakly_memoize f =
	let tbl = Hashtbl.create 32
	fun x ->
		try
			match Weak.get (Hashtbl.find tbl x) 0 with
				| None ->
					Hashtbl.remove tbl x
					raise Not_found
				| Some y -> y
		with
			| Not_found ->
				let y = f x
				let box = Weak.create 1
				Weak.set box 0 (Some y)
				Hashtbl.replace tbl x box
				Gc.finalise (fun _ -> Hashtbl.remove tbl x) y
				y
				
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

