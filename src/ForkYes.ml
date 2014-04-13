open Batteries

let can_fork = true

module Enum = struct
	let filter_map ?(procs=1) f en =
		if procs=1 then Batteries.Enum.filter_map f en
		else
			Batteries.Enum.filter_map (fun x -> x)
				Array.enum (ForkWork.map_array ~maxprocs:procs f (Array.of_enum en))
