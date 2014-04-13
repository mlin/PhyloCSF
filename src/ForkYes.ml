open Batteries

let can_fork = true

let map ?(procs=1) f lst =
	if procs=1 || List.length lst = 1 then List.map f lst
	else
		ForkWork.map_list ~maxprocs:procs f lst
