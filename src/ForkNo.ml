let can_fork = false

module Enum = struct
	let filter_map ?procs f en = Batteries.Enum.filter_map f en
