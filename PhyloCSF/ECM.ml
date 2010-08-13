(*
Reads empirical codon model (ECM) parameters in the format provided in the supplemental material of:

Kosiol C, Holmes I and Goldman N. An Empirical Codon Model for Protein Sequence Evolution. Mol. Biol. Evol. 2007 24(7):1464-1479; doi:10.1093/molbev/msm064

However, note that Kosiol et al. presented a 61x61 model while PhyloCSF uses a 64x64 model (including stop codons)
*)

open Batteries_uni
open Printf
open CamlPaml
open Expr

module PM = PhyloModel
module Codon = Code.Codon64

let re_sp = Str.regexp " "	
let import_parameters fn_ecm =
	let lines = Array.of_enum (File.lines_of fn_ecm)
	let raw_sij =
		Array.map
			fun line ->
				Array.of_list
					List.map Option.get
						List.map 
							fun s ->
								let s = String.trim s
								if s <> "" then Some (float_of_string s) else None
							Str.split re_sp line
			Array.sub lines 0 (Codon.dim - 1)
	let sij =
		Array.init Codon.dim
			fun i ->
				Array.init Codon.dim
					fun j ->
						if i = 0 && j = 0 then
							0.
						else if i > j then
							raw_sij.(i-1).(j)
						else if i < j then
							raw_sij.(j-1).(i)
						else
							0.

	if String.trim lines.(Codon.dim - 1) <> "" then failwith "ECM.import_parameters"

	let ecm_freqs =
		Array.of_list
			List.map Option.get
				List.map
					fun s ->
						let s = String.trim s
						if s <> "" then Some (float_of_string s) else None
					Str.split re_sp lines.(Codon.dim)
					
	let codons =
		Array.of_list
			List.map Option.get
				List.map
					fun s ->
						let s = String.trim s
						if s <> "" then
							assert (String.length s = 3)
							Some s
						else
							None
					List.flatten (List.map (fun i -> Str.split re_sp lines.(i)) (List.map ((+) Codon.dim) [3; 4; 5; 6]))

	if Array.length codons <> Codon.dim then failwith "ECM.import_parameters: incorrect codon order"
	codons |> Array.iteri (fun i s -> if Codon.index (s.[0],s.[1],s.[2]) <> i then failwith "ECM.import_parameters: incorrect codon order")
	
	sij, ecm_freqs
