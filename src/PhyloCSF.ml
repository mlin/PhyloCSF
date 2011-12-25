(* Wish list: PHYLIP alignment format *)

open Batteries_uni
open OptParse
open Printf
open CamlPaml

Gsl_error.init ()

module SMap = Map.StringMap
module SSet = Set.StringSet
module Codon = CamlPaml.Code.Codon64
type strategy = PhyloCSF of [`MaxLik | `FixedLik] | OmegaTest | Nop
type reading_frame = One | Three | Six
type orf_mode = AsIs | ATGStop | StopStop | StopStop3 | ToFirstStop | FromLastStop | ToOrFromStop

let opt_parser = OptParser.make ~usage:"%prog parameter_set [file1 file2 ...]\ninput will be read from stdin if no filenames are given." ()
let opt ?group ?h ?hide ?s ?short_names ?l ?long_names x = OptParser.add opt_parser ?group ?help:h ?hide ?short_name:s ?long_name:l x; x

let strategy = (opt ~l:"strategy" ~h:"evaluation strategy (default mle)"
                   (Opt.value_option "mle|fixed|omega" (Some (PhyloCSF `MaxLik))
                       (fun s ->
                           match String.lowercase s with
                               | "mle" -> PhyloCSF `MaxLik
                               | "fixed" -> PhyloCSF `FixedLik
                               | "omega" -> OmegaTest
							   | "nop" -> Nop
                               | x -> invalid_arg x)
                       (fun _ s -> sprintf "invalid strategy %s" s)))

let group = OptParser.add_group opt_parser "input interpretation"
let filenames = opt ~group ~l:"files" ~h:"input list(s) of alignment filenames instead of individual alignment(s)" (StdOpt.store_true ())
let remove_ref_gaps = opt ~group ~l:"removeRefGaps" ~h:"automatically remove any alignment columns that are gapped in the reference sequence (nucleotide columns are removed individually; be careful about reading frame). By default, it is an error for the reference sequence to contain gaps" (StdOpt.store_true ())
let allow_ref_gaps = opt ~hide:true ~group ~l:"allowRefGaps" ~h:"allow the reference sequence to contain gaps (each group of three nucleotide columns in the consensus alignment is treated as a codon site; be careful about reading frame)" (StdOpt.store_true ())
let desired_species = opt ~group ~l:"species" ~h:"hint that only this subset of species will be used in any of the alignments; this does not change the calculation mathematically, but can speed it up" (StdOpt.str_option ~metavar:"Species1,Species2,..." ())

let group = OptParser.add_group opt_parser "searching mulitple reading frames and ORFs"
let reading_frame = opt ~group ~l:"frames" ~s:'f' ~h:"how many reading frames to search (default 1)" (Opt.value_option "1|3|6" (Some One) (function "1" -> One | "3" -> Three | "6" -> Six | x -> invalid_arg x) (fun _ s -> sprintf "invalid reading frame %s" s))
let orf_mode = opt ~group ~l:"orf" ~h:"search for ORFs (default AsIs)" (Opt.value_option "AsIs|ATGStop|StopStop|StopStop3|ToFirstStop|FromLastStop|ToOrFromStop" (Some AsIs)  (fun s -> match String.lowercase s with "asis" -> AsIs | "atgstop" -> ATGStop | "stopstop" -> StopStop | "stopstop3" -> StopStop3 | "tofirststop" -> ToFirstStop | "fromlaststop" -> FromLastStop | "toorfromstop" -> ToOrFromStop| x -> invalid_arg x) (fun _ s -> sprintf "invalid ORF search mode %s"s))
let min_codons = opt ~group ~l:"minCodons" ~h:"minimum ORF length for searching over ORFs (default 25 codons)" (StdOpt.int_option ~default:25 ())
let print_orfs = opt ~group ~l:"allScores" ~h:"report scores of all regions evaluated, not just the max" (StdOpt.store_true ())

let group = OptParser.add_group opt_parser "output control"
let print_bls = opt ~group ~l:"bls" ~h:"include alignment branch length score (BLS) for the reported region in output" (StdOpt.store_true ())
let print_anc_comp_score = opt ~group ~l:"ancComp" ~h:"include ancestral sequence composition score in output" (StdOpt.store_true ())
let print_dna = opt ~group ~l:"dna" ~h:"include DNA sequence in output, the part of the reference (first) sequence whose score is reported" (StdOpt.store_true ())
let print_aa = opt ~group ~l:"aa" ~h:"include amino acid translation in output" (StdOpt.store_true ())
let debug = opt ~l:"debug" ~group ~h:"print extra information about parameters and errors" (StdOpt.store_true ())

let cmd = OptParser.parse_argv opt_parser

if List.length cmd < 1 then
	OptParser.usage opt_parser ()
	exit (-1)
	
let paramset = List.hd cmd
let fns_input = List.tl cmd

if Opt.get orf_mode <> AsIs && Opt.get allow_ref_gaps then
	eprintf "--allowRefGaps should not be used with --orf\n"
	OptParser.usage opt_parser ()
	exit (-1)
	
if Opt.get orf_mode <> AsIs && Opt.get reading_frame = One then
	eprintf "Warning: --orf with --frames=1; are you sure you don't want to search for ORFs in three or six frames?\n"
	flush stderr

if Opt.get debug then Printexc.record_backtrace true

(******************************************************************************)

let input_mfa lines =
	try
		let rec seq sofar =
			let buf = match sofar with ((_,buf) :: _) -> buf | _ -> assert false
			if not (Enum.is_empty lines) then
				let line = String.trim (Option.get (peek lines))
				if String.length line = 0 then
					seq sofar
				else if line.[0] <> '>' then
					Buffer.add_string buf line
					junk lines
					seq sofar
				else
					hdr sofar
			else
				let enum = List.rev sofar |> Array.of_list
				let species = Array.map fst enum
				let seqs = Array.map (fun (_,buf) -> Buffer.contents buf) enum
				let seqlen = String.length seqs.(0)
				if seqlen = 0 || exists ((=) "") (Array.enum species) || exists (fun s -> String.length s <> seqlen) (Array.enum seqs) then
					failwith "empty species name or sequence, or sequence length mismatch"
				species, seqs
		and hdr sofar =
			let line = Option.get (get lines)
			if String.length line = 0 then
				hdr sofar
			else if line.[0] <> '>' then
				failwith "bad header"
			else
				let hdr = String.lchop line
				let sp = String.trim (try String.left hdr (String.index hdr '|') with Not_found -> hdr)
				seq ((sp,Buffer.create 256) :: sofar)
		hdr []
	with
		| Failure msg -> failwith ("invalid MFA alignment: " ^ msg)
		| exn -> failwith ("invalid MFA alignment: " ^ (Printexc.to_string exn))

(******************************************************************************)

let maybe_remove_ref_gaps (species,aln) =
	if Opt.get remove_ref_gaps then
		let bufs = Array.map (fun _ -> Buffer.create 256) aln
		let refseq = aln.(0)
		for j = 0 to String.length refseq - 1 do
			if refseq.[j] <> '-' then
				for i = 0 to Array.length aln - 1 do
					Buffer.add_char bufs.(i) aln.(i).[j]
		species, (Array.map Buffer.contents bufs)
	else
		species, aln

let find_orfs ?(ofs=0) dna =
	let atg = Opt.get orf_mode = ATGStop
	let codon pos = Char.uppercase dna.[pos], Char.uppercase dna.[pos+1], Char.uppercase dna.[pos+2]
	let is_start pos = match codon pos with 'A','T','G' -> true | _ -> false
	let is_stop pos = match codon pos with 'T','A','A' | 'T','A','G' | 'T','G','A' -> true | _ -> false
	
	let len = String.length dna
	let orfs = ref []
	let starts = ref []
	
	for codon_lo = ofs to len-3 do
		if (codon_lo-ofs) mod 3 = 0 then
			let codon_hi = codon_lo+2
			if (not atg && !starts = [] && not (is_stop codon_lo)) || (atg && is_start codon_lo) then starts := codon_lo :: !starts
			if codon_hi+3 < len && is_stop (codon_hi+1) then
				!starts |> List.iter
					fun start ->
						if codon_hi>start+2 then
							assert ((start-ofs) mod 3 = 0)
							assert ((codon_hi-start+1) mod 3 = 0)
							orfs := (start,codon_hi) :: !orfs
				starts := []
				
	if not atg then
		(* add stuff that goes off the end of the alignment *)
		!starts |> List.iter
			fun start ->
				let rem = len-start
				orfs := (start,start+(rem/3)*3-1) :: !orfs
			
	if Opt.get orf_mode = StopStop3 then
		(* look at shorter sub-ORFs of each stop-stop ORF *)
		let all_suborfs =
			!orfs |> List.map
				fun (lo,hi) ->
					let suborfs = ref [lo,hi]
					let codons = (hi-lo+1)/3
					let lo2 = lo+(codons/3)*3
					if lo2 > lo then suborfs := (lo2,hi) :: !suborfs
					let lo3 = lo+(2*codons/3)*3
					if lo3 > lo2 && lo3 > lo then suborfs := (lo3,hi) :: !suborfs
					!suborfs
		orfs := List.flatten all_suborfs
	
	if Opt.get orf_mode = ToFirstStop && !orfs <> [] then
		let ((firstorflo,_) as firstorf) = List.hd (List.rev !orfs)
		orfs := if firstorflo = ofs then [firstorf] else []

	if Opt.get orf_mode = FromLastStop && !orfs <> [] then
		let ((_,lastorfhi) as lastorf) = List.hd !orfs
		orfs := if len - lastorfhi <= 3 then [lastorf] else []
		
	if Opt.get orf_mode = ToOrFromStop && !orfs <> [] then
		let ((firstorflo,_) as firstorf) = List.hd (List.rev !orfs)
		let ((_,lastorfhi) as lastorf) = List.hd !orfs
		orfs := if firstorflo = ofs then [firstorf] else []		
		orfs := if len - lastorfhi <= 3 && firstorf <> lastorf then lastorf :: !orfs else !orfs
		
	!orfs |> List.rev |> List.enum |> filter
		fun (lo,hi) ->
			assert ((hi-lo+1) mod 3 = 0)
			assert ((lo-ofs) mod 3 = 0)
			(hi-lo+1)/3 >= Opt.get min_codons
	
let candidate_regions dna rcdna =
	if Opt.get orf_mode = AsIs then
		let hi = String.length dna - 1
		[ [false,0,hi];
			(if Opt.get reading_frame <> One then
				[ (false,1,hi); (false,2,hi) ] else []);
			(if Opt.get reading_frame = Six then
				[ (true,0,hi); (true,1,hi); (true,2,hi) ] else []) ] |> List.flatten |> List.enum
	else
		let aoeu a = map (fun (b,c) -> a,b,c)
		let all_orfs = ref [aoeu false (find_orfs ~ofs:0 dna)]
		if Opt.get reading_frame <> One then
			all_orfs := aoeu false (find_orfs ~ofs:1 dna) :: !all_orfs
			all_orfs := aoeu false (find_orfs ~ofs:2 dna) :: !all_orfs
		if Opt.get reading_frame = Six then
			let rcdna = Code.DNA.revcomp dna
			all_orfs := aoeu true (find_orfs ~ofs:0 rcdna) :: !all_orfs
			all_orfs := aoeu true (find_orfs ~ofs:1 rcdna) :: !all_orfs
			all_orfs := aoeu true (find_orfs ~ofs:2 rcdna) :: !all_orfs
		concat (List.enum (List.rev !all_orfs))

let pleaves ?(lo=0) ?hi t leaf_ord aln =
	let hi = match hi with Some x -> x | None -> String.length aln.(0) - 1
	
	(* construct enumeration of leaf probability vectors *)
	let rec enum starting_pos =
		let pos = ref starting_pos
		let next () =
			if !pos+2 > hi then raise Enum.No_more_elements
			let lvs =
				Array.init (T.leaves t)
					fun t_row ->
						match leaf_ord.(t_row) with
							| None -> `Marginalize
							| Some aln_row ->
								let n1 = aln.(aln_row).[!pos]
								let n2 = aln.(aln_row).[!pos+1]
								let n3 = aln.(aln_row).[!pos+2]
								let codon = n1, n2, n3
								if not (Codon.is codon) then
									`Marginalize
								else
									`Certain (Codon.index codon)
			pos := !pos + 3
			lvs
		let count () = (hi - !pos + 1)/3
		let clone () = enum !pos
		Enum.make ~next:next ~count:count ~clone:clone
	enum lo

(******************************************************************************)

(* based on the Newick tree, compute the average branch length score for each alignment column in
the interval [lo,hi] *)
let bls nt aln which_row lo hi =
	assert (hi >= lo)
	let pres = function '-' | '.' | 'N' -> false | _ -> true
	let bl_total = ref 0.
	for i = lo to hi do
		let subtree =
			Newick.subtree
				fun sp -> try pres aln.(SMap.find sp which_row).[i] with Not_found -> false
				nt
		bl_total := !bl_total +. Option.map_default Newick.total_length 0. subtree
	!bl_total /. (Newick.total_length nt *. float (hi-lo+1))

(******************************************************************************)

let u2t = function 'u' -> 't' | 'U' -> 'T' | x -> x
let fn2id fn = if fn = "" then "(STDIN)" else fn
let translate dna =
	let aalen = String.length dna / 3
	let pp = String.create aalen
	for i = 0 to aalen - 1 do
		let codon = dna.[3*i], dna.[3*i+1], dna.[3*i+2]
		if Codon.is codon then
			let aa = Codon.translate codon
			pp.[i] <- aa
		else
			pp.[i] <- '?'
	pp
	
let process_alignment (nt,t,evaluator) fn =
	try
		(* load alignment from file *)
		let lines = if fn = "" then IO.lines_of stdin else File.lines_of fn
		let (species,aln) = maybe_remove_ref_gaps (input_mfa lines)
		let aln = Array.map (String.map u2t) aln
		
		(* sanity checks *)
		if not (Opt.get allow_ref_gaps) && aln.(0) |> String.enum |> exists ((=) '-') then
			failwith "the reference sequence (first alignment row) must be ungapped"
		let rc_aln = Array.map Code.DNA.revcomp aln (* will check that everything is valid nucleotide or gap. *)

		(* determine how the alignment rows correspond to the tree leaves *)
		let t_species = (0 -- (T.leaves t - 1)) /@ T.label t |> fold (fun set sp -> SSet.add sp set) SSet.empty
		let aln_species = Array.fold_right SSet.add species SSet.empty
		let wtf = SSet.diff aln_species t_species
		if SSet.cardinal wtf > 0 then failwith ("parameters not available for species: " ^ (String.join " " (SSet.elements wtf)))
		let which_row = Array.fold_lefti (fun m i sp -> SMap.add sp i m) SMap.empty species
		let leaf_ord = Array.init (T.leaves t) (fun i -> try Some (SMap.find (T.label t i) which_row) with Not_found -> None)

		(* generate list of candidate regions within the alignment *)
		let rgns = candidate_regions aln.(0) rc_aln.(0)
		
		try
			if Enum.is_empty rgns then
				assert (Opt.get orf_mode <> AsIs)
				failwith "no sufficiently long ORFs found"
			
			(* evaluate each candidate region *)
			let rgns_scores =
				rgns |> Enum.filter_map
					fun (rc,lo,hi) ->
						try
							let aln_leaves = Array.of_enum (pleaves ~lo:lo ~hi:hi t leaf_ord (if rc then rc_aln else aln))
							let rslt = evaluator aln_leaves
							Some (rslt,rc,lo,hi)
						with
							| exn ->
								(* problem evaluating an individual region within the
								   alignment: complain, but proceed, as maybe some
								   other region will succeed. *)
								printf "%s\texception\t%d\t%d%s\t%s\n"
									fn2id fn
									lo
									hi
									if Opt.get reading_frame = Six then (if rc then "\t-" else "\t+") else ""
									Printexc.to_string exn
								if Opt.get debug then
									flush stdout
									eprintf "%s" (Printexc.get_backtrace ())
									flush stderr
								None

			if Enum.is_empty rgns_scores then failwith "no regions successfully evaluated"

			let report_score ty (rslt,rc,lo,hi) =
				printf "%s\t%s\t%.4f" (fn2id fn) ty rslt.PhyloCSFModel.score
				if Opt.get reading_frame <> One || Opt.get orf_mode <> AsIs then
					printf "\t%d\t%d" lo hi
				if Opt.get reading_frame = Six then
					printf "\t%c" (if rc then '-' else '+')
				if Opt.get print_bls then
					let rgn_bls = bls nt (if rc then rc_aln else aln) which_row lo hi
					printf "\t%.4f" rgn_bls
				if Opt.get print_anc_comp_score then
					printf "\t%.4f" rslt.PhyloCSFModel.anc_comp_score
				let refdna = String.sub (if rc then rc_aln.(0) else aln.(0)) lo (hi-lo+1)
				if Opt.get print_dna then
					printf "\t%s" refdna
				if Opt.get print_aa then
					printf "\t%s" (translate refdna)
				if Opt.get debug then
					printf "\t#"
					foreach (List.enum rslt.PhyloCSFModel.diagnostics) (fun (k,v) -> printf " %s=%s" k v)
				printf "\n"
						
			if Opt.get print_orfs then
				Enum.clone rgns_scores |> iter (report_score "orf_score(decibans)")
			reduce max rgns_scores |> report_score (if Opt.get orf_mode <> AsIs || Opt.get reading_frame <> One then "max_score(decibans)" else "score(decibans)")
		with
			| ((Assert_failure _) as exn) -> raise exn
			(* move on to the next alignment: convergence problems, no ORFs found, etc. *)
			| exn -> printf "%s\tfailure\t%s\n" (fn2id fn) (Printexc.to_string exn)
	with		(* stop everything: file not found, improperly formatted alignment, etc. *)
		| exn ->
			printf "%s\tabort\t%s\n" (fn2id fn) (Printexc.to_string exn)
			if Opt.get debug then
				flush stdout
				eprintf "%s" (Printexc.get_backtrace ())
				flush stderr
			exit (-1)
	flush stdout
	
(******************************************************************************)

(* parse a key<tab>value\n file format *)
let parse_kv lines =
	let f kv line =
		let line = String.strip line
		if String.length line = 0 || line.[0] = '#' then
			kv
		else
			try
				let (k,v) = String.split line "\t"
				PMap.add k v kv
			with Not_found -> PMap.add line "" kv
	fold f PMap.empty lines

let initialize_strategy () =
	let paramfile ?(required=true) suffix =
		let fn =
			if (try ignore (String.index paramset '/'); false with Not_found -> true) then
				try
					let base = Sys.getenv "PHYLOCSF_BASE"
					if not (Sys.file_exists base && Sys.is_directory base) then
						raise Not_found
					Filename.concat (Filename.concat base "PhyloCSF_Parameters") (paramset ^ suffix)
				with
					| _ -> failwith "PHYLOCSF_BASE environment variable must be set to the root directory of the source or executable distribution"
			else
				(* paramset is a relative or absolute path *)
				paramset ^ suffix
		if required && not (Sys.file_exists fn) then failwith (sprintf "could not find required parameter file %s" fn)
		fn
	
	let fn_tree = paramfile ".nh"
	let nt = File.with_file_in fn_tree (fun input -> NewickParser.parse NewickLexer.token (Lexing.from_input input))

	let desired_species = if Opt.is_set desired_species then String.nsplit (Opt.get desired_species) "," |> List.fold_left (flip SSet.add) SSet.empty else SSet.empty
	let snt =
		if desired_species = SSet.empty then
			nt
		else
			match Newick.subtree (flip SSet.mem desired_species) nt with
				| Some snt when Newick.leaves snt > 1 -> snt
				| _ -> failwith "specify at least two available --species"

	let t = T.of_newick snt
	
	let evaluator =
		match Opt.get strategy with
			| PhyloCSF substrategy ->
				let fn_ecm_c = paramfile "_coding.ECM"
				let fn_ecm_nc = paramfile "_noncoding.ECM"	

				let s1, pi1 = ECM.import_parameters fn_ecm_c
				let s2, pi2 = ECM.import_parameters fn_ecm_nc
				let model = PhyloCSFModel.make s1 pi1 s2 pi2 t

				fun leaves ->
					PhyloCSFModel.score
						match substrategy with `MaxLik -> PhyloCSFModel.MaxLik | `FixedLik -> PhyloCSFModel.FixedLik
						model
						leaves
			| OmegaTest ->
				let fn_config = paramfile ~required:false "_omega"
				
				let omega_H1, sigma_H1 =
					if Sys.file_exists fn_config then
						try
							let kv = parse_kv (File.lines_of fn_config)
							Some (float_of_string (PMap.find "omega_H1" kv)), Some (float_of_string (PMap.find "sigma_H1" kv))
						with
							| exn -> failwith (sprintf "configuration file %s exists but does not specify valid omega_H1 and sigma_H1 values" fn_config)
					else
						None, None
			
				fun leaves -> (OmegaModel.score ?omega_H1 ?sigma_H1 t leaves)
			| Nop -> fun leaves -> { PhyloCSFModel.score = 0.0; anc_comp_score = 0.0; diagnostics = [] }
	snt, t, evaluator
	
let main () =
	let strategy = initialize_strategy ()
	
	let fns_alignments =
		if Opt.get filenames then
			if fns_input = [] then
				if Unix.isatty Unix.stdin then
					printf "Input alignment filename(s) then EOF (Ctrl+D): "
					flush stdout
				IO.lines_of stdin
			else
				(List.enum fns_input) /@ File.lines_of |> Enum.concat
		else if fns_input = [] then
			if Unix.isatty Unix.stdin then
				printf "Input alignment then EOF (Ctrl+D):\n"
				flush stdout
			Enum.singleton ""
		else
			List.enum fns_input
		
	foreach fns_alignments (process_alignment strategy)

main ()
