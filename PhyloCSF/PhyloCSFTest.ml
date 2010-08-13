(* Test harness for PhyloCSF executable...feeds it simulated alignments of 5'UTR+ORF+3'UTR, to make sure we can recover the ORF. 
_build/PhyloCSFTest.native "_build/PhyloCSF.native --frames=6 --orfs=ATGStop --strategy=fixed" ../PhyloCSFParameters/12flies
*)

open Batteries_uni
open OptParse
open Printf
open CamlPaml

Gsl_error.init ()
Random.self_init ()

module Codon = CamlPaml.Code.Codon64

let opt_parser = OptParser.make ~usage:"%prog \"/path/to/PhyloCSF --phyloCSF --options\" /path/to/parameters" ()
let opt ?group ?h ?hide ?s ?short_names ?l ?long_names x = OptParser.add opt_parser ?group ?help:h ?hide ?short_name:s ?long_name:l x; x

let min_utr = opt ~l:"minUTR" (StdOpt.int_option ~default:10 ())
let max_utr = opt ~l:"maxUTR" (StdOpt.int_option ~default:50 ())
let min_cds = opt ~l:"minCDS" (StdOpt.int_option ~default:40 ())
let max_cds = opt ~l:"maxCDS" (StdOpt.int_option ~default:200 ())
let randomize_frame = opt ~l:"constantFrame" (StdOpt.store_false ())
let randomize_strand = opt ~l:"constantStrand" (StdOpt.store_false ())

let cmd = OptParser.parse_argv opt_parser

if List.length cmd < 2 then
	OptParser.usage opt_parser ()
	exit (-1) 

let fn_exe = List.nth cmd 0
let fp_params = List.nth cmd 1
	
(******************************************************************************)

let load_parameters () =
	let fn_tree = fp_params ^ ".nh"
	let fn_ecm_c = fp_params ^ "_coding.ECM"
	let fn_ecm_nc = fp_params ^ "_noncoding.ECM"
	if not (List.for_all Sys.file_exists [fn_tree; fn_ecm_c; fn_ecm_nc]) then
		failwith (sprintf "Could not find required parameter files prefixed by %s\n" fp_params)

	let nt = File.with_file_in fn_tree (fun input -> NewickParser.parse NewickLexer.token (Lexing.from_input input))
	let t = T.of_newick nt
	let s1, pi1 = ECM.import_parameters fn_ecm_c
	let s2, pi2 = ECM.import_parameters fn_ecm_nc
	let model = PhyloCSFModel.make s1 pi1 s2 pi2 t
	nt, t, model
	
let alignment_column_producer model_instance =
	let m = PhyloModel.P14n.model model_instance
	let rec next () =
		let leaves = PhyloModel.simulate m ()
		if Codon.is_stop (Codon.of_index leaves.(0)) then
			next ()
		else
			Array.init (T.leaves (PhyloModel.tree (PhyloModel.component m 0)))
				fun i ->
					let which_codon = leaves.(i)
					let n1,n2,n3 = Codon.of_index which_codon
					String.of_list [n1; n2; n3]
	let rec enum () = Enum.make ~next:next ~count:(fun () -> raise Enum.Infinite_enum) ~clone:(fun () -> enum ())
	enum ()
	
let full_alignment_columns { PhyloCSFModel.coding_model = coding_model; noncoding_model = noncoding_model } =
	let leaves = (T.leaves (PhyloModel.tree (PhyloModel.component (PhyloModel.P14n.model coding_model) 0)))
	
	let sites5' = if Opt.get max_utr > 0 then (Opt.get min_utr) + (Random.int (Opt.get max_utr - Opt.get min_utr)) else 0
	let sites3' = if Opt.get max_utr > 0 then (Opt.get min_utr) + (Random.int (Opt.get max_utr - Opt.get min_utr)) else 0
	let sitesCDS = if Opt.get max_cds > 0 then (Opt.get min_cds) + (Random.int (Opt.get max_cds - Opt.get min_cds)) else 0
	
	let utr_producer = alignment_column_producer noncoding_model
	let cds_producer = alignment_column_producer coding_model
	
	let site1 =
		if Opt.get max_utr = 0 then
			None
		else
			let site = Option.get (get utr_producer)
			if Opt.get randomize_frame then
				let trim = Random.int 3
				site |> Array.iteri (fun i s -> site.(i) <- String.sub s 0 trim)
			Some site
	
	let atg = Array.make leaves "ATG"
	let stop = Array.init leaves (fun _ -> match Random.int 3 with 0 -> "TAA" | 1 -> "TAG" | 2 -> "TGA" | _ -> assert false)
	
	let orf_lo = (match site1 with Some site -> String.length site.(0) | None -> 0) + (sites5' * 3)
	let orf_hi = orf_lo + 2 + (sitesCDS * 3)
	printf "\t%d\t%d" orf_lo orf_hi

	Enum.concat (List.enum [ (match site1 with Some site -> Enum.singleton site | None -> Enum.empty ());
								Enum.take sites5' utr_producer;
								Enum.singleton atg;
								Enum.take sitesCDS cds_producer;
								Enum.singleton stop;
								Enum.take sites3' utr_producer ])

let stringify_alignment_columns cols =
	let col0 = Option.get (Enum.peek cols)
	let rows = Array.length col0
	let bufs = Array.init rows (fun _ -> Buffer.create 1000)
	cols |> iter
		fun col ->
			for i = 0 to rows-1 do
				Buffer.add_string bufs.(i) col.(i)
	Array.map Buffer.contents bufs
	
let maybe_revcomp seqs =
	if Opt.get randomize_strand && Random.bool () then
		printf "\t-"
		Array.map Code.DNA.revcomp seqs
	else
		printf "\t+"
		seqs

let mfa headers seqs = 
	assert (Array.length headers = Array.length seqs)
	let buf = Buffer.create 1000
	Array.iter2 
		fun header seq ->
			Buffer.add_string buf (sprintf ">%s\n" header)
			Buffer.add_string buf seq
			Buffer.add_char buf '\n'
		headers
		seqs
	buf
	
let open_phylocsf () =
	let cmd = sprintf "%s %s" fn_exe fp_params
	Unix.open_process_out ~cleanup:true cmd
	
let main () =
	let _, t, model = load_parameters ()
	
	let headers = Array.init (T.leaves t) (T.label t)
	
	while true do
		printf "\ttruth\t\t\t"
		let aln = mfa headers (maybe_revcomp (stringify_alignment_columns (full_alignment_columns model)))
		printf "\n"
		flush stdout
		flush stderr
		let out = open_phylocsf ()
		Buffer.output_buffer out aln
		flush out
		ignore (Unix.close_process_out out)
		printf "\n"
		flush stdout

main ()
