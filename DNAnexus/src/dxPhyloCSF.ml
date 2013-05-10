open Batteries
open Printf
open JSON.Operators
open DNAnexus

Printexc.record_backtrace true

(* abbreviated JSON casts *)
module J = struct
  include JSON
  let of_str_list lst = of_list (List.map (fun x -> `String x) lst)
  let of_str_assoc lst = of_assoc (List.map (fun (k,v) -> k, `String v) lst)

(* return a list of (starting,limit) tuples *)
let partition_rows length rows_per_job =  
  if length = 0 then []
  else
    (* make distribution more uniform *)
    let jobs = 1 + (length - 1) / rows_per_job
    let k = int_of_float (ceil (float length /. float jobs))
    (0 --^ jobs) /@ (fun i -> (i*k), (min k (length-i*k))) |> List.of_enum

let make_scores_gtable name app_id input =
  let optional_cols = []
  GTable.make_new (J.of_assoc [
    "name", `String name;
    "columns", GTable.json_of_columns (Array.of_list (("name",`String) :: ("score",`Float) :: optional_cols));
    "indices", GTable.json_of_indices ["name", `Lexicographic ["name",`Asc,None]];
    "details", J.of_assoc [
      "PhyloCSF", make_link app_id;
      "PhyloCSF_input", input
    ]
  ])

let get_species_names alignments_table = 
  try
    let alignments_cols = GTable.columns alignments_table
    let n = Array.length alignments_cols
    assert (n >= 4)
    Array.sub alignments_cols 2 (n-2) |> Array.map
      function
        | (sp,`String) when sp <> "" -> sp
        | _ -> assert false
  with 
    | exn ->
      eprintf "%s\n%s\n" (Printexc.to_string exn) (Printexc.get_backtrace ())
      raise (DNAnexus.AppInternalError "Could not understand input alignments table. Please make sure it was generated with the maf_stitcher app.")

let run_PhyloCSF species_set species_names alignment_row =
  let n = Array.length species_names
  assert (Array.length alignment_row = n+3)
  let alignment_name = match alignment_row.(1) with
    | `String nm -> nm
    | _ -> assert false
  let seqs =
    Array.sub alignment_row 3 n |> Array.map
      function
        | `String sp -> sp
        | _ -> assert false

  (* TODO: keep PhyloCSF online to avoid reinstantiating phylo-models for every alignment *)

  let phylocsf_out, phylocsf_in = Unix.open_process ("/PhyloCSF " ^ species_set ^ " --removeRefGaps")
  for i = 0 to n-1 do
    fprintf phylocsf_in ">%s\n" species_names.(i)
    fprintf phylocsf_in "%s\n" seqs.(i)
  flush phylocsf_in
  close_out phylocsf_in

  let result = IO.read_all phylocsf_out

  match Unix.close_process (phylocsf_out, phylocsf_in) with 
    | Unix.WEXITED 0 -> ()
    | _ ->
        printf "%s\n" result
        raise (DNAnexus.AppInternalError ("PhyloCSF exited abnormally; see the job log for more details."))

  List.of_enum
    IO.lines_of (IO.input_string result) |> Enum.filter_map
      fun line ->
        if String.length line = 0 || line.[0] = '#' then
          print_endline line
          None
        else
          printf "%s\n" line
          let flds = Array.of_list (String.nsplit line "\t")
          assert (Array.length flds >= 2)
          assert (flds.(1) = "score(decibans)")
          Some [| `String alignment_name; `Float (float_of_string flds.(2)) |]

(* main job: instantiate the job tree *)
let main input =
  let alignments_table = GTable.bind_link (input$"alignments")
  let alignments_desc = GTable.describe alignments_table

  (* resolve which parameter set to use *)
  let species_set =
    if input$?"species_set" then J.string (input$"species_set")
    else
      let alignments_deets = GTable.get_details alignments_table
      let alignments_species =
        try
          J.string (alignments_deets$"maf_stitcher_input"$"species_set")
        with _ ->
          raise (AppError "Please provide the species_set input (could not detect which PhyloCSF parameters to use for these alignments)")
      match alignments_species with
        | "hg19_33mammals" -> "33mammals"
        | "dm3_12flies" -> "12flies"
        | _ -> raise (AppError (sprintf "No PhyloCSF parameters available for maf_stitcher species set: %s. Regenerate the alignments with a supported species set, or override this check by providing the species_set input" alignments_species))

  (* verify the species names can be extracted from the alignments table.
     TODO: also verify that the obtained species names match those supported by species_set *)
  ignore (get_species_names alignments_table)

  (* get own ID, for details breadcrumbs *)
  let job_desc = api_call [job_id (); "describe"] J.empty
  let app_id =
    if job_desc$?"app" then J.string (job_desc$"app")
    else J.string (job_desc$"applet")

  (* initialize output GTable *)
  let output_name =
    if input$?"output_name" then J.string (input$"output_name")
    else J.string (alignments_desc$"name") ^ " PhyloCSF"
  let output_gtable = make_scores_gtable output_name app_id input

  (* partition alignments and launch subjobs *)
  let subjobs =
    partition_rows (JSON.int (alignments_desc$"length")) (JSON.int (input$"alignments_per_job")) |> List.map
      fun (starting, limit) ->
        let subjob = J.of_assoc [
          "species_set", `String species_set;
          "starting", `Int starting;
          "limit", `Int limit;
          "scores", `String (GTable.id output_gtable)
        ]
        new_job "process" (input $+ ("process",subjob))

  (* set up postprocessing job and output *)
  let postprocess_input =
    (input
      $+ ("postprocess", J.of_assoc [
            "scores", `String (GTable.id output_gtable);
            "errors", J.of_list (List.map (fun subjob -> make_jobref subjob "errors") subjobs)]))
  let postprocess = new_job "postprocess" postprocess_input

  J.of_assoc [
    "scores", make_link (GTable.id output_gtable);
    "errors", make_jobref postprocess "errors"
  ]

(* subjob: process some of the alignments *)
let process input =
  let alignments_table = GTable.bind_link (input$"alignments")
  let species_names = get_species_names alignments_table
  let starting = J.int (input$"process"$"starting")
  let limit = J.int (input$"process"$"limit")
  let scores_table = GTable.bind (None,J.string (input$"process"$"scores"))

  let worker_process alignment =
    try
      run_PhyloCSF (J.string (input$"process"$"species_set")) species_names alignment
    with
      | DNAnexus.AppError msg ->
          raise
            ForkWork.ChildExn [
              "AppError"; J.string ((GTable.json_of_row alignment)$@1); msg
            ]
      | exn ->
          raise
            ForkWork.ChildExn [
              "AppInternalError";
              J.string ((GTable.json_of_row alignment)$@1);
              Printexc.to_string exn;
              Printexc.get_backtrace ()
            ]

  let fw = ForkWork.manager ()
  let rows = GTable.iterate_rows ~starting ~limit alignments_table
  let promises = rows /@ (ForkWork.fork fw worker_process)

  let errors =
    promises |> Enum.fold
      fun n ans ->
        match ForkWork.await_result fw ans with
          | `OK [] -> n+1
          | `OK scores ->
            scores |> List.iter (fun score_row -> GTable.add_row scores_table score_row)
            n
          | `Exn ["AppError"; name; msg] ->
              eprintf "AppError on %s: %s\n" name msg
              raise (DNAnexus.AppError (sprintf "%s: %s" name msg))
          | `Exn ["AppInternalError"; name; msg; bt] ->
              eprintf "AppInternalError on %s: %s\n" name msg
              if bt <> "" then eprintf "%s\n" bt
              raise (DNAnexus.AppInternalError (sprintf "%s: %s" name msg))
          | _ -> assert false
      0

  GTable.flush_rows scores_table
  J.of_assoc ["errors", `Int errors]

(* finish-up job: collect stats & initiate GTable closure *)
let postprocess input =
  let gtable = GTable.bind (None, JSON.string (input$"postprocess"$"scores"))

  let stats =
    Enum.fold
      fun failures j -> failures + JSON.int j
      0
      Vect.enum (JSON.array (input$"postprocess"$"errors"))
  let stats_json = J.of_assoc ["errors", `Int stats]

  GTable.set_details gtable (GTable.get_details gtable $+ ("PhyloCSF_stats", stats_json))

  GTable.close gtable
  J.of_assoc ["scores", make_link (GTable.id gtable); "errors", `Int stats]

(* entry point *)
job_main
  fun input ->
    if input $? "process" then process input
    else if input $? "postprocess" then postprocess input
    else main input
