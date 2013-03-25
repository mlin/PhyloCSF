open Printf
open Batteries
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
  GTable.make_new (J.of_assoc [
    "name", `String name;
    "columns", GTable.json_of_columns (Array.of_list (("name",`String) :: ("score",`Float) :: species_cols));
    "indices", GTable.json_of_indices ["name", `Lexicographic ["name",`Asc,None]];
    "details", J.of_assoc [
      "PhyloCSF", make_link app_id;
      "PhyloCSF_input", input
    ]
  ])

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
    "alignments", make_link None (GTable.id output_gtable);
    "errors", make_jobref postprocess "errors"
  ]

(* subjob: process some of the alignments *)
let process input =
  let alignments_table = GTable.bind_link (input$"alignments")
  let starting = J.int (input$"subjob"$"starting")
  let limit = J.int (input$"subjob"$"limit")
  let scores_table = GTable.bind (J.string (input$"subjob"$"scores"))

  let worker_process alignment =
    try
      run_PhyloCSF alignment
    with
      | DNAnexus.AppError msg ->
          raise
            ForkWork.ChildExn [
              "AppError"; region.rgn_id; msg
            ]
      | exn ->
          raise
            ForkWork.ChildExn [
              "AppInternalError";
              J.string (alignment$@0);
              Printexc.to_string exn;
              Printexc.get_backtrace ()
            ]

  let fw = ForkWork.manager ()
  let rows = GTable.iterate_rows ~starting ~limit alignments_table
  let promises = rows /@ (ForkWork.fork mgr worker_process)

  let errors =
    promises |> Enum.fold
      fun n ans ->
        match ans with
          | `OK (Some score_row) -> GTable.add_rows scores_table score_row; n
          | `OK None -> n+1
          | `Exn ["AppError"; name; msg] ->
              eprintf "AppError on %s: %s\n" name msg
              raise (DNAnexus.AppError (sprintf "%s: %s" name msg))
          | `Exn ["AppInternalError"; name; msg; bt]
              eprintf "AppInternalError on %s: %s\n" name msg
              if bt <> "" then eprintf "%s\n" bt
              raise (DNAnexus.AppInternalError (sprintf "%s: %s" name msg))
          | _ -> assert false
      0

  GTable.flush_rows scores_table
  J.of_assoc ["errors", `Int errors]

(* finish-up job: collect stats & initiate GTable closure *)
let finishjob input =
  let gtable = GTable.bind (None, JSON.string (input$"finishjob"$"alignments"))

  let stats =
    Enum.fold
      fun failures j -> failures + JSON.int (j$"failures")
      0
      Vect.enum (JSON.array (input$"finishjob"$"stats"))
  let stats_json = J.of_assoc ["failures", `Int stats]

  GTable.set_details gtable (GTable.get_details gtable $+ ("stitch_stats", stats_json))

  GTable.close gtable
  J.of_assoc ["alignments", make_link (GTable.id gtable)]

(* entry point *)
job_main
  fun input ->
    if input $? "subjob" then subjob input
    else if input $? "finishjob" then finishjob input
    else mainjob input
