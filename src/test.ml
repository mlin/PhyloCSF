open Batteries
open Printf
open OUnit
open Should

let dn_here = Filename.dirname Sys.argv.(0)

let fn_PhyloCSF = Filename.concat dn_here "PhyloCSF.native"

let slow () = skip_if (try ignore (Sys.getenv "SKIP_SLOW"); true with Not_found -> false) "SKIP_SLOW"

(* test results on the three bundled example alignments *)

let run_PhyloCSF species params =
  let cmd = sprintf "%s %s %s" fn_PhyloCSF (Filename.concat dn_here ("../PhyloCSF_Parameters/" ^ species)) params

  let phylocsf_in = Unix.open_process_in ~cleanup:true cmd

  let phylocsf_answer = input_line phylocsf_in
  match Unix.close_process_in phylocsf_in with
    | Unix.WEXITED 0 -> String.nsplit phylocsf_answer "\t"
    | _ -> raise Exit

let talAA () =
  let ans = run_PhyloCSF "12flies" "../PhyloCSF_Examples/tal-AA.fa"
  print_endline (String.join "\t" ans)
  List.nth ans 1 $hould # equal "score(decibans)"
  float_of_string (List.nth ans 2) $hould # be # within (297.62,297.63)

let aldh2_ex5_out () =
  let ans = run_PhyloCSF "29mammals" "../PhyloCSF_Examples/ALDH2.exon5.fa"
  print_endline (String.join "\t" ans)
  List.nth ans 1 $hould # equal "score(decibans)"
  float_of_string (List.nth ans 2) $hould # be # within (-178.93,-178.92)

let aldh2_ex5_in () =
  let abbreviate = (try ignore (Sys.getenv "SKIP_SLOW"); true with Not_found -> false)
  let ans = run_PhyloCSF "29mammals" ("../PhyloCSF_Examples/ALDH2.exon5.fa --frames=" ^ (if abbreviate then "3" else "6"))
  print_endline (String.join "\t" ans)
  List.nth ans 1 $hould # equal "max_score(decibans)"
  float_of_string (List.nth ans 2) $hould # be # within (218.26,218.27)
  int_of_string (List.nth ans 3) $hould # equal 1
  int_of_string (List.nth ans 4) $hould # equal 111
  if not abbreviate then List.nth ans 5 $hould # equal "+"

let aldh2_mRNA () =
  slow ()
  let ans = run_PhyloCSF "29mammals" "../PhyloCSF_Examples/Aldh2.mRNA.fa --orf=ATGStop --frames=3 --removeRefGaps --aa"
  print_endline (String.join "\t" ans)
  List.nth ans 1 $hould # equal "max_score(decibans)"
  float_of_string (List.nth ans 2) $hould # be # within (2013.92,2013.93)
  int_of_string (List.nth ans 3) $hould # equal 343
  int_of_string (List.nth ans 4) $hould # equal 1899
  string (List.nth ans 5) $hould # be # matching (Str.regexp "^MLRAALTTVRRGPRLSRLLSAAA.*")

let example_tests = "examples" >::: [
    "tal-AA" >:: talAA;
    "ALDH2 ex5 out-of-frame" >:: aldh2_ex5_out;
    "ALDH2 ex5 in-frame" >:: aldh2_ex5_in;
    "Aldh2 mRNA" >:: aldh2_mRNA
  ]

(* simulations
   FIXME: currently there is no verification of the results except that PhyloCSF
          successfully exits; you have to eyeball the test stdout. *)

let fn_testSim = Filename.concat dn_here "testSim.native"

let check cmd = match Sys.command cmd with
  | 0 -> ()
  | c when (c=130 || c=255) -> exit c (* ctrl-C *)
  | c -> failwith (sprintf "Child process exit code %d" c)

let sim_AsIs_fixed () =
  check (sprintf "%s --maxUTR=0 --constantFrame --constantStrand 12flies ' --strategy=fixed'" fn_testSim)

let sim_AsIs_mle () =
  slow ()
  check (sprintf "%s --maxUTR=0 --constantFrame --constantStrand 12flies ' --strategy=mle'" fn_testSim)

let sim_ATGStop_fixed () =
  check (sprintf "%s 12flies ' --orf=ATGStop --frames=6 --strategy=fixed'" fn_testSim)

let sim_ATGStop_mle () =
  slow ()
  check (sprintf "%s 12flies ' --orf=ATGStop --frames=6 --strategy=mle'" fn_testSim)


let sim_tests = "simulations" >::: [
    "AsIs_fixed" >:: sim_AsIs_fixed;
    "AsIs_mle" >:: sim_AsIs_mle;
    "ATGStop_fixed" >:: sim_ATGStop_fixed;
    "ATGStop_mle" >:: sim_ATGStop_mle
  ]

let all_tests =
    "PhyloCSF" >::: [
      sim_tests;
      example_tests
    ]

run_test_tt_main all_tests
