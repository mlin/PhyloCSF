open OUnit
open Printf

let fs = sprintf "%.4f"

let assert_equal_vector ?msg ?epsilon = assert_equal ~cmp:(fun ar1 ar2 -> List.for_all (fun (x,y) -> cmp_float ?epsilon x y) (List.combine (Array.to_list ar1) (Array.to_list ar2))) ~printer:(fun ar -> String.concat " " (Array.to_list (Array.map fs ar))) ?msg

module TestInfer = struct
	let nt = NewickParser.parse NewickLexer.token (Lexing.from_string "(A,(B,C))")
	let t = T.of_newick nt
	
	let sm0 = Gsl.Matrix.of_arrays [| [| 0.8; 0.2; |]; [| 0.25; 0.75 |] |]
	let sm1 = Gsl.Matrix.of_arrays [| [| 0.9; 0.1; |]; [| 0.85; 0.15 |] |]
	let sms = [| sm0; sm1; sm1; sm1 |]
	let prior = [| 0.6; 0.4 |]

	let ( + ) = ( +. )
	let ( * ) = ( *. )
	let ( / ) = ( /. )
	
	let case lvs () =
		let alf i = [| if lvs.(i) = 0 then 1. else 0.; if lvs.(i) = 1 then 1. else 0. |]
		let a0 = alf 0
		let a1 = alf 1
		let a2 = alf 2		
		let a3 = [| (sm1.{0,0} * a1.(0) + sm1.{0,1} * a1.(1)) * (sm1.{0,0} * a2.(0) + sm1.{0,1} * a2.(1));
					(sm1.{1,0} * a1.(0) + sm1.{1,1} * a1.(1)) * (sm1.{1,0} * a2.(0) + sm1.{1,1} * a2.(1)) |]		
		let a4 = [| (sm0.{0,0} * a0.(0) + sm0.{0,1} * a0.(1)) * (sm1.{0,0} * a3.(0) + sm1.{0,1} * a3.(1));
					(sm0.{1,0} * a0.(0) + sm0.{1,1} * a0.(1)) * (sm1.{1,0} * a3.(0) + sm1.{1,1} * a3.(1)) |]
		let z = prior.(0) * a4.(0) + prior.(1) * a4.(1)
		
		let inter = PhyloLik.prepare t sms prior (Array.map (fun x -> `Certain x) lvs)
		assert_equal ~cmp:(cmp_float ~epsilon:0.001) ~printer:fs ~msg:"likelihood" z (PhyloLik.likelihood inter)
		
		let b4 = prior
		let b3 = [| (b4.(0) * sm0.{0,0} * a0.(0) + b4.(0) * sm0.{0,1} * a0.(1)) * sm1.{0,0} +
						(b4.(1) * sm0.{1,0} * a0.(0) + b4.(1) * sm0.{1,1} * a0.(1)) * sm1.{1,0};
					(b4.(0) * sm0.{0,0} * a0.(0) + b4.(0) * sm0.{0,1} * a0.(1)) * sm1.{0,1} +
						(b4.(1) * sm0.{1,0} * a0.(0) + b4.(1) * sm0.{1,1} * a0.(1)) * sm1.{1,1}
					 |]
		
		assert_equal_vector ~msg:"posterior(root)" ~epsilon:0.001 [| b4.(0) * a4.(0) / z; b4.(1) * a4.(1) / z |] (PhyloLik.node_posterior inter 4)
		assert_equal_vector ~msg:"posterior(internal)" ~epsilon:0.001 [| b3.(0) * a3.(0) / z; b3.(1) * a3.(1) / z |] (PhyloLik.node_posterior inter 3)
		
		
	let tests = "PhyloLik" >::: [	"000" >:: case [| 0; 0; 0 |];
									"001" >:: case [| 0; 0; 1 |];
									"010" >:: case [| 0; 1; 0 |];
									"011" >:: case [| 0; 1; 1 |];
									"100" >:: case [| 1; 0; 0 |];
									"101" >:: case [| 1; 0; 1 |];
									"110" >:: case [| 1; 1; 0 |];
									"111" >:: case [| 1; 1; 1 |];
								]

module BeagleTests = struct
	(* ports of examples found in the BEAGLE source: http://code.google.com/p/beagle-lib/source/browse/ *)
	
	module TinyTest = struct
		let human = "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCCCCTTATTTCTACTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTG-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCGGAGCTTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCT-AAATCCTATATATCTTA-CACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGA"
		let chimp = "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTCAAATCCCCTTATTTCTACTAGGACTATAAGAATCGAACTCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTAAGCACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCTCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AACCCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCGACACCTCTTTACAGTGA"
		let gorilla = "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGAGGTTTAAACCCCCTTATTTCTACTAGGACTATGAGAATTGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTGTCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTCACATCCTTCCCGTACTAAGAAATTTAGGTTAAACATAGACCAAGAGCCTTCAAAGCCCTTAGTAAGTTA-CAACACTTAATTTCTGTAAGGACTGCAAAACCCTACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCAATGGGACTCAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAT-TCACCTCGGAGCTTGGTAAAAAGAGGCCCAGCCTCTGTCTTTAGATTTACAGTCCAATGCCTTA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGACCTTCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAGGTTAAATTACGGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCGTTAACCTTTTAAGTTAAAGATTAAGAGTATCGGCACCTCTTTGCAGTGA"
		let t = T.of_newick (NewickParser.parse NewickLexer.token (Lexing.from_string "((human:0.1,chimp:0.1):0.1,gorilla:0.2);"))
		let q =
			Q.Diag.scale
				Q.Diag.of_Q [|	[| (-3.); 1.; 1.; 1.; |];
								[| 1.; (-3.); 1.; 1.; |];
								[| 1.; 1.; (-3.); 1.; |];
								[| 1.; 1.; 1.; (-3.); |] |]
				(1. /. 3.)

		let qdiagtests = "Q.Diag" >::: [TestCase (fun () -> assert_bool "Q.Diag.reversible" (Q.Diag.reversible q));
										TestCase (fun () -> assert_equal_vector ~msg:"Q.Diag.equilibrium" (Q.Diag.equilibrium q) [| 0.25; 0.25; 0.25; 0.25 |])]
										
		let likelihood which_q () =
			let m = PhyloModel.make t [| which_q |]
			let ll = ref 0.
			for i = 0 to String.length human - 1 do
				let lf = Array.map (fun ch -> if Code.DNA.is ch then `Certain (Code.DNA.index ch) else `Marginalize) [| human.[i]; chimp.[i]; gorilla.[i] |]
				ll := !ll +. log (PhyloLik.likelihood (PhyloModel.prepare_lik m lf))
			assert_equal ~cmp:(cmp_float ~epsilon:0.001) ~printer:string_of_float ~msg:"mismatch" (-1574.63623) !ll

		let qc =
			Q.Diag.scale
				Q.Diag.of_Q ~force_complex:true [|
						[| (-3.); 1.; 1.; 1.; |];
						[| 1.; (-3.); 1.; 1.; |];
						[| 1.; 1.; (-3.); 1.; |];
						[| 1.; 1.; 1.; (-3.); |]
					|]
				(1. /. 3.)

		let qcdiagtests = "Q.Diag (complex)" >::: [
				TestCase (fun () -> assert_bool "Q.Diag.reversible" (Q.Diag.reversible qc));
				TestCase (fun () -> assert_equal_vector ~msg:"Q.Diag.equilibrium" (Q.Diag.equilibrium qc) [| 0.25; 0.25; 0.25; 0.25 |])
			]
					
		let tests = "TinyTest" >::: [ qdiagtests; "likelihood" >:: likelihood q;
									  qcdiagtests; "likelihood (complex)" >:: likelihood qc ]

	let tests = "BeagleTest" >::: [ TinyTest.tests ]

let all_tests = TestList [TestInfer.tests; BeagleTests.tests]

ignore (run_test_tt_main all_tests)
