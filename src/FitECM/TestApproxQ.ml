open OUnit

let test1 () = todo "test1 todo"

let all_tests = "test1" >:: test1

run_test_tt_main all_tests
