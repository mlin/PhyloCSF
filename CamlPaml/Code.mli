(** abstractions for nucleotides, amino acids, and codons *)

(** common signature for a code *)
module type S = sig
	type t
	
	(** alphabet size (e.g. 4 for nucleotides, 20 for amino acids) *)
	val dim : int
	
	(** return the index of the given codeword (between 0 and dim-1) *)
	val index : t -> int
	
	(** return a codeword from its index *)
	val of_index : int -> t
	
	val compare : t -> t -> int

(** 4 DNA nucleotides *)
module DNA : sig
	type t = char
	val dim : int
	val index : t -> int
	val of_index : int -> t
	val compare : t -> t -> int

	val is : char -> bool
	
	val comp : t -> t
	val revcomp : string -> string

(** 20 amino acids *)
module AA : sig
	type t = char
	val dim : int
	val index : t -> int
	val of_index : int -> t
	val compare : t -> t -> int
	
	val is : char -> bool
	
	val blosum62 : t -> t -> int

(** 64 codons (alphabetical, including stop codons) *)
module Codon64 : sig
	type t = DNA.t*DNA.t*DNA.t
	val dim : int
	val index : t -> int
	val of_index : int -> t
	val compare : t -> t -> int

	val is : char*char*char -> bool
	
	val is_stop : t -> bool
	val is_stop_index : int -> bool
	
	val translate : t -> AA.t
	val translate_index : int -> AA.t
	
(** @deprecated old name for Codon64 *)
module Codon1 : sig
	type t = DNA.t*DNA.t*DNA.t
	val dim : int
	val index : t -> int
	val of_index : int -> t
	val compare : t -> t -> int

	val is : char*char*char -> bool
	
	val is_stop : t -> bool
	val is_stop_index : int -> bool
	
	val translate : t -> AA.t
	val translate_index : int -> AA.t


(** 61 sense codons in a distance-based order *)
module Codon61 : sig
	type t = DNA.t*DNA.t*DNA.t
	val dim : int
	val index : t -> int
	val of_index : int -> t
	val compare : t -> t -> int

	val is : char*char*char -> bool
	
	val translate : t -> AA.t
	val translate_index : int -> AA.t

(** @deprecated old name for Codon61 *)
module Codon2 : sig
	type t = DNA.t*DNA.t*DNA.t
	val dim : int
	val index : t -> int
	val of_index : int -> t
	val compare : t -> t -> int

	val is : char*char*char -> bool
	
	val translate : t -> AA.t
	val translate_index : int -> AA.t
