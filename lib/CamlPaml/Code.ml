open Printf
open List

module type S = sig
	type t
	val dim : int
	val index : t -> int
	val of_index : int -> t
	val compare : t -> t -> int

module DNA = struct
	type t = char

	let dim = 4

	let compare = compare

	let is = function
		| 'A' | 'a'
		| 'C' | 'c'
		| 'G' | 'g'
		| 'T' | 't' -> true
		| _ -> false

	let index = function
		| 'A' | 'a' -> 0
		| 'C' | 'c' -> 1
		| 'G' | 'g' -> 2
		| 'T' | 't' -> 3
		| _ as nuc -> invalid_arg (sprintf "unrecognized nucleotide %c" nuc)

	let of_index = function
		| 0 -> 'A'
		| 1 -> 'C'
		| 2 -> 'G'
		| 3 -> 'T'
		| _ -> invalid_arg "CamlPaml.Code.DNA.of_index"

	let comp = function
		| 'A' -> 'T'
		| 'G' -> 'C'
		| 'C' -> 'G'
		| 'T' -> 'A'
		| 'a' -> 't'
		| 'g' -> 'c'
		| 'c' -> 'g'
		| 't' -> 'a'
		| 'N' -> 'N'
		| 'n' -> 'n'
		| '-' -> '-'
		| _ as nuc -> invalid_arg (sprintf "unrecognized nucleotide %c" nuc)

	let revcomp s =
		let l = String.length s
		let s' = String.create l
		for i = 0 to l-1 do
			s'.[l-i-1] <- comp s.[i]
		s'

let ord a =
	let t = Hashtbl.create 16
	Array.iteri (fun i x -> Hashtbl.add t x i) a
	fun x -> Hashtbl.find t x

module AA = struct
	type t = char
	let dim = 20
	let compare = compare

	let aa2twenty = [|
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;

      -1; 0;-1; 1; 2; 3; 4; 5; 6; 7;-1; 8; 9;10;11;-1;
      12;13;14;15;16;-1;17;18;-1;19;-1;-1;-1;-1;-1;-1;

      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;
      -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1 |]


	let twenty2aa = ord aa2twenty

	let is c =
		let ci = Char.code c
		ci <= (Array.length aa2twenty) && aa2twenty.(ci) >= 0

	let index c = aa2twenty.(Char.code c)

	let of_index ci = Char.chr (twenty2aa ci)

	let blosum62_matrix = [|
	   4; 0;-2;-1;  -2; 0;-2;-1;  -1;-1;-1;-2;  -1;-1;-1;-1;   0; 0;-3;-2;
       0; 9;-3;-4;  -2;-3;-3;-1;  -3;-1;-1;-3;  -3;-3;-3;-1;  -1;-1;-2;-2;
       -2;-3; 6; 2;  -3;-1;-1;-3;  -1;-4;-3; 1;  -1; 0;-2; 0;  -1;-3;-4;-3;
       -1;-4; 2; 5;  -3;-2; 0;-3;   1;-3;-2; 0;  -1; 2; 0; 0;  -1;-2;-3;-2;

       -2;-2;-3;-3;   6;-3;-1; 0;  -3; 0; 0;-3;  -4;-3;-3;-2;  -2;-1; 1; 3;
       0;-3;-1;-2;  -3; 6;-2;-4;  -2;-4;-3; 0;  -2;-2;-2; 0;  -2;-3;-2;-3;
       -2;-3;-1; 0;  -1;-2; 8;-3;  -1;-3;-2; 1;  -2; 0; 0;-1;  -2;-3;-2; 2;
       -1;-1;-3;-3;   0;-4;-3; 4;  -3; 2; 1;-3;  -3;-3;-3;-2;  -1; 3;-3;-1;

       -1;-3;-1; 1;  -3;-2;-1;-3;   5;-2;-1; 0;  -1; 1; 2; 0;  -1;-2;-3;-2;
       -1;-1;-4;-3;   0;-4;-3; 2;  -2; 4; 2;-3;  -3;-1;-1;-1;  -1; 1;-2;-1;
       -1;-1;-3;-2;   0;-3;-2; 1;  -1; 2; 5;-2;  -2; 0;-1;-1;  -1; 1;-1;-1;
       -2;-3; 1; 0;  -3; 0; 1;-3;   0;-3;-2; 6;  -2; 0; 0; 1;   0;-3;-4;-2;

       -1;-3;-1;-1;  -4;-2;-2;-3;  -1;-3;-2;-2;   7;-1;-2;-1;  -1;-2;-4;-3;
       -1;-3; 0; 2;  -3;-2; 0;-3;   1;-1; 0; 0;  -1; 5; 1; 0;  -1;-2;-2;-1;
       -1;-3;-2; 0;  -3;-2; 0;-3;   2;-1;-1; 0;  -2; 1; 5;-1;  -1;-3;-3;-2;
       -1;-1; 0; 0;  -2; 0;-1;-2;   0;-1;-1; 1;  -1; 0;-1; 4;   1;-2;-3;-2;

       0;-1;-1;-1;  -2;-2;-2;-1;  -1;-1;-1; 0;  -1;-1;-1; 1;   5; 0;-2;-2;
       0;-1;-3;-2;  -1;-3;-3; 3;  -2; 1; 1;-3;  -2;-2;-3;-2;   0; 4;-3;-1;
       -3;-2;-4;-3;   1;-2;-2;-3;  -3;-2;-1;-4;  -4;-2;-3;-3;  -2;-3;11; 2;
       -2;-2;-3;-2;   3;-3; 2;-1;  -2;-1;-1;-2;  -3;-1;-2;-2;  -2;-1; 2 ;7 |]

	let blosum62 aa1 aa2 =
		let idx1 = aa2twenty.(Char.code aa1)
		let idx2 = aa2twenty.(Char.code aa2)
		blosum62_matrix.(20*idx1+idx2)


module Codon64 = struct
	type t = DNA.t*DNA.t*DNA.t
	let dim = 64
	let compare = compare

	let is (nuc1,nuc2,nuc3) = DNA.is nuc1 && DNA.is nuc2 && DNA.is nuc3

	let index (nuc1,nuc2,nuc3) =
		16 * (DNA.index nuc1) + 4 * (DNA.index nuc2) + (DNA.index nuc3)

	let of_index idx =
		let sixteens = idx / 16
		let fours = (idx - 16 * sixteens) / 4
		let ones = (idx - 16 * sixteens - 4 * fours)
		(DNA.of_index sixteens,
    	 DNA.of_index fours,
    	 DNA.of_index ones)

	let stop1 = index ('T','A','A')
	let stop2 = index ('T','A','G')
	let stop3 = index ('T','G','A')

	let is_stop_index ci = ci = stop1 || ci = stop2 || ci = stop3

	let is_stop c = is_stop_index (index c)

	let translation_table = [|
     'K';'N';'K';'N';
	 'T';'T';'T';'T';
	 'R';'S';'R';'S';
	 'I';'I';'M';'I';

	 'Q';'H';'Q';'H';
	 'P';'P';'P';'P';
	 'R';'R';'R';'R';
	 'L';'L';'L';'L';

	 'E';'D';'E';'D';
	 'A';'A';'A';'A';
	 'G';'G';'G';'G';
	 'V';'V';'V';'V';

	 '*';'Y';'*';'Y';
	 'S';'S';'S';'S';
	 '*';'C';'W';'C';
	 'L';'F';'L';'F' |]

	let translate c = translation_table.(index c)
	let translate_index idx = translation_table.(idx)
module Codon1 = Codon64

module Codon61 = struct
	type t = DNA.t*DNA.t*DNA.t
	let dim = 61

	let is c = Codon1.is c && not (Codon1.is_stop c)

	let order = [| ('T','T','T'); ('T','T','C'); ('T','T','A'); ('T','T','G'); ('T','C','T'); ('T','C','C'); ('T','C','A'); ('T','C','G'); ('T','A','T'); ('T','A','C'); ('T','G','T'); ('T','G','C'); ('T','G','G'); ('C','T','T'); ('C','T','C'); ('C','T','A'); ('C','T','G'); ('C','C','T'); ('C','C','C'); ('C','C','A'); ('C','C','G'); ('C','A','T'); ('C','A','C'); ('C','A','A'); ('C','A','G'); ('C','G','T'); ('C','G','C'); ('C','G','A'); ('C','G','G'); ('A','T','T'); ('A','T','C'); ('A','T','A'); ('A','T','G'); ('A','C','T'); ('A','C','C'); ('A','C','A'); ('A','C','G'); ('A','A','T'); ('A','A','C'); ('A','A','A'); ('A','A','G'); ('A','G','T'); ('A','G','C'); ('A','G','A'); ('A','G','G'); ('G','T','T'); ('G','T','C'); ('G','T','A'); ('G','T','G'); ('G','C','T'); ('G','C','C'); ('G','C','A'); ('G','C','G'); ('G','A','T'); ('G','A','C'); ('G','A','A'); ('G','A','G'); ('G','G','T'); ('G','G','C'); ('G','G','A'); ('G','G','G') |]
	let lookup = ord order

	let index (n1,n2,n3) = lookup (Char.uppercase n1,Char.uppercase n2,Char.uppercase n3)
	let of_index idx = order.(idx)

	let compare c1 c2 = compare (index c1) (index c2)

	let translate c = Codon1.translate c
	let translate_index idx = translate (of_index idx)
module Codon2 = Codon61
