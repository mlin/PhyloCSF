type t =
	| Val of float
	| Var of int
	| Add of t*t
	| Sub of t*t
	| Mul of t*t
	| Div of t*t

let rec eval' v = function
	| Val x -> x
	| Var (k) -> v.(k)
	| Add (l,r) -> eval' v l +. eval' v r
	| Sub (l,r) -> eval' v l -. eval' v r
	| Mul (l,r) -> eval' v l *. eval' v r
	| Div (t,b) -> eval' v t /. eval' v b
let eval e v = eval' v e

let rec deriv' k = function
	| Val _ -> Val 0.
	| Var (k') when k = k' -> Val 1.
	| Var (k') -> Val 0.
	| Add (l,r) -> Add ((deriv' k l),(deriv' k r))
	| Sub (l,r) -> Sub ((deriv' k l),(deriv' k r))
	| Mul (l,r) -> Add ((Mul ((deriv' k l),r)), (Mul (l,(deriv' k r))))
	| Div (t,b) -> Div ((Sub ((Mul ((deriv' k t),b)), (Mul (t,(deriv' k b))))), (Mul (b,b)))
let deriv e k = deriv' k e

let rec fixed_point f init =
	let rslt = f init
	if rslt = init then
		rslt
	else
		fixed_point f rslt

let rec simplify_step = function
	| Add (Val 0.,x) | Add (x,Val 0.)
	| Sub (x,Val 0.)
	| Mul (Val 1.,x) | Mul(x,Val 1.)
	| Div (x,Val 1.) -> simplify_step x

	| Sub (e,e') when e=e' -> Val 0.
	| Mul (Val 0.,x) | Mul(x,Val 0.) -> Val 0.

	| Div (e,e') when e=e' -> Val 1.

	| Add (Val l,Val r) -> Val (l +. r)
	| Sub (Val l,Val r) -> Val (l -. r)
	| Mul (Val l,Val r) -> Val (l *. r)
	| Div (Val l,Val r) -> Val (l /. r)

	| Add (l,r) -> Add (simplify_step l, simplify_step r)
	| Sub (l,r) -> Sub (simplify_step l, simplify_step r)
	| Mul (l,r) -> Mul (simplify_step l, simplify_step r)
	| Div (l,r) -> Div (simplify_step l, simplify_step r)

	| x -> x

let simplify = fixed_point simplify_step

let to_string ?(fmt=string_of_float) x =
	let b = Buffer.create 10
	let rec f = function
		| Val x -> Buffer.add_string b (fmt x)
		| Var k ->
			Buffer.add_char b 'v'
			Buffer.add_string b (string_of_int k)
		| Add (l,r) ->
			Buffer.add_char b '('
			f l
			Buffer.add_string b " + "
			f r
			Buffer.add_char b ')'
		| Sub (l,r) ->
			Buffer.add_char b '('
			f l
			Buffer.add_string b " - "
			f r
			Buffer.add_char b ')'
		| Mul (l,r) ->
			Buffer.add_char b '('
			f l
			Buffer.add_string b " * "
			f r
			Buffer.add_char b ')'
		| Div (l,r) ->
			Buffer.add_char b '('
			f l
			Buffer.add_string b " / "
			f r
			Buffer.add_char b ')'
	f x
	Buffer.contents b
