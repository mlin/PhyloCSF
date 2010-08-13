(** symbolic representation of arithmetic expressions with variables *)

type t =
	| Val of float	(** Constant value *)
	| Var of int	(** Variable (identified by an integer) *)
	| Add of t*t	(** Addition *)
	| Sub of t*t	(** Subtraction *)
	| Mul of t*t	(** Multiplication *)
	| Div of t*t	(** Division *)

(** [eval expr assignments] evaluates [expr] with respect to the given assignment of the variables (i.e. [assignments.(i)] is substituted for variable [i]) *)
val eval : t -> float array -> float

(** [deriv expr v] computes the derivative of [expr] with respect to variable [v]. *)
val deriv : t -> int -> t

(** simplify the expression by evaluating operations on constants, removing additions of 0, multiplications by 1, etc.*)
val simplify : t -> t

(** make a nice infix notation string of the expression
@param fmt if you want to control the precision, etc. (defaults to [string_of_float])
 *)
val to_string : ?fmt:(float->string) -> t -> string
