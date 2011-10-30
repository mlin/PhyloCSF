(** low-level routines for numerical optimization *)

(** {1 Single-dimensional maximization } *)

(** a [maximizer] iteratively narrows down a local maximum of a one-dimensional function [f(x)] *)
class type maximizer = object
	(** current estimate of the locally maximizing setting, [argmax_x f(x)] *)
	method maximum : unit -> float

	(** current lower and upper bounds on [x] around the local maximum *)
	method interval : unit -> (float*float)

	(** perform the next iteration *)
	method iterate : unit -> unit

(** create a [maximizer] for the function [f(x)]
@param init initial estimate of [x]. It is required that [f(lo) < f(init) > f(hi)].
@param lo lower bound on [x] for the search
@param hi upper bound on [x] for the search
*)
val make_maximizer : f:(float -> float) -> init:float -> lo:float -> hi:float -> maximizer

(** perform a random search within an interval [(lo,hi)] for a suitable initialization point for the one-dimensional maximizer, i.e. a point [x] such that [f(lo) < f(x) > f(hi)].
	@param init any initial guess between [lo] and [hi] (does not have to satisfy any other conditions, but a good guess would make this go faster)
	@param maxtries maximum number of random points to try
	@param logspace perform the random search over the interval between [lo] and [hi] in log space, with the result that the portion of the interval close to [lo] is searched more thoroughly than the portion close to [hi]. [lo] must be positive.
	@raise Out_of_range if a suitable starting point could not be found. The parameter to the exception is [lo] if [f lo > f hi], which suggests that the maximum is at a point less than [lo], or [hi] if [f hi > f lo], which suggests that the maximum is at a point greater than [hi].
	@raise Failure if a suitable starting point could not be found and [f lo = f hi] in fewer than [maxtries] attempts.
*)
val find_init : ?maxtries:int -> ?logspace:bool -> f:(float -> float) -> init:float -> lo:float -> hi:float -> unit -> float
exception Out_of_range of float

(** {1 Multi-dimensional gradient ascent } *)

(** a multi_maximizer performs gradient ascent on a multidimensional function *)
class type multi_maximizer = object
	(** current estimate of the maximum [f(x)], the corresponding location [x], and the gradient [df(x)] *)
	method maximum : unit -> float*(float array)*(float array)

	(** perform the next iteration *)
	method iterate : unit -> unit

(** create a maximizer for the function [f(x)] given its gradient [df]
	@param init initial point from which to begin the gardient ascent
*)
val make_multi_maximizer : f:(float array -> float) -> df:(float array -> float array) -> init:(float array) -> multi_maximizer

type domain = Real | Pos | Neg | NonPos | NonNeg | Probability

(** [check_domain domain x] returns [true] if [x] is in [domain], false otherwise *)
val check_domain : domain -> float -> bool

(** wrap the function [f(x)] to prevent the [multi_maximizer] from searching outside of the specified domain for each element of the vector input. *)
val enforce_domain : ?rejection_value:float -> (float array -> float) -> (domain array) -> (float array -> float)
