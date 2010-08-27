#include <stddef.h>
#include <string.h>
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/alloc.h>
#include <caml/fail.h>
#include <caml/callback.h>
#include <math.h>
#include <xmmintrin.h>

value dot(value a, value b, value len)
{
	int size = Int_val(len);
	double z = 0.0, fres = 0.0;
	double ftmp[2] __attribute__ ((aligned (16))) = { 0.0, 0.0 };
	__m128d mres;
	double *pa = (double*) a;
	double *pb = (double*) b;

	if ((size / 2) != 0) {
		unsigned int i;
		unsigned int idx = 0;
		mres = _mm_load_sd(&z);
		for (i = 0; i < size / 2; i++) {
			mres = _mm_add_pd(mres, _mm_mul_pd(_mm_loadu_pd(&pa[idx]),
			_mm_loadu_pd(&pb[idx])));
			idx += 2;
		}

		_mm_store_pd(ftmp, mres);                

		fres = ftmp[0] + ftmp[1];
	}

	if ((size % 2) != 0) {
		unsigned int i;
		for (i = size - size % 2; i < size; i++)
			fres += pa[i] * pb[i];
	}

	return copy_double(fres);
}
