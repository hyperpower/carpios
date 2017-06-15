/************************
 //  \file   Qsort.h
 //  \brief
 // 
 //  \author zhou
 //  \date   12 mai 2014 
 ***********************/
#ifndef _SORT_H_
#define _SORT_H_

#include "array_list.hpp"

namespace carpio
{
template<typename TYPE>
inline void _swap_(TYPE& a, TYPE& b)
{
	TYPE temp = a;
	a = b;
	b = temp;
}

template<typename TYPE>
inline bool _compare_less_(TYPE& a, TYPE& b)
{
	return a < b;
}

template<typename TYPE>
inline bool _compare_great_(TYPE& a, TYPE& b)
{
	return a > b;
}

template<typename TYPE>
void Sort_Bubble(ArrayListT<TYPE> & arr,
		bool (*compare)(const TYPE&, const TYPE&))
{
	int num = arr.size();
	int i, j;
	for (i = 0; i < num; i++) {
		for (j = 1; j < num - i; j++) {
			if (compare(arr[j - 1], arr[j]))
				_swap_(arr[j - 1], arr[j]);
		}
	}
}

template<typename TYPE>
void Sort_Insert(ArrayListT<TYPE> & arr,
		bool (*compare)(const TYPE&, const TYPE&))
{
	int num = arr.size();
	TYPE temp;
	int i, j;
	for (i = 1; i < num; i++) {
		temp = arr[i];
		for (j = i; j > 0 && compare(arr[j - 1], temp); j--)
			arr[j] = arr[j - 1];
		arr[j] = temp;
	}
}

template<typename TYPE>
int QSort(ArrayListV<TYPE>& v, int base_ptr, int total_elems,
		bool (*compare)(const TYPE&, const TYPE&))
{
	double pivot_buffer;
	const int MAX_THRESH = 4;
	const int BYTES_PER_WORD = 8;
	typedef struct
	{
		int lo;
		int hi;
	} stack_node;

	if (total_elems > MAX_THRESH) {

		int lo = base_ptr;
		int hi = lo + total_elems - 1;

		stack_node stack[BYTES_PER_WORD * sizeof(long)]; /* Largest size needed for 32-bit int!!! */
		stack_node *top = stack + 1;

		while (stack < top) {
			int left_ptr;
			int right_ptr;
			{
				{
					/* Select median value from among LO, MID, and HI. Rearrange
					 LO and HI so the three values are sorted. This lowers the
					 probability of picking a pathological pivot value and
					 skips a comparison for both the LEFT_PTR and RIGHT_PTR. */

					int mid = lo + (hi - lo) / 2;

					if (compare(v[mid], v[lo]))
						_swap_(v[mid], v[lo]);
					if (compare(v[hi], v[mid]))
						_swap_(v[hi], v[mid]);
					else
						goto jump_over;

					if (compare(v[mid], v[lo]))
						_swap_(v[mid], v[lo]);

					jump_over:

					pivot_buffer = v[mid];
				}

				left_ptr = lo + 1;
				right_ptr = hi - 1;

				/* Here's the famous ``collapse the walls'' section of quicksort.
				 Gotta like those tight inner loops!  They are the main reason
				 that this algorithm runs much faster than others. */
				do {
					while (compare(v[left_ptr], pivot_buffer))
						left_ptr++;

					while (compare(pivot_buffer, v[right_ptr]))
						right_ptr--;

					if (left_ptr < right_ptr) {
						_swap_(v[left_ptr], v[right_ptr]);
						left_ptr++;
						right_ptr--;
					} else if (left_ptr == right_ptr) {
						left_ptr++;
						right_ptr--;
						break;
					}
				} while (left_ptr <= right_ptr);
			}

			/* Set up pointers for next iteration.  First determine whether
			 left and right partitions are below the threshold size. If so,
			 ignore one or both.  Otherwise, push the larger partition's
			 bounds on the stack and continue sorting the smaller one. */
#define _push(LOW,HIGH) do {top->lo = LOW;top++->hi = HIGH;} while (0)
#define _pop(LOW,HIGH)  do {LOW = (--top)->lo;HIGH = top->hi;} while (0)

			if ((right_ptr - lo) <= MAX_THRESH) {
				if ((hi - left_ptr) <= MAX_THRESH)
					_pop(lo, hi);
				else
					lo = left_ptr;
			} else if ((hi - left_ptr) <= MAX_THRESH)
				hi = right_ptr;
			else if ((right_ptr - lo) > (hi - left_ptr)) {
				_push(lo, right_ptr);
				lo = left_ptr;
			} else {
				_push(left_ptr, hi);
				hi = right_ptr;
			}
		}
	}

	/* Once the BASE_PTR array is partially sorted by quicksort the rest
	 is completely sorted using insertion sort, since this is efficient
	 for partitions below MAX_THRESH size. BASE_PTR points to the beginning
	 of the array to sort, and END_PTR points at the very last element in
	 the array (*not* one beyond it!). */

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

	{
		int end_ptr = base_ptr + total_elems - 1;
		int run_ptr;
		int tmp_ptr = base_ptr;
		int thresh = MIN(end_ptr, base_ptr + MAX_THRESH);

		for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++)
			if (compare(v[run_ptr], v[tmp_ptr]))
				tmp_ptr = run_ptr;

		if (tmp_ptr != base_ptr)
			_swap_(v[tmp_ptr], v[base_ptr]);

		for (run_ptr = base_ptr + 1; (tmp_ptr = run_ptr += 1) <= end_ptr;) {

			while (compare(v[run_ptr], v[tmp_ptr -= 1]))
				;

			if ((tmp_ptr += 1) != run_ptr) {
				int trav;

				for (trav = run_ptr + 1; --trav >= run_ptr;) {
					double c;
					c = v[trav];
					int hi, lo;

					for (hi = lo = trav; (lo -= 1) >= tmp_ptr; hi = lo)
						v[hi] = v[lo];
					v[hi] = c;
				}
			}
		}
	}
	return 1;
}

template<typename TYPE>
int QSort_ascend(ArrayListV<TYPE>& v, int base_ptr, int total_elems)
{
	return QSort(v, base_ptr, total_elems, _compare_less_);
}

template<typename TYPE>
int QSort_decrease(ArrayListV<TYPE>& v, int base_ptr, int total_elems)
{
	return QSort(v, base_ptr, total_elems, _compare_great_);
}

template<typename TYPE>
int QSort(ArrayListV<St>& v,
		ArrayListV<TYPE>& x,
		int base_ptr,
		int total_elems,
		bool (*compare)(const TYPE&, const TYPE&))
{
	int pivot_buffer;
	//double pixot_buffer;
	const int MAX_THRESH = 4;
	const int BYTES_PER_WORD = 8;
	typedef struct
	{
		int lo;
		int hi;
	} stack_node;
	if (total_elems > MAX_THRESH) {

		int lo = base_ptr;
		int hi = lo + total_elems - 1;

		stack_node stack[BYTES_PER_WORD * sizeof(long)]; /* Largest size needed for 32-bit int!!! */
		stack_node *top = stack + 1;

		while (stack < top) {
			int left_ptr;
			int right_ptr;
			{
				{
					/* Select median value from among LO, MID, and HI. Rearrange
					 LO and HI so the three values are sorted. This lowers the
					 probability of picking a pathological pivot value and
					 skips a comparison for both the LEFT_PTR and RIGHT_PTR. */

					int mid = lo + (hi - lo) / 2;

					if (compare(v[mid], v[lo])) {
						_swap_(v[mid], v[lo]);
						_swap_(x[mid], x[lo]);
					}
					if (compare(v[hi], v[mid])) {
						_swap_(v[hi], v[mid]);
						_swap_(x[hi], x[mid]);
					} else
						goto jump_over;

					if (compare(v[mid], v[lo])) {
						_swap_(v[mid], v[lo]);
						_swap_(x[mid], x[lo]);
					}

					jump_over:

					pivot_buffer = v[mid];
					//pixot_buffer = x[mid];
				}

				left_ptr = lo + 1;
				right_ptr = hi - 1;

				/* Here's the famous ``collapse the walls'' section of quicksort.
				 Gotta like those tight inner loops!  They are the main reason
				 that this algorithm runs much faster than others. */
				do {
					while (compare(v[left_ptr], pivot_buffer))
						left_ptr++;

					while (compare(pivot_buffer, v[right_ptr]))
						right_ptr--;

					if (left_ptr < right_ptr) {
						_swap_(v[left_ptr], v[right_ptr]);
						_swap_(x[left_ptr], x[right_ptr]);
						left_ptr++;
						right_ptr--;
					} else if (left_ptr == right_ptr) {
						left_ptr++;
						right_ptr--;
						break;
					}
				} while (left_ptr <= right_ptr);
			}

			/* Set up pointers for next iteration.  First determine whether
			 left and right partitions are below the threshold size. If so,
			 ignore one or both.  Otherwise, push the larger partition's
			 bounds on the stack and continue sorting the smaller one. */
#define _push(LOW,HIGH) do {top->lo = LOW;top++->hi = HIGH;} while (0)
#define _pop(LOW,HIGH)  do {LOW = (--top)->lo;HIGH = top->hi;} while (0)

			if ((right_ptr - lo) <= MAX_THRESH) {
				if ((hi - left_ptr) <= MAX_THRESH)
					_pop(lo, hi);
				else
					lo = left_ptr;
			} else if ((hi - left_ptr) <= MAX_THRESH)
				hi = right_ptr;
			else if ((right_ptr - lo) > (hi - left_ptr)) {
				_push(lo, right_ptr);
				lo = left_ptr;
			} else {
				_push(left_ptr, hi);
				hi = right_ptr;
			}
		}
	}

	/* Once the BASE_PTR array is partially sorted by quicksort the rest
	 is completely sorted using insertion sort, since this is efficient
	 for partitions below MAX_THRESH size. BASE_PTR points to the beginning
	 of the array to sort, and END_PTR points at the very last element in
	 the array (*not* one beyond it!). */

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

	{
		int end_ptr = base_ptr + total_elems - 1;
		int run_ptr;
		int tmp_ptr = base_ptr;
		int thresh = MIN(end_ptr, base_ptr + MAX_THRESH);

		for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++)
			if (compare(v[run_ptr], v[tmp_ptr]))
				tmp_ptr = run_ptr;

		if (tmp_ptr != base_ptr) {
			_swap_(v[tmp_ptr], v[base_ptr]);
			_swap_(x[tmp_ptr], x[base_ptr]);
		}

		for (run_ptr = base_ptr + 1; (tmp_ptr = run_ptr += 1) <= end_ptr;) {

			while (compare(v[run_ptr], v[tmp_ptr -= 1]))
				;

			if ((tmp_ptr += 1) != run_ptr) {
				int trav;

				for (trav = run_ptr + 1; --trav >= run_ptr;) {
					int c;
					double d;
					c = v[trav];
					d = x[trav];
					int hi, lo;

					for (hi = lo = trav; (lo -= 1) >= tmp_ptr; hi = lo) {
						v[hi] = v[lo];
						x[hi] = x[lo];
					}
					v[hi] = c;
					x[hi] = d;
				}
			}
		}
	}
	return 1;
}

template<typename TYPE>
int QSort_ascend(ArrayListV<St>& v, ArrayListV<TYPE>& x, int base_ptr,
		int total_elems)
{
	return QSort(v, x, base_ptr, total_elems, _compare_less_);
}

}

#endif /* QSORT_H_ */
