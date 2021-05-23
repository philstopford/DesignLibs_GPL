using System;
using System.Linq;

namespace Latinizer
{
	public static class Latinize
	{
		//****************************************************************************80
		public static double[] r8mat_latinize ( int m, int n, double[] table )
		//****************************************************************************80
		//
		//  Purpose:
		//
		//    R8MAT_LATINIZE "latinizes" a real table.
		//
		//  Discussion:
		//
		//    On output, each row of the table will have the properties that:
		//    1) the minimum and maximum row values are the same as on input;
		//    2) the row contains N evenly spaced values between the
		//       minimum and maximum;
		//    3) in each row, the elements retain their ordering.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    04 February 2012
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, int M, the spatial dimension.
		//
		//    Input, int N, the number of columns.
		//
		//    Input/output, double TABLE[M*N].  On input, the dataset to
		//    be "latinized".  On output, the latinized dataset.
		//
		{
			int i;
			int[] indx;
			int j;

			double[] v = new double[n];

			for ( i = 0; i < m; i++ )
			{
				for ( j = 0; j < n; j++ )
				{
					v[j] = table[i+j*m];
				}
				double v_min = v.Min();
				double v_max = v.Max();
				indx = r8vec_sort_heap_index_a_new ( n, v );

				for ( j = 0; j < n; j++ )
				{
					table[i+indx[j]*m] =  ( ( double ) ( n - j - 1 ) * v_min   
											+ ( double ) (     j     ) * v_max ) 
										  / ( double ) ( n     - 1 );

				}
			}

			return table;
		}        
		
		static int[] r8vec_sort_heap_index_a_new ( int n, double[] a )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    R8VEC_SORT_HEAP_INDEX_A_NEW does an indexed heap ascending sort of an R8VEC.
		//
		//  Discussion:
		//
		//    The sorting is not actually carried out.  Rather an index array is
		//    created which defines the sorting.  This array may be used to sort
		//    or index the array, or to sort or index related arrays keyed on the
		//    original array.
		//
		//    Once the index array is computed, the sorting can be carried out
		//    "implicitly:
		//
		//      A(INDX(I)), I = 1 to N is sorted,
		//
		//    after which A(I), I = 1 to N is sorted.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    30 March 2004
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, int N, the number of entries in the array.
		//
		//    Input, double A[N], an array to be index-sorted.
		//
		//    Output, int R8VEC_SORT_HEAP_INDEX_A_NEW[N], contains the sort index.  The
		//    I-th element of the sorted array is A(INDX(I)).
		//
		{
			int i;

			int[] indx = new int[n];

			for (i = 1; i <= n; i++ )
			{
				indx[i-1] = i;
			}

			int l = n / 2 + 1;
			int ir = n;

			for ( ; ; )
			{
				double aval;
				int indxt;
				if ( 1 < l )
				{
					l = l - 1;
					indxt = indx[l-1];
					aval = a[indxt-1];
				}
				else
				{
					indxt = indx[ir-1];
					aval = a[indxt-1];
					indx[ir-1] = indx[0];
					ir = ir - 1;

					if ( ir == 1 )
					{
						indx[0] = indxt;
						for ( i = 0; i < n; i++ )
						{
							indx[i] = indx[i] - 1;
						}
						break;
					}
				}

				i = l;
				int j = l + l;

				while ( j <= ir )
				{
					if ( j < ir )
					{
						if ( a[indx[j-1]-1] < a[indx[j]-1] )
						{
							j = j + 1;
						}
					}

					if ( aval < a[indx[j-1]-1] )
					{
						indx[i-1] = indx[j-1];
						i = j;
						j = j + j;
					}
					else
					{
						j = ir + 1;
					}
				}
				indx[i-1] = indxt;
			}
			return indx;
		}
	}
}