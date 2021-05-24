using System;
using System.Linq;
using Burkardt.Table;

namespace Burkardt.Latinizer
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
			double[] v = new double[n];

			for (int i = 0; i < m; i++ )
			{
				for (int j = 0; j < n; j++ )
				{
					v[j] = table[i+j*m];
				}
				double v_min = v.Min();
				double v_max = v.Max();
				int[] indx = TableMisc.r8vec_sort_heap_index_a_new ( n, v );

				for (int j = 0; j < n; j++ )
				{
					table[i+indx[j]*m] =  ( ( double ) ( n - j - 1 ) * v_min   
											+ ( double ) (     j     ) * v_max ) 
										  / ( double ) ( n     - 1 );

				}
			}

			return table;
		}        
		
	}
}