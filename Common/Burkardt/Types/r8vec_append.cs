using System.Collections.Generic;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_append(ref int n, ref double[][] a, double value)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_APPEND appends an entry to an R8VEC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, int *N, the current size of the array.  On output,
            //    the array is one entry longer.
            //
            //    Input/output, double **A, the array.  On output, the array has had 
            //    VALUE appended.
            //
            //    Input, double VALUE, a value to be appended to A.
            //
        {
            for (int r = 0; r < a.Length; r++)
            {
                double[] temp = new double[n + 1];
                for (int c = 0; c < n; c++)
                {
                    temp[c] = a[r][c];
                }

                temp[n] = value;
                a[r] = temp;
            }
            //
            //  Increase N.
            //
            n = n + 1;
        }

        public static double[] r8vec_append_new(int n, double[] a, double value )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_APPEND_NEW appends a value to an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 May 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the input vector.
        //
        //    Input, double A[N], the vector to be modified.  On output, the vector
        //    has been reallocated, has one more entry than on input, and that last
        //    entry is VALUE.
        //
        //    Input, double VALUE, the value to be appended to the vector.
        //
        //    Output, double R8VEC_APPEND[N+1], a copy of the vector
        //    with one more entry than on input, and that last
        //    entry is VALUE.
        //
        {
            double[] b;
            int i;

            b = new double[n + 1];

            for (i = 0; i < n; i++)
            {
                b[i] = a[i];
            }

            b[n] = value;

            return b;
        }
    }
}