namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8vec_convolution(int m, double[] x, int n, double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CONVOLUTION returns the convolution of two R8VEC's.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The I-th entry of the convolution can be formed by summing the products 
        //    that lie along the I-th diagonal of the following table:
        //
        //    Y3 | 3   4   5   6   7
        //    Y2 | 2   3   4   5   6
        //    Y1 | 1   2   3   4   5
        //       +------------------
        //        X1  X2  X3  X4  X5
        //
        //    which will result in:
        //
        //    Z = ( X1 * Y1,
        //          X1 * Y2 + X2 * Y1,
        //          X1 * Y3 + X2 * Y2 + X3 * Y1,
        //                    X2 * Y3 + X3 * Y2 + X4 * Y1,
        //                              X3 * Y3 + X4 * Y2 + X5 * Y1,
        //                                        X4 * Y3 + X5 * Y2,
        //                                                  X5 * Y3 )
        //            
        //  Example:
        //
        //    Input:
        //
        //      X = (/ 1, 2, 3, 4 /)
        //      Y = (/ -1, 5, 3 /)
        //
        //    Output:
        //
        //      Z = (/ -1, 3, 10, 17, 29, 12 /)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of X.
        //
        //    Input, double X[M], the first vector to be convolved.
        //
        //    Input, int N, the dimension of Y.
        //
        //    Input, double Y[N], the second vector to be convolved.
        //
        //    Output, double R8VEC_CONVOLUTION[M+N-1], the convolution of X and Y.
        //
        {
            int i;
            int j;
            double[] z;

            z = new double[m + n - 1];

            for (i = 0; i < m + n - 1; i++)
            {
                z[i] = 0.0;
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    z[j + i] = z[j + i] + x[i] * y[j];
                }
            }

            return z;
        }

        public static double[] r8vec_convolution_circ(int n, double[] x, double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CONVOLUTION_CIRC returns the discrete circular convolution of two R8VEC's.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    z(1+m) = xCCy(m) = sum ( 0 <= k <= n-1 ) x(1+k) * y(1+m-k)
        //
        //    Here, if the index of Y becomes nonpositive, it is "wrapped around"
        //    by having N added to it.
        //
        //    The circular convolution is equivalent to multiplication of Y by a
        //    circulant matrix formed from the vector X.
        //
        //  Example:
        //
        //    Input:
        //
        //      X = (/ 1, 2, 3, 4 /)
        //      Y = (/ 1, 2, 4, 8 /)
        //
        //    Output:
        //
        //      Circulant form:
        //
        //      Z = ( 1 4 3 2 )   ( 1 )
        //          ( 2 1 4 3 )   ( 2 )
        //          ( 3 2 1 4 ) * ( 4 )
        //          ( 4 3 2 1 )   ( 8 )
        //
        //      The formula:
        //
        //      Z = (/ 1*1 + 2*8 + 3*4 + 4*2,
        //             1*2 + 2*1 + 3*8 + 4*4,
        //             1*4 + 2*2 + 3*1 + 4*8,
        //             1*8 + 2*4 + 3*2 + 4*1 /)
        //
        //      Result:
        //
        //      Z = (/ 37, 44, 43, 26 /)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vectors.
        //
        //    Input, double X[N], Y[N], the vectors to be convolved.
        //
        //    Output, double R8VEC_CONVOLVE_CIRC[N], the circular convolution of X and Y.
        //
        {
            int i;
            int m;
            double[] z;

            z = new double[n];

            for (m = 1; m <= n; m++)
            {
                z[m - 1] = 0.0;
                for (i = 1; i <= m; i++)
                {
                    z[m - 1] = z[m - 1] + x[i - 1] * y[m - i];
                }

                for (i = m + 1; i <= n; i++)
                {
                    z[m - 1] = z[m - 1] + x[i - 1] * y[n + m - i];
                }
            }

            return z;
        }

    }
}