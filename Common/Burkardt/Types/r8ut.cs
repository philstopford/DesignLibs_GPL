using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8ut_inverse ( int n, double[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8UT_INVERSE computes the inverse of a R8UT matrix.
            //
            //  Discussion:
            //
            //    The R8UT storage format is used for an M by N upper triangular matrix,
            //    and allocates space even for the zero entries.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double A[N*N], the R8UT matrix.
            //
            //    Output, double R8UT_INVERSE[N*N], the inverse of the upper
            //    triangular matrix.
            //
        {
            double[] b;
            int i;
            int j;
            int k;
            //
            //  Check.
            //
            for ( i = 0; i < n; i++ )
            {
                if ( a[i+i*n] == 0.0 )
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8UT_INVERSE - Fatal error!");
                    Console.WriteLine("  Zero diagonal element.");
                    return null;
                }
            }

            b = new double[n*n];

            for ( j = n-1; 0 <= j; j-- )
            {
                for ( i = n-1; 0 <= i; i-- )
                {
                    if ( j < i )
                    {
                        b[i+j*n] = 0.0;
                    }
                    else if ( i == j )
                    {
                        b[i+j*n] = 1.0 / a[i+j*n];
                    }
                    else if ( i < j )
                    {
                        b[i+j*n] = 0.0;

                        for ( k = i+1; k <= j; k++ )
                        {
                            b[i+j*n] = b[i+j*n] - a[i+k*n] * b[k+j*n];
                        }
                        b[i+j*n] = b[i+j*n] / a[i+i*n];
                    }
                }
            }

            return b;
        }
    }
}