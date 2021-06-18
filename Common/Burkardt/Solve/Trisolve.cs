using System;

namespace Burkardt
{
    public static class Trisolve
    {
        public static double[] trisolve(int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRISOLVE factors and solves a tridiagonal system.
        //
        //  Discussion:
        //
        //    The three nonzero diagonals of the N by N matrix are stored as 3
        //    columns of an N by 3 matrix.
        //
        //  Example:
        //
        //    Here is how a tridiagonal matrix of order 5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //      A43 A44 A45
        //      A54 A55  *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input/output, double A[N*3].
        //    On input, the tridiagonal matrix.
        //    On output, the data in these vectors has been overwritten
        //    by factorization information.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double TRISOLVE[N], the solution of the linear system.
        //
        {
            int i;
            double[] x;
            double xmult;
            //
            //  The diagonal entries can't be zero.
            //
            for (i = 0; i < n; i++)
            {
                if (a[i + 1 * n] == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("TRISOLVE - Fatal error!");
                    Console.WriteLine("  A(" + i + ",2) = 0.");
                    return null;
                }
            }

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = b[i];
            }

            for (i = 1; i < n; i++)
            {
                xmult = a[i + 0 * n] / a[i - 1 + 1 * n];
                a[i + 1 * n] = a[i + 1 * n] - xmult * a[i - 1 + 2 * n];
                x[i] = x[i] - xmult * x[i - 1];
            }

            x[n - 1] = x[n - 1] / a[n - 1 + 1 * n];
            for (i = n - 2; 0 <= i; i--)
            {
                x[i] = (x[i] - a[i + 2 * n] * x[i + 1]) / a[i + 1 * n];
            }

            return x;
        }
    }
}