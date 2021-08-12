namespace Burkardt.SolveNS
{
    public static class Tridisolve
    {
        public static double[] tridisolve(int n, double[] a, double[] b, double[] c, double[] d )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIDISOLVE solves a tridiagonal system of linear equations.
        //
        //  Discussion:
        //
        //    We can describe an NxN tridiagonal matrix by vectors A, B, and C, where
        //    A and C are of length N-1.  In that case, a linear system can be
        //    represented as
        //                        b(1) * x(1) + c(1) * x(2)   = d(1),
        //      a(j-1) * x(j-1) + b(j) * x(j) + c(j) * x(j+1) = d(j), j = 2:n-1,
        //      a(n-1) * x(n-1) + b(n) * x(n)                 = d(n)
        //
        //    This function produces the solution vector X.
        //
        //    This function is derived from Cleve Moler's Matlab suite.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input, double A(N-1), B(N), C(N-1), the matrix entries.
        //
        //    Input, double D(N), the right hand side.
        //
        //    Output, double TRIDISOLVE[N], the solution.
        //
        {
            double[] bi;
            int j;
            double mu;
            double[] x;

            x = new double[n];

            for (j = 0; j < n; j++)
            {
                x[j] = d[j];
            }

            bi = new double[n];
            for (j = 0; j < n; j++)
            {
                bi[j] = 1.0 / b[j];
            }

            for (j = 0; j < n - 1; j++)
            {
                mu = a[j] * bi[j];
                b[j + 1] = b[j + 1] - mu * c[j];
                x[j + 1] = x[j + 1] - mu * x[j];
            }

            x[n - 1] = x[n - 1] * bi[n - 1];
            for (j = n - 2; 0 <= j; j--)
            {
                x[j] = (x[j] - c[j] * x[j + 1]) * bi[j];
            }

            return x;
        }
    }
}