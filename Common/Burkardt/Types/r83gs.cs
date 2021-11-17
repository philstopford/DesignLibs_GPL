using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r83_gs_sl(int n, double[] a, double[] b, ref double[] x, int it_max,
            int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_GS_SL solves a R83 system using Gauss-Seidel iteration.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:N), the diagonal in
        //    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
        //    original matrix is "collapsed" vertically into the array.
        //
        //    This routine simply applies a given number of steps of the
        //    iteration to an input approximate solution.  On first call, you can
        //    simply pass in the zero vector as an approximate solution.  If
        //    the returned value is not acceptable, you may call again, using
        //    it as the starting point for additional iterations.
        //
        //  Example:
        //
        //    Here is how a R83 matrix of order 5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 2.
        //
        //    Input, double A[3*N], the R83 matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Input/output, double X[N], an approximate solution to the system.
        //
        //    Input, int IT_MAX, the maximum number of iterations to take.
        //
        //    Input, int JOB, specifies the system to solve.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
    {
        int i;
        int it_num;
        //
        //  No diagonal matrix entry can be zero.
        //
        for (i = 0; i < n; i++)
        {
            switch (a[1 + i * 3])
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R83_GS_SL - Fatal error!");
                    Console.WriteLine("  Zero diagonal entry, index = " + i + "");
                    return;
            }
        }

        switch (job)
        {
            case 0:
            {
                for (it_num = 1; it_num <= it_max; it_num++)
                {
                    x[0] = (b[0] - a[2 + 0 * 3] * x[1]) / a[1 + 0 * 3];
                    for (i = 1; i < n - 1; i++)
                    {
                        x[i] = (b[i] - a[0 + i * 3] * x[i - 1] - a[2 + i * 3] * x[i + 1]) / a[1 + i * 3];
                    }

                    x[n - 1] = (b[n - 1] - a[0 + (n - 1) * 3] * x[n - 2]) / a[1 + (n - 1) * 3];
                }

                break;
            }
            default:
            {
                for (it_num = 1; it_num <= it_max; it_num++)
                {
                    x[0] = (b[0] - a[0 + 1 * 3] * x[1])
                           / a[1 + 0 * 3];
                    for (i = 1; i < n - 1; i++)
                    {
                        x[i] = (b[i] - a[2 + (i - 1) * 3] * x[i - 1] - a[0 + (i + 1) * 3] * x[i + 1])
                               / a[1 + i * 3];
                    }

                    x[n - 1] = (b[n - 1] - a[2 + (n - 2) * 3] * x[n - 2])
                               / a[1 + (n - 1) * 3];
                }

                break;
            }
        }
    }
}