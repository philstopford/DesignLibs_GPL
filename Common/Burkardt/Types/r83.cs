using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int r83_np_fa(int n, double[] a)
//****************************************************************************80
//
//  Purpose:
//
//    R83_NP_FA factors a R83 system without pivoting.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    Because this routine does not use pivoting, it can fail even when
//    the matrix is not singular, and it is liable to make larger
//    errors.
//
//    R83_NP_FA and R83_NP_SL may be preferable to the corresponding
//    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
//    in one step, and does not save the factorization.
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
//    11 January 2004
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
//    Input/output, double A[3*N].
//    On input, the tridiagonal matrix.  On output, factorization information.
//
//    Output, int R83_NP_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
        {
            for (int i = 1; i <= n - 1; i++)
            {
                if (a[1 + (i - 1) * 3] == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R83_NP_FA - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + i + "");
                    return i;
                }

//
//  Store the multiplier in L.
//
                a[2 + (i - 1) * 3] = a[2 + (i - 1) * 3] / a[1 + (i - 1) * 3];
//
//  Modify the diagonal entry in the next column.
//
                a[1 + i * 3] = a[1 + i * 3] - a[2 + (i - 1) * 3] * a[0 + i * 3];
            }

            if (a[1 + (n - 1) * 3] == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R83_NP_FA - Fatal error!");
                Console.WriteLine("  Zero pivot on step " + n + "");
                return n;
            }

            return 0;
        }

        public static double[] r83_np_sl(int n, double[] a_lu, double[] b, int job )
//****************************************************************************80
//
//  Purpose:
//
//    R83_NP_SL solves a R83 system factored by R83_NP_FA.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
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
//    12 January 2004
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
//    Input, double A_LU[3*N], the LU factors from R83_NP_FA.
//
//    Input, double B[N], the right hand side of the linear system.
//    On output, B contains the solution of the linear system.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, double R83_NP_SL[N], the solution of the linear system.
//
        {
            double[] x = new double[n];

            for (int i = 0; i < n; i++)
            {
                x[i] = b[i];
            }

            if (job == 0)
            {
//
//  Solve L * Y = B.
//
                for (int i = 1; i < n; i++)
                {
                    x[i] = x[i] - a_lu[2 + (i - 1) * 3] * x[i - 1];
                }

//
//  Solve U * X = Y.
//
                for (int i = n; 1 <= i; i--)
                {
                    x[i - 1] = x[i - 1] / a_lu[1 + (i - 1) * 3];
                    if (1 < i)
                    {
                        x[i - 2] = x[i - 2] - a_lu[0 + (i - 1) * 3] * x[i - 1];
                    }
                }
            }
            else
            {
//
//  Solve U' * Y = B
//
                for (int i = 1; i <= n; i++)
                {
                    x[i - 1] = x[i - 1] / a_lu[1 + (i - 1) * 3];
                    if (i < n)
                    {
                        x[i] = x[i] - a_lu[0 + i * 3] * x[i - 1];
                    }
                }

//
//  Solve L' * X = Y.
//
                for (int i = n - 1; 1 <= i; i--)
                {
                    x[i - 1] = x[i - 1] - a_lu[2 + (i - 1) * 3] * x[i];
                }
            }

            return x;
        }
    }
}