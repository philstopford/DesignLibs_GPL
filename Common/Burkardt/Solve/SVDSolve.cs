using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class SVDSolve
    {
        public static double[] svd_solve(int m, int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SVD_SOLVE solves a linear system in the least squares sense.
        //
        //  Discussion:
        //
        //    The vector X returned by this routine should always minimize the 
        //    Euclidean norm of the residual ||A*x-b||.
        //
        //    If the matrix A does not have full column rank, then there are multiple
        //    vectors that attain the minimum residual.  In that case, the vector
        //    X returned by this routine is the unique such minimizer that has the 
        //    the minimum possible Euclidean norm, that is, ||A*x-b|| and ||x||
        //    are both minimized.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of A.
        //
        //    Input, int N, the number of columns of A.
        //
        //    Input, double A[M*N], the matrix.
        //
        //    Input, double B[M], the right hand side.
        //
        //    Output, double SVD_SOLVE[N], the least squares solution.
        //
        {
            double[] a_copy;
            double[] e;
            int i;
            int info;
            int lda;
            int ldu;
            int ldv;
            int job;
            double[] sdiag;
            double smax;
            double stol;
            double[] sub;
            double[] u;
            double[] ub;
            double[] v;
            double[] work;
            double[] x;
            //
            //  Get the SVD.
            //
            a_copy = typeMethods.r8mat_copy_new(m, n, a);
            lda = m;
            sdiag = new double[Math.Max(m + 1, n)];
            e = new double[Math.Max(m + 1, n)];
            u = new double[m * m];
            ldu = m;
            v = new double[n * n];
            ldv = n;
            work = new double[m];
            job = 11;

            info = DSVDC.dsvdc(ref a_copy, lda, m, n, ref sdiag, ref e, ref u, ldu, ref v, ldv, work, job);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("SVD_SOLVE - Failure!");
                Console.WriteLine("  The SVD could not be calculated.");
                Console.WriteLine("  LINPACK routine DSVDC returned a nonzero");
                Console.WriteLine("  value of the error flag, INFO = " + info + "");
                return null;
            }

            ub = typeMethods.r8mat_mtv_new(m, m, u, b);
            //
            //  For singular problems, there may be tiny but nonzero singular values
            //  that should be ignored.  This is a reasonable attempt to avoid such 
            //  problems, although in general, the user might wish to control the tolerance.
            //
            smax = typeMethods.r8vec_max(n, sdiag);
            if (smax <= double.Epsilon)
            {
                smax = 1.0;
            }

            stol = double.Epsilon * smax;

            sub = new double[n];

            for (i = 0; i < n; i++)
            {
                sub[i] = 0.0;
                if (i < m)
                {
                    if (stol <= sdiag[i])
                    {
                        sub[i] = ub[i] / sdiag[i];
                    }
                }
            }

            x = typeMethods.r8mat_mv_new(n, n, v, sub);

            return x;
        }
    }
}