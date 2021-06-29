using System;

namespace Burkardt.MinpackNS
{
    public static partial class Minpack
    {
        public static void dogleg(int n, double[] r, int lr, double[] diag, double[] qtb,
                double delta, ref double[] x, double[] wa1, double[] wa2, int rIndex = 0, int wa1Index = 0, int xIndex = 0, int wa2Index = 0, int qtbIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    dogleg() combines Gauss-Newton and gradient for a minimizing step.
            //
            //  Discussion:
            //
            //    Given an M by N matrix A, an n by n nonsingular diagonal
            //    matrix d, an m-vector b, and a positive number delta, the
            //    problem is to determine the convex combination x of the
            //    gauss-newton and scaled gradient directions that minimizes
            //    (a*x - b) in the least squares sense, subject to the
            //    restriction that the euclidean norm of d*x be at most delta.
            //
            //    This function completes the solution of the problem
            //    if it is provided with the necessary information from the
            //    qr factorization of a. 
            //
            //    That is, if a = q*r, where q has orthogonal columns and r is an upper 
            //    triangular matrix, then dogleg expects the full upper triangle of r and
            //    the first n components of Q'*b.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 April 2010
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Jorge More, Burton Garbow, Kenneth Hillstrom,
            //    User Guide for MINPACK-1,
            //    Technical Report ANL-80-74,
            //    Argonne National Laboratory, 1980.
            //
            //  Parameters:
            //
            //    Input, int N, the order of R.
            //
            //    Input, double R[LR], the upper triangular matrix R stored by rows.
            //
            //    Input, int LR, the size of the storage for R, which should be at
            //    least (n*(n+1))/2.
            //
            //    Input, double DIAG[N], the diagonal elements of the matrix D.
            //
            //    Input, double QTB[N], the first n elements of the vector 
            //    (q transpose)*b.
            //
            //    Input, double DELTA, an upper bound on the euclidean norm of d*x.
            //
            //    Output, double X[N], contains the desired convex combination of the 
            //    gauss-newton direction and the scaled gradient direction.
            //
            //    Workspace, WA1[N].
            //
            //    Workspace, WA2[N].
            //
        {
            double alpha;
            double bnorm;
            double epsmch;
            double gnorm;
            int i;
            int j;
            int jj;
            int jp1;
            int k;
            int l;
            double qnorm;
            double sgnorm;
            double sum;
            double temp;
            //
            //  EPSMCH is the machine precision.
            //
            epsmch = double.Epsilon;
            //
            //  Calculate the Gauss-Newton direction.
            //
            jj = (n * (n + 1)) / 2 + 1;

            for (k = 1; k <= n; k++)
            {
                j = n - k + 1;
                jp1 = j + 1;
                jj = jj - k;
                l = jj + 1;
                sum = 0.0;
                for (i = jp1; i <= n; i++)
                {
                    sum = sum + r[rIndex + (l - 1)] * x[xIndex + (i - 1)];
                    l = l + 1;
                }

                temp = r[rIndex + (jj - 1)];
                if (temp == 0.0)
                {
                    l = j;
                    for (i = 1; i <= j; i++)
                    {
                        temp = Math.Max(temp, Math.Abs(r[rIndex + (l - 1)]));
                        l = l + n - i;
                    }

                    temp = epsmch * temp;
                    if (temp == 0.0)
                    {
                        temp = epsmch;
                    }
                }

                x[xIndex + (j - 1)] = (qtb[qtbIndex + (j - 1)] - sum) / temp;
            }

            //
            //  Test whether the Gauss-Newton direction is acceptable.
            //
            for (j = 0; j < n; j++)
            {
                wa1[wa1Index + (j)] = 0.0;
                wa2[wa2Index + (j)] = diag[j] * x[xIndex + (j)];
            }

            qnorm = Helpers.enorm(n, wa2, wa2Index);

            if (qnorm <= delta)
            {
                return;
            }

            //
            //  The Gauss-Newton direction is not acceptable.
            //  Calculate the scaled gradient direction.
            //
            l = 0;
            for (j = 0; j < n; j++)
            {
                temp = qtb[qtbIndex + (j)];
                for (i = j; i < n; i++)
                {
                    wa1[wa1Index + (i - 1)] = wa1[wa1Index + (i - 1)] + r[rIndex + (l - 1)] * temp;
                    l = l + 1;
                }

                wa1[wa1Index + (j)] = wa1[wa1Index + (j)] / diag[j];
            }

            //
            //  Calculate the norm of the scaled gradient and test for
            //  the special case in which the scaled gradient is zero.
            //
            gnorm = Helpers.enorm(n, wa1, xIndex:wa1Index);
            sgnorm = 0.0;
            alpha = delta / qnorm;
            //
            //  Calculate the point along the scaled gradient
            //  at which the quadratic is minimized.
            //
            if (gnorm != 0.0)
            {
                for (j = 0; j < n; j++)
                {
                    wa1[wa1Index + (j)] = (wa1[wa1Index + (j)] / gnorm) / diag[j];
                }

                l = 0;
                for (j = 0; j < n; j++)
                {
                    sum = 0.0;
                    for (i = j; i < n; i++)
                    {
                        sum = sum + r[rIndex + (l)] * wa1[wa1Index + (i)];
                        l = l + 1;
                    }

                    wa2[wa2Index + (j)] = sum;
                }

                temp = Helpers.enorm(n, wa2, wa2Index);
                sgnorm = (gnorm / temp) / temp;
                alpha = 0.0;
                //
                //  If the scaled gradient direction is not acceptable,
                //  calculate the point along the dogleg at which the quadratic is minimized.
                //
                if (sgnorm < delta)
                {
                    bnorm = Helpers.enorm(n, qtb, qtbIndex);
                    temp = (bnorm / gnorm) * (bnorm / qnorm) * (sgnorm / delta);
                    temp = temp - (delta / qnorm) * (sgnorm / delta) * (sgnorm / delta)
                           + Math.Sqrt(Math.Pow(temp - (delta / qnorm), 2)
                                       + (1.0 - (delta / qnorm) * (delta / qnorm))
                                       * (1.0 - (sgnorm / delta) * (sgnorm / delta)));
                    alpha = ((delta / qnorm)
                             * (1.0 - (sgnorm / delta) * (sgnorm / delta))) / temp;
                }
            }

            //
            //  Form appropriate convex combination of the Gauss-Newton
            //  direction and the scaled gradient direction.
            //
            temp = (1.0 - alpha) * Math.Min(sgnorm, delta);
            for (j = 0; j < n; j++)
            {
                x[xIndex + (j)] = temp * wa1[wa1Index + (j)] + alpha * x[xIndex + (j)];
            }
        }
    }
}