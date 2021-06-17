using System;

namespace Burkardt
{
    public static partial class Matrix
    {
        public static double[] kahan(double alpha, int m, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KAHAN returns the KAHAN matrix.
            //
            //  Formula:
            //
            //    if ( I = J )
            //      A(I,I) =  sin(ALPHA)^(I)
            //    else if ( I < J )
            //      A(I,J) = - sin(ALPHA)^(I) * cos(ALPHA)
            //    else
            //      A(I,J) = 0
            //
            //  Example:
            //
            //    ALPHA = 0.25, N = 4
            //
            //    S  -C*S    -C*S      -C*S
            //    0     S^2  -C*S^2    -C*S^2
            //    0     0       S^3    -C*S^3
            //    0     0       0         S^4
            //
            //    where
            //
            //      S = sin(ALPHA), C=COS(ALPHA)
            //
            //  Properties:
            //
            //    A is upper triangular.
            //
            //    A = B * C, where B is a diagonal matrix and C is unit upper triangular.
            //    For instance, for the case M = 3, N = 4:
            //
            //    A = | S 0    0    |  * | 1 -C -C  -C |
            //        | 0 S^2  0    |    | 0  1 -C  -C |
            //        | 0 0    S^3  |    | 0  0  1  -C |
            //
            //    A is generally not symmetric: A' /= A.
            //
            //    A has some interesting properties regarding estimation of
            //    condition and rank.
            //
            //    det ( A ) = sin(ALPHA^(N*(N+1)/2).
            //
            //    LAMBDA(I) = sin ( ALPHA )^I
            //
            //    A is nonsingular if and only if sin ( ALPHA ) =/= 0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 June 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Nicholas Higham,
            //    A survey of condition number estimation for triangular matrices,
            //    SIAM Review,
            //    Volume 9, 1987, pages 575-596.
            //
            //    W Kahan,
            //    Numerical Linear Algebra,
            //    Canadian Mathematical Bulletin,
            //    Volume 9, 1966, pages 757-801.
            //
            //  Parameters:
            //
            //    Input, double ALPHA, the scalar that defines A.  A typical
            //    value is 1.2.  The "interesting" range of ALPHA is 0 < ALPHA < PI.
            //
            //    Input, int M, N, the order of the matrix.
            //
            //    Output, double KAHAN[M*N], the matrix.
            //
        {
            double[] a;
            double csi;
            int i;
            int j;
            double si;

            a = new double[m * n];

            for (i = 0; i < m; i++)
            {
                si = Math.Pow(Math.Sin(alpha), i + 1);
                csi = -Math.Cos(alpha) * si;
                for (j = 0; j < n; j++)
                {
                    if (j < i)
                    {
                        a[i + j * m] = 0.0;
                    }
                    else if (j == i)
                    {
                        a[i + j * m] = si;
                    }
                    else
                    {
                        a[i + j * m] = csi;
                    }
                }
            }

            return a;
        }

        public static double[] kahan_inverse(double alpha, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KAHAN_INVERSE returns the inverse of the KAHAN matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 June 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double ALPHA, the scalar that defines A.  A typical 
            //    value is 1.2.  The "interesting" range of ALPHA is 0 < ALPHA < PI.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, double KAHAN_INVERSE[N*N], the matrix.
            //
        {
            double[] a;
            double ci;
            int i;
            int j;
            double si;

            a = new double[n * n];

            ci = Math.Cos(alpha);

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    if (i == j)
                    {
                        a[i + j * n] = 1.0;
                    }
                    else if (i == j - 1)
                    {
                        a[i + j * n] = ci;
                    }
                    else if (i < j)
                    {
                        a[i + j * n] = ci * Math.Pow(1.0 + ci, j - i - 1);
                    }
                    else
                    {
                        a[i + j * n] = 0.0;
                    }
                }
            }

            //
            //  Scale the columns.
            //
            for (j = 0; j < n; j++)
            {
                si = Math.Pow(Math.Sin(alpha), j + 1);
                for (i = 0; i < n; i++)
                {
                    a[i + j * n] = a[i + j * n] / si;
                }
            }

            return a;
        }
    }
}