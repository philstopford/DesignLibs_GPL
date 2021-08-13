using System;
using Burkardt.Types;

namespace Burkardt.Function
{
    public static class Cardinal
    {
        public static double[] cardinal_cos(int j, int m, int n, double[] t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CARDINAL_COS evaluates the J-th cardinal cosine basis function.
            //
            //  Discussion:
            //
            //    The base points are T(I) = pi * I / ( M + 1 ), 0 <= I <= M + 1.
            //    Basis function J is 1 at T(J), and 0 at T(I) for I /= J
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    John Boyd,
            //    Exponentially convergent Fourier-Chebyshev quadrature schemes on
            //    bounded and infinite intervals,
            //    Journal of Scientific Computing,
            //    Volume 2, Number 2, 1987, pages 99-109.
            //
            //  Parameters:
            //
            //    Input, int J, the index of the basis function.
            //    0 <= J <= M + 1.
            //
            //    Input, int M, indicates the size of the basis set.
            //
            //    Input, int N, the number of sample points.
            //
            //    Input, double T[N], one or more points in [0,pi] where the
            //    basis function is to be evaluated.
            //
            //    Output, double CARDINAL_COS[N], the value of the function at T.
            //
        {
            double[] c;
            double cj;
            int i;
            const double r8_eps = 2.220446049250313E-016;
            const double r8_pi = 3.141592653589793;
            double tj;

            c = new double[n];

            if ((j % (m + 1)) == 0)
            {
                cj = 2.0;
            }
            else
            {
                cj = 1.0;
            }

            tj = r8_pi * (double)(j) / (double)(m + 1);

            for (i = 0; i < n; i++)
            {
                if (Math.Abs(t[i] - tj) <= r8_eps)
                {
                    c[i] = 1.0;
                }
                else
                {
                    c[i] = typeMethods.r8_mop(j + 1)
                           * Math.Sin(t[i])
                           * Math.Sin((double)(m + 1) * t[i])
                           / cj
                           / (double)(m + 1)
                           / (Math.Cos(t[i]) - Math.Cos(tj));
                }
            }

            return c;
        }

        public static double[] cardinal_sin(int j, int m, int n, double[] t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CARDINAL_SIN evaluates the J-th cardinal sine basis function.
            //
            //  Discussion:
            //
            //    The base points are T(I) = pi * I / ( M + 1 ), 0 <= I <= M + 1.
            //    Basis function J is 1 at T(J), and 0 at T(I) for I /= J
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    John Boyd,
            //    Exponentially convergent Fourier-Chebyshev quadrature schemes on
            //    bounded and infinite intervals,
            //    Journal of Scientific Computing,
            //    Volume 2, Number 2, 1987, pages 99-109.
            //
            //  Parameters:
            //
            //    Input, int J, the index of the basis function.
            //    0 <= J <= M + 1.
            //
            //    Input, int M, indicates the size of the basis set.
            //
            //    Input, int N, the number of sample points.
            //
            //    Input, double T[N], one or more points in [0,pi] where the
            //    basis function is to be evaluated.
            //
            //    Output, double CARDINAL_SIN[N], the value of the function at T.
            //
        {
            int i;
            const double r8_eps = 2.220446049250313E-016;
            const double r8_pi = 3.141592653589793;
            double[] s;
            double tj;

            s = new double[n];

            tj = r8_pi * (double)(j) / (double)(m + 1);

            for (i = 0; i < n; i++)
            {
                if (Math.Abs(t[i] - tj) <= r8_eps)
                {
                    s[i] = 1.0;
                }
                else
                {
                    s[i] = typeMethods.r8_mop(j + 1)
                           * Math.Sin(tj)
                           * Math.Sin((double)(m + 1) * t[i])
                           / (double)(m + 1)
                           / (Math.Cos(t[i]) - Math.Cos(tj));
                }
            }

            return s;
        }
    }
}