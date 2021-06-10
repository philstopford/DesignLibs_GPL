using System;

namespace Burkardt.ChebyshevNS
{
    public static class Chebyshev
    {
        public static double[] chebyshev_coefficients(double a, double b, int n,
                Func<double,double> f )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBYSHEV_COEFFICIENTS determines Chebyshev interpolation coefficients.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 September 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Roger Broucke,
            //    Algorithm 446:
            //    Ten Subroutines for the Manipulation of Chebyshev Series,
            //    Communications of the ACM,
            //    Volume 16, Number 4, April 1973, pages 254-256.
            //
            //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
            //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
            //    Second Edition,
            //    Cambridge University Press, 1992,
            //    ISBN: 0-521-43064-X,
            //    LC: QA297.N866.
            //
            //  Parameters:
            //
            //    Input, double A, B, the domain of definition.
            //
            //    Input, int N, the order of the interpolant.
            //
            //    Input, double F ( double X ), an external function.
            //
            //    Output, double CHEBYSHEV_COEFFICIENTS[N], the Chebyshev coefficients.
            //
        {
            double[] fx = new double[n];

            for (int i = 0; i < n; i++)
            {
                double angle = (double) (2 * i + 1) * Math.PI / (double) (2 * n);
                double x = Math.Cos(angle);
                x = 0.5 * (a + b) + x * 0.5 * (b - a);
                fx[i] = f(x);
            }

            double[] c = new double[n];

            for (int i = 0; i < n; i++)
            {
                c[i] = 0.0;
                for (int j = 0; j < n; j++)
                {
                    double angle = (double) (i * (2 * j + 1)) * Math.PI / (double) (2 * n);
                    c[i] = c[i] + fx[j] * Math.Cos(angle);
                }
            }

            for (int i = 0; i < n; i++)
            {
                c[i] = 2.0 * c[i] / (double) (n);
            }

            return c;
        }

        public static double[] chebyshev_interpolant(double a, double b, int n, double[] c, int m,
        double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV_INTERPOLANT evaluates a Chebyshev interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    Volume 16, Number 4, April 1973, pages 254-256.
        //
        //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
        //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
        //    Second Edition,
        //    Cambridge University Press, 1992,
        //    ISBN: 0-521-43064-X,
        //    LC: QA297.N866.
        //
        //  Parameters:
        //
        //    Input, double A, B, the domain of definition.
        //
        //    Input, int N, the order of the polynomial.
        //
        //    Input, double C[N], the Chebyshev coefficients.
        //
        //    Input, int M, the number of points.
        //
        //    Input, double X[M], the point at which the polynomial is
        //    to be evaluated.
        //
        //    Output, double CHEBYSHEF_INTERPOLANT[M], the value of the Chebyshev
        //    polynomial at X.
        //
        {
            double[] cf = new double[m];

            for (int j = 0; j < m; j++)
            {
                double dip1 = 0.0;
                double di = 0.0;
                double y = (2.0 * x[j] - a - b) / (b - a);

                for (int i = n - 1; 1 <= i; i--)
                {
                    double dip2 = dip1;
                    dip1 = di;
                    di = 2.0 * y * dip1 - dip2 + c[i];
                }

                cf[j] = y * di - dip1 + 0.5 * c[0];
            }

            return cf;
        }

        public static double[] chebyshev_zeros(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBYSHEV_ZEROS returns zeroes of the Chebyshev polynomial T(N)(X).
            //
            //  Discussion:
            //
            //    We produce the Chebyshev zeros in ascending order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Output, double CHEBYSHEV_ZEROS[N], the zeroes of T(N)(X).
            //
        {
            double[] x = new double[n];

            for (int i = 0; i < n; i++)
            {
                double angle = (double) (2 * (n - i) - 1) * Math.PI / (double) (2 * n);
                x[i] = Math.Cos(angle);
            }

            return x;
        }

    }
}