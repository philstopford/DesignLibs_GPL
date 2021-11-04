using System;
using Burkardt.Types;

namespace Burkardt.ChebyshevNS
{
    public static class Chebyshev
    {
        public static double[] cheby(int nf, int npl, Func<double, double[]> functn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY carries out the Chebyshev analysis of one or more functions.
            //
            //  Discussion:
            //
            //    This routine carries out the simultaneous Chebyshev analysis of 
            //    NF functions.
            //
            //    The output is a matrix containing one Chebyshev series per column.
            //
            //    An example of a routine to compute the function values is:
            //
            //      double *functn ( double a )
            //      {
            //        double *val;
            //        val = new double[2];
            //        val[0] = sin(a);
            //        val[1] = cos(a);
            //        return val;
            //      }
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 September 2011
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Roger Broucke.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Roger Broucke,
            //    Algorithm 446:
            //    Ten Subroutines for the Manipulation of Chebyshev Series,
            //    Communications of the ACM,
            //    October 1973, Volume 16, Number 4, pages 254-256.
            //
            //  Parameters:
            //
            //    Input, int NF, the number of functions to be analyzed.
            //
            //    Input, int NPL, the number of terms in the 
            //    Chebyshev series.
            //
            //    Input, int NPLMAX, the leading dimension of X.
            //
            //    Input, external FUNCTN, the name of a routine which computes
            //    the function values at any given point.
            //
            //    Output, double CHEBY[NPL*NF], the Chebyshev series.
            //
        {
            double enn;
            double fk;
            double[] fxj;
            double[] gc;
            int j;
            int k;
            int l;
            int lm;
            int n;
            double pen;
            double[] x;
            double xj;

            x = typeMethods.r8vec_zero_new(npl * nf);

            n = npl - 1;
            enn = (double) (n);
            pen = 3.1415926535897932 / enn;

            gc = new double[2 * n];
            for (k = 1; k <= 2 * n; k++)
            {
                fk = (double) (k - 1);
                gc[k - 1] = Math.Cos(fk * pen);
            }

            for (j = 0; j < npl; j++)
            {
                xj = gc[j];
                fxj = functn(xj);

                if (j == 0 || j == npl - 1)
                {
                    typeMethods.r8vec_scale(0.5, nf, ref fxj);
                }

                for (l = 0; l < npl; l++)
                {
                    lm = (l * j) % (2 * n);
                    for (k = 0; k < nf; k++)
                    {
                        x[l + k * npl] = x[l + k * npl] + fxj[k] * gc[lm];
                    }
                }
            }

            typeMethods.r8vec_scale(2.0 / enn, npl * nf, ref x);

            return x;
        }

        public static double[] chebyshev_even1(int n, Func<double, double> f)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBYSHEV_EVEN1 returns the even Chebyshev coefficients of F.
            //
            //  Discussion:
            //
            //    The coefficients are calculated using the extreme points of Tn(x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    D B Hunter, H V Smith,
            //    A quadrature formula of Clenshaw-Curtis type for the Gegenbauer 
            //    weight function,
            //    Journal of Computational and Applied Mathematics,
            //    Volume 177, 2005, pages 389-400.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points to use.
            //    1 <= N.
            //
            //    Input, double F ( double x ), the function to be 
            //    integrated with the Gegenbauer weight.
            //
            //    Output, double CHEBYSHEV_EVEN1[1+N/2], the even Chebyshev coefficients of F.
            //
        {
            double[] a2;
            int j;
            int r;
            int rh;
            double r8_n;
            
            int s;
            double total;

            s = (n / 2);

            r8_n = (double) (n);

            a2 = new double[s + 1];

            for (r = 0; r <= 2 * s; r = r + 2)
            {
                total = 0.5 * f(1.0);
                for (j = 1; j < n; j++)
                {
                    total = total + f(Math.Cos(j * Math.PI / r8_n))
                        * Math.Cos(r * j * Math.PI / r8_n);
                }

                total = total + 0.5 * typeMethods.r8_mop(r) * f(-1.0);
                rh = r / 2;
                a2[rh] = (2.0 / r8_n) * total;
            }

            return a2;
        }

        public static double[] chebyshev_even2(int n, Func<double, double> f)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBYSHEV_EVEN2 returns the even Chebyshev coefficients of F.
            //
            //  Discussion:
            //
            //    The coefficients are calculated using the zeros of Tn(x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    D B Hunter, H V Smith,
            //    A quadrature formula of Clenshaw-Curtis type for the Gegenbauer 
            //    weight function,
            //    Journal of Computational and Applied Mathematics,
            //    Volume 177, 2005, pages 389-400.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points to use.
            //    1 <= N.
            //
            //    Input, double F ( double x ), the function to be 
            //    integrated with the Gegenbauer weight.
            //
            //    Output, double CHEBYSHEV_EVEN2(0:N/2), the even Chebyshev coefficients of F.
            //
        {
            double[] b2;
            int j;
            int r;
            
            int rh;
            int s;
            double total;
            double x1;
            double x2;

            s = (n / 2);

            b2 = new double[s + 1];

            for (r = 0; r <= 2 * s; r = r + 2)
            {
                total = 0.0;
                for (j = 0; j <= n; j++)
                {
                    x1 = (double) (2 * j + 1) * Math.PI / 2.0 / (double) (n + 1);
                    x2 = (double) (r * (2 * j + 1)) * Math.PI
                         / 2.0 / (double) (n + 1);
                    total = total + f(Math.Cos(x1)) * Math.Cos(x2);
                }

                rh = r / 2;
                b2[rh] = (2.0 / (double) (n + 1)) * total;
            }

            return b2;
        }

        public static double[] chebyshev_coefficients(double a, double b, int n,
                Func<double, double> f)

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
                double[] x)

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