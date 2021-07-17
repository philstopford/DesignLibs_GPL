using System;
using Burkardt.Types;

namespace Burkardt.Legendre
{
    public static class QuadratureRule
    {
        public static void legendre_compute(int order, ref double[] xtab, ref double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_COMPUTE computes a Gauss-Legendre quadrature rule.
        //
        //  Discussion:
        //
        //    The integration interval is [ -1, 1 ].
        //
        //    The weight function is w(x) = 1.0.
        //
        //    The integral to approximate:
        //
        //      Integral ( -1 <= X <= 1 ) F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 April 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //    ORDER must be greater than 0.
        //
        //    Output, double XTAB[ORDER], the abscissas of the rule.
        //
        //    Output, double WEIGHT[ORDER], the weights of the rule.
        //    The weights are positive, symmetric, and should sum to 2.
        //
        {
            double d1;
            double d2pn;
            double d3pn;
            double d4pn;
            double dp;
            double dpn;
            double e1;
            double fx;
            double h;
            int i;
            int iback;
            int k;
            int m;
            int mp1mi;
            int ncopy;
            int nmove;
            double p;
            double pk;
            double pkm1;
            double pkp1;
            const double r8_pi = 3.1415926535897932385;
            double t;
            double u;
            double v;
            double x0;
            double xtemp;

            if (order < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_COMPUTE - Fatal error!");
                Console.WriteLine("  Illegal value of ORDER = " + order + "");
                return;
            }

            e1 = (double) (order * (order + 1));

            m = (order + 1) / 2;

            for (i = 1; i <= m; i++)
            {
                mp1mi = m + 1 - i;

                t = (double) (4 * i - 1) * r8_pi / (double) (4 * order + 2);

                x0 = Math.Cos(t) * (1.0 - (1.0 - 1.0 / (double) (order))
                    / (double) (8 * order * order));

                pkm1 = 1.0;
                pk = x0;

                for (k = 2; k <= order; k++)
                {
                    pkp1 = 2.0 * x0 * pk - pkm1 - (x0 * pk - pkm1) / (double) (k);
                    pkm1 = pk;
                    pk = pkp1;
                }

                d1 = (double) (order) * (pkm1 - x0 * pk);

                dpn = d1 / (1.0 - x0 * x0);

                d2pn = (2.0 * x0 * dpn - e1 * pk) / (1.0 - x0 * x0);

                d3pn = (4.0 * x0 * d2pn + (2.0 - e1) * dpn) / (1.0 - x0 * x0);

                d4pn = (6.0 * x0 * d3pn + (6.0 - e1) * d2pn) / (1.0 - x0 * x0);

                u = pk / dpn;
                v = d2pn / dpn;
                //
                //  Initial approximation H:
                //
                h = -u * (1.0 + 0.5 * u * (v + u * (v * v - d3pn / (3.0 * dpn))));
                //
                //  Refine H using one step of Newton's method:
                //
                p = pk + h * (dpn + 0.5 * h * (d2pn + h / 3.0
                    * (d3pn + 0.25 * h * d4pn)));

                dp = dpn + h * (d2pn + 0.5 * h * (d3pn + h * d4pn / 3.0));

                h = h - p / dp;

                xtemp = x0 + h;

                xtab[mp1mi - 1] = xtemp;

                fx = d1 - h * e1 * (pk + 0.5 * h * (dpn + h / 3.0
                    * (d2pn + 0.25 * h * (d3pn + 0.2 * h * d4pn))));

                weight[mp1mi - 1] = 2.0 * (1.0 - xtemp * xtemp) / (fx * fx);
            }

            if ((order % 2) == 1)
            {
                xtab[0] = 0.0;
            }

            //
            //  Shift the data up.
            //
            nmove = (order + 1) / 2;
            ncopy = order - nmove;

            for (i = 1; i <= nmove; i++)
            {
                iback = order + 1 - i;
                xtab[iback - 1] = xtab[iback - ncopy - 1];
                weight[iback - 1] = weight[iback - ncopy - 1];
            }

            //
            //  Reflect values for the negative abscissas.
            //
            for (i = 1; i <= order - nmove; i++)
            {
                xtab[i - 1] = -xtab[order - i];
                weight[i - 1] = weight[order - i];
            }

            return;
        }

        public static void legendre_ek_compute(int n, ref double[] x, ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_EK_COMPUTE: Legendre quadrature rule by the Elhay-Kautsky method.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 April 2011
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Sylvan Elhay, Jaroslav Kautsky,
            //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
            //    Interpolatory Quadrature,
            //    ACM Transactions on Mathematical Software,
            //    Volume 13, Number 4, December 1987, pages 399-415.
            //
            //  Parameters:
            //
            //    Input, int N, the order.
            //
            //    Output, double X[N], the abscissas.
            //
            //    Output, double W[N], the weights.
            //
        {
            double[] bj;
            int i;
            double zemu;
            //
            //  Define the zero-th moment.
            //
            zemu = 2.0;
            //
            //  Define the Jacobi matrix.
            //
            bj = new double[n];

            for (i = 0; i < n; i++)
            {
                bj[i] = (double) ((i + 1) * (i + 1))
                        / (double) (4 * (i + 1) * (i + 1) - 1);
                bj[i] = Math.Sqrt(bj[i]);
            }

            for (i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            w[0] = Math.Sqrt(zemu);

            for (i = 1; i < n; i++)
            {
                w[i] = 0.0;
            }

            //
            //  Diagonalize the Jacobi matrix.
            //
            IMTQLX.imtqlx(n, ref x, ref bj, ref w);

            for (i = 0; i < n; i++)
            {
                w[i] = w[i] * w[i];
            }

            return;
        }
    }
}