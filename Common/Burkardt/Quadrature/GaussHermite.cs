using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class GaussHermite
    {
        public static void hermite_compute(int order, ref double[] xtab, ref double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
        //
        //  Discussion:
        //
        //    The abscissas are the zeros of the N-th order Hermite polynomial.
        //
        //    The integration interval is ( -oo, +oo ).
        //
        //    The weight function is w(x) = exp ( - x*x ).
        //
        //    The integral to approximate:
        //
        //      Integral ( -oo < X < +oo ) exp ( - X*X ) * F(X) dX
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
        //    30 April 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the formula to be computed.
        //
        //    Output, double XTAB[ORDER], the Gauss-Hermite abscissas.
        //
        //    Output, double WEIGHT[ORDER], the Gauss-Hermite weights.
        //
        {
            double cc;
            double dp2 = 0;
            int i;
            double p1 = 0;
            double s;
            double temp;
            double x = 0;

            cc = 1.7724538509 * typeMethods.r8_gamma((double) (order))
                 / Math.Pow(2.0, order - 1);

            s = Math.Pow(2.0 * (double) (order) + 1.0, 1.0 / 6.0);

            for (i = 0; i < (order + 1) / 2; i++)
            {
                if (i == 0)
                {
                    x = s * s * s - 1.85575 / s;
                }
                else if (i == 1)
                {
                    x = x - 1.14 * Math.Pow((double) (order), 0.426) / x;
                }
                else if (i == 2)
                {
                    x = 1.86 * x - 0.86 * xtab[0];
                }
                else if (i == 3)
                {
                    x = 1.91 * x - 0.91 * xtab[1];
                }
                else
                {
                    x = 2.0 * x - xtab[i - 2];
                }

                PolynomialNS.Hermite.hermite_root(ref x, order, ref dp2, ref p1);

                xtab[i] = x;
                weight[i] = (cc / dp2) / p1;

                xtab[order - i - 1] = -x;
                weight[order - i - 1] = weight[i];
            }

            //
            //  Reverse the order of the XTAB values.
            //
            typeMethods.r8vec_reverse(order, ref xtab);

            if (false)
            {
                for (i = 0; i < order / 2; i++)
                {
                    temp = xtab[i];
                    xtab[i] = xtab[order - 1 - i];
                    xtab[order - 1 - i] = temp;
                }
            }
        }
    }
}