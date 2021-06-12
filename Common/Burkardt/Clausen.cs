using System;
using Burkardt.Chebyshev;
using Burkardt.CorrelationNS;
using Burkardt.Types;

namespace Burkardt
{
    public static class Clausen
    {
        public static double clausen(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CLAUSEN evaluates the Clausen function Cl2(x).
            //
            //  Discussion:
            //
            //    Note that the first coefficient, a0 in Koelbig's paper, 
            //    is doubled here, to account for a different convention in
            //    Chebyshev coefficients.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 December 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Kurt Koelbig,
            //    Chebyshev coefficients for the Clausen function Cl2(x),
            //    Journal of Computational and Applied Mathematics,
            //    Volume 64, Number 3, 1995, pages 295-297.
            //
            //  Parameters:
            //
            //    Input, double X, the evaluation point.
            //
            //    Output, double CLAUSEN, the value of the function.
            //
        {
            //
            //  Chebyshev expansion for -pi/2 < x < +pi/2.
            //
            double[] c1 = {
                0.05590566394715132269,
                0.00000000000000000000,
                0.00017630887438981157,
                0.00000000000000000000,
                0.00000126627414611565,
                0.00000000000000000000,
                0.00000001171718181344,
                0.00000000000000000000,
                0.00000000012300641288,
                0.00000000000000000000,
                0.00000000000139527290,
                0.00000000000000000000,
                0.00000000000001669078,
                0.00000000000000000000,
                0.00000000000000020761,
                0.00000000000000000000,
                0.00000000000000000266,
                0.00000000000000000000,
                0.00000000000000000003
            }
            ;

            //
            //  Chebyshev expansion for pi/2 < x < 3pi/2.
            //
            double[] c2 = {
                0.00000000000000000000,
                -0.96070972149008358753,
                0.00000000000000000000,
                0.04393661151911392781,
                0.00000000000000000000,
                0.00078014905905217505,
                0.00000000000000000000,
                0.00002621984893260601,
                0.00000000000000000000,
                0.00000109292497472610,
                0.00000000000000000000,
                0.00000005122618343931,
                0.00000000000000000000,
                0.00000000258863512670,
                0.00000000000000000000,
                0.00000000013787545462,
                0.00000000000000000000,
                0.00000000000763448721,
                0.00000000000000000000,
                0.00000000000043556938,
                0.00000000000000000000,
                0.00000000000002544696,
                0.00000000000000000000,
                0.00000000000000151561,
                0.00000000000000000000,
                0.00000000000000009172,
                0.00000000000000000000,
                0.00000000000000000563,
                0.00000000000000000000,
                0.00000000000000000035,
                0.00000000000000000000,
                0.00000000000000000002
            }
            ;
            const double r8_pi = 3.141592653589793;
            int n1 = 19;
            int n2 = 30;
            double value;
            double x2;
            double x3;
            double xa;
            double xb;
            double xc;
            //
            //  The function is periodic.  Wrap X into [-pi/2, 3pi/2].
            //
            xa = -0.5 * r8_pi;
            xb = 0.5 * r8_pi;
            xc = 1.5 * r8_pi;
            x2 = typeMethods.r8_wrap(x, xa, xc);
            //
            //  Choose the appropriate expansion.
            //
            if (x2 < xb)
            {
                x3 = 2.0 * x2 / r8_pi;
                value = x2 - x2 * Math.Log(Math.Abs(x2))
                        + 0.5 * Math.Pow(x2, 3) * Correlation.r8_csevl(x3, c1, n1);
            }
            else
            {
                x3 = 2.0 * x2 / r8_pi - 2.0;
                value = Correlation.r8_csevl(x3, c2, n2);
            }

            return value;
        }

        public static void clausen_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLAUSEN_VALUES returns some values of the Clausen's integral.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      CLAUSEN(x) = integral ( 0 <= t <= x ) -ln ( 2 * sin ( t / 2 ) ) dt
        //
        //    The data was reported by McLeod.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Allan McLeod,
        //    Algorithm 757:
        //    MISCFUN: A software package to compute uncommon special functions,
        //    ACM Transactions on Mathematical Software,
        //    Volume 22, Number 3, September 1996, pages 288-301.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.14137352886760576684E-01,
                0.13955467081981281934E+00,
                -0.38495732156574238507E+00,
                0.84831187770367927099E+00,
                0.10139591323607685043E+01,
                -0.93921859275409211003E+00,
                0.72714605086327924743E+00,
                0.43359820323553277936E+00,
                -0.98026209391301421161E-01,
                -0.56814394442986978080E+00,
                -0.70969701784448921625E+00,
                0.99282013254695671871E+00,
                -0.98127747477447367875E+00,
                -0.64078266570172320959E+00,
                0.86027963733231192456E+00,
                0.39071647608680211043E+00,
                0.47574793926539191502E+00,
                0.10105014481412878253E+01,
                0.96332089044363075154E+00,
                -0.61782699481929311757E+00
            }
            ;

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                -0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                -1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                -3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                -5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                -10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                -30.0000000000E+00,
                50.0000000000E+00
            }
            ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }
    }
}