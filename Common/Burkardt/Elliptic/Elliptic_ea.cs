using System;

namespace Burkardt.Elliptic
{
    public static class EA
    {
        public static double evaluate(double a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_EA evaluates the complete elliptic integral E(A).
            //
            //  Discussion:
            //
            //    The value is computed using Carlson elliptic integrals:
            //
            //      E(a) = RF ( 0, 1-sin^2(a), 1 ) - 1/3 sin^2(a) RD ( 0, 1-sin^2(a), 1 ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, the argument.
            //
            //    Output, double ELLIPTIC_EA, the function value.
            //
        {
            double errtol;
            int ierr = 0;
            double k;
            
            double value;
            double x;
            double y;
            double z;

            k = Math.Sin(a * Math.PI / 180.0);

            x = 0.0;
            y = (1.0 - k) * (1.0 + k);
            z = 1.0;
            errtol = 1.0E-03;

            value = Integral.rf(x, y, z, errtol, ref ierr)
                    - k * k * Integral.rd(x, y, z, errtol, ref ierr) / 3.0;

            return value;
        }

        public static void values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_EA_VALUES returns values of the complete elliptic integral E(A).
        //
        //  Discussion:
        //
        //    This is one form of what is sometimes called the complete elliptic
        //    integral of the second kind.
        //
        //    The function is defined by the formula:
        //
        //      E(A) = integral ( 0 <= T <= PI/2 )
        //        sqrt ( 1 - sin ( A )^2 * sin ( T )^2 ) dT
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      EllipticE[(Sin[Pi*a/180])^2]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2004
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
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &X, the argument of the function, measured
        //    in degrees.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 18;

            double[] fx_vec =
            {
                1.570796326794897E+00,
                1.567809073977622E+00,
                1.558887196601596E+00,
                1.544150496914673E+00,
                1.523799205259774E+00,
                1.498114928422116E+00,
                1.467462209339427E+00,
                1.432290969306756E+00,
                1.393140248523812E+00,
                1.350643881047676E+00,
                1.305539094297794E+00,
                1.258679624779997E+00,
                1.211056027568459E+00,
                1.163827964493139E+00,
                1.118377737969864E+00,
                1.076405113076403E+00,
                1.040114395706010E+00,
                1.012663506234396E+00
            }
            ;

            double[] x_vec =
            {
                0.0E+00,
                5.0E+00,
                10.0E+00,
                15.0E+00,
                20.0E+00,
                25.0E+00,
                30.0E+00,
                35.0E+00,
                40.0E+00,
                45.0E+00,
                50.0E+00,
                55.0E+00,
                60.0E+00,
                65.0E+00,
                70.0E+00,
                75.0E+00,
                80.0E+00,
                85.0E+00
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