namespace Burkardt.Elliptic
{
    public class FM
    {
        public static double evaluate(double m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_FM evaluates the complete elliptic integral F(M).
            //
            //  Discussion:
            //
            //    The value is computed using Carlson elliptic integrals:
            //
            //      F(m) = RF ( 0, 1-m, 1 ).
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
            //    Input, double M, the argument.
            //
            //    Output, double ELLIPTIC_FM, the function value.
            //
        {
            double errtol;
            int ierr = 0;
            double value;
            double x;
            double y;
            double z;

            x = 0.0;
            y = 1.0 - m;
            z = 1.0;
            errtol = 1.0E-03;

            value = Integral.rf(x, y, z, errtol, ref ierr);

            return value;
        }

        public static void values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_FM_VALUES returns values of the complete elliptic integral F(M).
        //
        //  Discussion:
        //
        //    This is one form of what is sometimes called the complete elliptic
        //    integral of the first kind.
        //
        //    The function is defined by the formula:
        //
        //      F(M) = integral ( 0 <= T <= PI/2 )
        //        dT / sqrt ( 1 - M * sin ( T )^2 )
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      EllipticK[m]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 August 2004
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
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                1.570796326794897E+00,
                1.591003453790792E+00,
                1.612441348720219E+00,
                1.635256732264580E+00,
                1.659623598610528E+00,
                1.685750354812596E+00,
                1.713889448178791E+00,
                1.744350597225613E+00,
                1.777519371491253E+00,
                1.813883936816983E+00,
                1.854074677301372E+00,
                1.898924910271554E+00,
                1.949567749806026E+00,
                2.007598398424376E+00,
                2.075363135292469E+00,
                2.156515647499643E+00,
                2.257205326820854E+00,
                2.389016486325580E+00,
                2.578092113348173E+00,
                2.908337248444552E+00
            }
            ;

            double[] x_vec =
            {
                0.00E+00,
                0.05E+00,
                0.10E+00,
                0.15E+00,
                0.20E+00,
                0.25E+00,
                0.30E+00,
                0.35E+00,
                0.40E+00,
                0.45E+00,
                0.50E+00,
                0.55E+00,
                0.60E+00,
                0.65E+00,
                0.70E+00,
                0.75E+00,
                0.80E+00,
                0.85E+00,
                0.90E+00,
                0.95E+00
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