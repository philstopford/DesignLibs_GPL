using System;

namespace Burkardt.Elliptic;

public static class FA
{
    public static double evaluate(double a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_FA evaluates the complete elliptic integral F(A).
        //
        //  Discussion:
        //
        //    The value is computed using Carlson elliptic integrals:
        //
        //      F(a) = RF ( 0, 1-sin^2(a), 1 ).
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
        //    Output, double ELLIPTIC_FA, the function value.
        //
    {
        int ierr = 0;

        double k = Math.Sin(a * Math.PI / 180.0);
        const double x = 0.0;
        double y = (1.0 - k) * (1.0 + k);
        const double z = 1.0;
        const double errtol = 1.0E-03;

        double value = Integral.rf(x, y, z, errtol, ref ierr);

        return value;
    }

    public static void values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_FA_VALUES returns values of the complete elliptic integral F(A).
        //
        //  Discussion:
        //
        //    This is one form of what is sometimes called the complete elliptic integral
        //    of the first kind.
        //
        //    The function is defined by the formula:
        //
        //      F(A) = integral ( 0 <= T <= PI/2 )
        //        dT / sqrt ( 1 - sin ( A )^2 * sin ( T )^2 )
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      EllipticK[(Sin[a*Pi/180])^2]
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
        const int N_MAX = 18;

        double[] fx_vec =
            {
                0.1570796326794897E+01,
                0.1573792130924768E+01,
                0.1582842804338351E+01,
                0.1598142002112540E+01,
                0.1620025899124204E+01,
                0.1648995218478530E+01,
                0.1685750354812596E+01,
                0.1731245175657058E+01,
                0.1786769134885021E+01,
                0.1854074677301372E+01,
                0.1935581096004722E+01,
                0.2034715312185791E+01,
                0.2156515647499643E+01,
                0.2308786798167196E+01,
                0.2504550079001634E+01,
                0.2768063145368768E+01,
                0.3153385251887839E+01,
                0.3831741999784146E+01
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

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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