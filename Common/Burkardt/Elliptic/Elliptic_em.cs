namespace Burkardt.Elliptic;

public class EM
{
    public static double evaluate(double m)

        //****************************************************************************80

        //
        //  Purpose:
        //
        //    ELLIPTIC_EM evaluates the complete elliptic integral E(M).
        //
        //  Discussion:
        //
        //    The value is computed using Carlson elliptic integrals:
        //
        //      E(m) = RF ( 0, 1-m, 1 ) - 1/3 m RD ( 0, 1-m, 1 ).
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
        //    Output, double ELLIPTIC_EM, the function value.
        //
    {
        double errtol;
        int ierr = 0;
        double value = 0;
        double x;
        double y;
        double z;

        x = 0.0;
        y = 1.0 - m;
        z = 1.0;
        errtol = 1.0E-03;

        value = Integral.rf(x, y, z, errtol, ref ierr)
                - m * Integral.rd(x, y, z, errtol, ref ierr) / 3.0;

        return value;
    }

    public static void values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_EM_VALUES returns values of the complete elliptic integral E(M).
        //
        //  Discussion:
        //
        //    This is one form of what is sometimes called the complete elliptic
        //    integral of the second kind.
        //
        //    The function is defined by the formula:
        //
        //      E(M) = integral ( 0 <= T <= PI/2 )
        //        sqrt ( 1 - M * sin ( T )^2 ) dT
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      EllipticE[m]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2004
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
        const int N_MAX = 20;

        double[] fx_vec =
            {
                1.570796326794897E+00,
                1.550973351780472E+00,
                1.530757636897763E+00,
                1.510121832092819E+00,
                1.489035058095853E+00,
                1.467462209339427E+00,
                1.445363064412665E+00,
                1.422691133490879E+00,
                1.399392138897432E+00,
                1.375401971871116E+00,
                1.350643881047676E+00,
                1.325024497958230E+00,
                1.298428035046913E+00,
                1.270707479650149E+00,
                1.241670567945823E+00,
                1.211056027568459E+00,
                1.178489924327839E+00,
                1.143395791883166E+00,
                1.104774732704073E+00,
                1.060473727766278E+00
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