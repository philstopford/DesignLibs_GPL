﻿namespace Burkardt.Probability;

public static class Psi
{
    public static void psi_values(ref int n_data, ref double x, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PSI_VALUES returns some values of the Psi or Digamma function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      PolyGamma[x]
        //
        //    or
        //
        //      Polygamma[0,x]
        //
        //    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
        //
        //    PSI(1) = -Euler's constant.
        //
        //    PSI(X+1) = PSI(X) + 1 / X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 June 2013
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
            -10.42375494041108E+00,
            -5.289039896592188E+00,
            -3.502524222200133E+00,
            -2.561384544585116E+00,
            -1.963510026021423E+00,
            -1.540619213893190E+00,
            -1.220023553697935E+00,
            -0.9650085667061385E+00,
            -0.7549269499470514E+00,
            -0.5772156649015329E+00,
            -0.4237549404110768E+00,
            -0.2890398965921883E+00,
            -0.1691908888667997E+00,
            -0.6138454458511615E-01,
            0.3648997397857652E-01,
            0.1260474527734763E+00,
            0.2085478748734940E+00,
            0.2849914332938615E+00,
            0.3561841611640597E+00,
            0.4227843350984671E+00
        };

        double[] x_vec =
        {
            0.1E+00,
            0.2E+00,
            0.3E+00,
            0.4E+00,
            0.5E+00,
            0.6E+00,
            0.7E+00,
            0.8E+00,
            0.9E+00,
            1.0E+00,
            1.1E+00,
            1.2E+00,
            1.3E+00,
            1.4E+00,
            1.5E+00,
            1.6E+00,
            1.7E+00,
            1.8E+00,
            1.9E+00,
            2.0E+00
        };

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