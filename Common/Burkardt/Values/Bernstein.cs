﻿namespace Burkardt.Values;

public static class Bernstein
{
    public static void bernstein_poly_values(ref int n_data, ref int n, ref int k, ref double x,
        ref double b)
    {
        bernstein_poly_01_values(ref n_data, ref n, ref k, ref x, ref b);
    }
    public static void bernstein_poly_01_values(ref int n_data, ref int n, ref int k, ref double x,
            ref double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_01_VALUES returns some values of the Bernstein polynomials.
        //
        //  Discussion:
        //
        //    The Bernstein polynomials are assumed to be based on [0,1].
        //
        //    The formula for the Bernstein polynomials is
        //
        //      B(N,I)(X) = [N!/(I//(N-I)!)] * (1-X)^(N-I) * X^I
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Binomial[n,i] * (1-x)^(n-i) * x^i
        //
        //  First values:
        //
        //    B(0,0)(X) = 1
        //
        //    B(1,0)(X) =      1-X
        //    B(1,1)(X) =               X
        //
        //    B(2,0)(X) =     (1-X)^2
        //    B(2,1)(X) = 2 * (1-X)   * X
        //    B(2,2)(X) =               X^2
        //
        //    B(3,0)(X) =     (1-X)^3
        //    B(3,1)(X) = 3 * (1-X)^2 * X
        //    B(3,2)(X) = 3 * (1-X)   * X^2
        //    B(3,3)(X) =               X^3
        //
        //    B(4,0)(X) =     (1-X)^4
        //    B(4,1)(X) = 4 * (1-X)^3 * X
        //    B(4,2)(X) = 6 * (1-X)^2 * X^2
        //    B(4,3)(X) = 4 * (1-X)   * X^3
        //    B(4,4)(X) =               X^4
        //
        //  Special values:
        //
        //    B(N,I)(X) has a unique maximum value at X = I/N.
        //
        //    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
        //
        //    B(N,I)(1/2) = C(N,K) / 2^N
        //
        //    For a fixed X and N, the polynomials add up to 1:
        //
        //      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
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
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref int N, the degree of the polynomial.
        //
        //    Output, ref int K, the index of the polynomial.
        //
        //    Output, ref double X, the argument of the polynomial.
        //
        //    Output, ref double B, the value of the polynomial B(N,K)(X).
        //
    {
        const int N_MAX = 15;

        double[] b_vec =
        {
            0.1000000000000000E+01,
            0.7500000000000000E+00,
            0.2500000000000000E+00,
            0.5625000000000000E+00,
            0.3750000000000000E+00,
            0.6250000000000000E-01,
            0.4218750000000000E+00,
            0.4218750000000000E+00,
            0.1406250000000000E+00,
            0.1562500000000000E-01,
            0.3164062500000000E+00,
            0.4218750000000000E+00,
            0.2109375000000000E+00,
            0.4687500000000000E-01,
            0.3906250000000000E-02
        };

        int[] k_vec =
        {
            0,
            0, 1,
            0, 1, 2,
            0, 1, 2, 3,
            0, 1, 2, 3, 4
        };

        int[] n_vec =
        {
            0,
            1, 1,
            2, 2, 2,
            3, 3, 3, 3,
            4, 4, 4, 4, 4
        };

        double[] x_vec =
        {
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00,
            0.25E+00
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
            n = 0;
            k = 0;
            x = 0.0;
            b = 0.0;
        }
        else
        {
            n = n_vec[n_data - 1];
            k = k_vec[n_data - 1];
            x = x_vec[n_data - 1];
            b = b_vec[n_data - 1];
        }
    }

}