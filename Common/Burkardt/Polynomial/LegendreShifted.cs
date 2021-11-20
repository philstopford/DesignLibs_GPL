﻿namespace Burkardt.PolynomialNS;

public static class LegendreShifted
{
    public static double[] p01_polynomial_value(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P01_POLYNOMIAL_VALUE evaluates the shifted Legendre polynomials P01(n,x).
        //
        //  Discussion:
        //
        //    The shifted Legendre polynomial P01(n,x) has the domain [0,1], and
        //    is related to the standard Legendre polynomial P(n,x) by
        //      P01(n,x) = P(n,2*x-1).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 March 2016
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
        //    Daniel Zwillinger, editor,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996.
        //
        //  Parameters:
        //
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest order polynomial to evaluate.
        //    Note that polynomials 0 through N will be evaluated.
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Output, double P01_POLYNOMIAL_VALUE[M*(N+1)], the values of the Legendre
        //    polynomials of order 0 through N.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] v = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            v[i + 0 * m] = 1.0;
        }

        switch (n)
        {
            case < 1:
                return v;
        }

        for (i = 0; i < m; i++)
        {
            v[i + 1 * m] = 2.0 * x[i] - 1.0;
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                v[i + j * m] = ((2 * j - 1) * (2.0 * x[i] - 1.0) * v[i + (j - 1) * m]
                                - (j - 1) * v[i + (j - 2) * m])
                               / j;
            }
        }

        return v;
    }

    public static void p01_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P01_POLYNOMIAL_VALUES: values of shifted Legendre polynomials.
        //
        //  Discussion:
        //
        //    If we denote the Legendre polynomial by P(n)(x), and the shifted 
        //    Legendre polynomial by P01(n)(x), then
        //
        //      P01(n)(x) = P(n)(2*x-1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &N, the order of the function.
        //
        //    Output, double &X, the point where the function is evaluated.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 22;

        double[] fx_vec =
        {
            0.1000000000000000E+01,
            0.2500000000000000E+00,
            -0.4062500000000000E+00,
            -0.3359375000000000E+00,
            0.1577148437500000E+00,
            0.3397216796875000E+00,
            0.2427673339843750E-01,
            -0.2799186706542969E+00,
            -0.1524540185928345E+00,
            0.1768244206905365E+00,
            0.2212002165615559E+00,
            0.0000000000000000E+00,
            -0.1475000000000000E+00,
            -0.2800000000000000E+00,
            -0.3825000000000000E+00,
            -0.4400000000000000E+00,
            -0.4375000000000000E+00,
            -0.3600000000000000E+00,
            -0.1925000000000000E+00,
            0.8000000000000000E-01,
            0.4725000000000000E+00,
            0.1000000000000000E+01
        };

        int[] n_vec =
        {
            0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3, 3,
            3
        };

        double[] x_vec =
        {
            0.625E+00,
            0.625E+00,
            0.625E+00,
            0.625E+00,
            0.625E+00,
            0.625E+00,
            0.625E+00,
            0.625E+00,
            0.625E+00,
            0.625E+00,
            0.625E+00,
            0.50E+00,
            0.55E+00,
            0.60E+00,
            0.65E+00,
            0.70E+00,
            0.75E+00,
            0.80E+00,
            0.85E+00,
            0.90E+00,
            0.95E+00,
            1.00E+00
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
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            n = n_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }
}