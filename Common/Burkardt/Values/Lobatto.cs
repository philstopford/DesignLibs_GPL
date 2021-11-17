namespace Burkardt.Values;

public static class Lobatto
{

    public static void lobatto_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOBATTO_POLYNOMIAL_VALUES returns values of the completed Lobatto polynomials.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      n * LegendreP [ n - 1, x ] - n * x * LegendreP [ n, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0
        //    before the first call.  On each call, the routine increments N_DATA by 1,
        //    and returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref int N, the order of the function.
        //
        //    Output, ref double X, the point where the function is evaluated.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 31;

        double[] fx_vec =
        {
            0.9375000000000000,
            0.7031250000000000,
            -0.9667968750000000,
            -1.501464843750000,
            0.3639221191406250,
            2.001914978027344,
            0.6597948074340820,
            -1.934441328048706,
            -1.769941113889217,
            1.215243665501475,
            0.000000000000000,
            0.8692500000000000,
            1.188000000000000,
            1.109250000000000,
            0.7680000000000000,
            0.2812500000000000,
            -0.2520000000000000,
            -0.7507500000000000,
            -1.152000000000000,
            -1.410750000000000,
            -1.500000000000000,
            -1.410750000000000,
            -1.152000000000000,
            -0.7507500000000000,
            -0.2520000000000000,
            0.2812500000000000,
            0.7680000000000000,
            1.109250000000000,
            1.188000000000000,
            0.8692500000000000,
            0.000000000000000
        };

        int[] n_vec =
        {
            1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3
        };

        double[] x_vec =
        {
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            -1.00,
            -0.90,
            -0.80,
            -0.70,
            -0.60,
            -0.50,
            -0.40,
            -0.30,
            -0.20,
            -0.10,
            0.00,
            0.10,
            0.20,
            0.30,
            0.40,
            0.50,
            0.60,
            0.70,
            0.80,
            0.90,
            1.00
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

    public static void lobatto_polynomial_derivatives(ref int n_data, ref int n, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOBATTO_POLYNOMIAL_DERIVATIVES: derivatives of the completed Lobatto polynomials.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      n * LegendreP [ n - 1, x ] - n * x * LegendreP [ n, x ]
        //
        //     In Mathematica, the completed Lobatto polynomial can be evaluated by:
        //
        //       n * LegendreP [ n - 1, x ] - n * x * LegendreP [ n, x ]
        //
        //     The derivative is:
        //
        //         n * D[LegendreP [ n - 1, x ], {x} ] 
        //       - n * LegendreP [ n, x ] 
        //       - n * x * D[LegendreP [ n, x ], {x}]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0
        //    before the first call.  On each call, the routine increments N_DATA by 1,
        //    and returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref int N, the order of the function.
        //
        //    Output, ref double X, the point where the function is evaluated.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 31;

        double[] fx_vec =
        {
            -0.5,
            2.437500000000000,
            4.031250000000000,
            -3.154296875000000,
            -10.19165039062500,
            -1.019622802734375,
            15.67544555664063,
            10.97668933868408,
            -15.91419786214828,
            -24.33202382177114,
            12.00000000000000,
            5.670000000000000,
            0.9600000000000000,
            -2.310000000000000,
            -4.320000000000000,
            -5.250000000000000,
            -5.280000000000000,
            -4.590000000000000,
            -3.360000000000000,
            -1.770000000000000,
            0.0,
            1.770000000000000,
            3.360000000000000,
            4.590000000000000,
            5.280000000000000,
            5.250000000000000,
            4.320000000000000,
            2.310000000000000,
            -0.9600000000000000,
            -5.670000000000000,
            -12.00000000000000
        };

        int[] n_vec =
        {
            1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3, 3,
            3, 3
        };

        double[] x_vec =
        {
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            -1.00,
            -0.90,
            -0.80,
            -0.70,
            -0.60,
            -0.50,
            -0.40,
            -0.30,
            -0.20,
            -0.10,
            0.00,
            0.10,
            0.20,
            0.30,
            0.40,
            0.50,
            0.60,
            0.70,
            0.80,
            0.90,
            1.00
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