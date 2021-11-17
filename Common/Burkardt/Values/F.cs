namespace Burkardt.Values;

public static class F
{
    public static void f_cdf_values(ref int n_data, ref int a, ref int b, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_CDF_VALUES returns some values of the F CDF test function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = FRatioDistribution [ dfn, dfd ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 August 2004
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
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref int A, ref int B, the parameters of the function.
        //
        //    Output, ref double X, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 20;

        int[] a_vec =
        {
            1,
            1,
            5,
            1,
            2,
            4,
            1,
            6,
            8,
            1,
            3,
            6,
            1,
            1,
            1,
            1,
            2,
            3,
            4,
            5
        };

        int[] b_vec =
        {
            1,
            5,
            1,
            5,
            10,
            20,
            5,
            6,
            16,
            5,
            10,
            12,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5
        };

        double[] fx_vec =
        {
            0.5000000000000000E+00,
            0.4999714850534485E+00,
            0.4996034370170990E+00,
            0.7496993658293228E+00,
            0.7504656462757382E+00,
            0.7514156325324275E+00,
            0.8999867031372156E+00,
            0.8997127554259699E+00,
            0.9002845660853669E+00,
            0.9500248817817622E+00,
            0.9500574946122442E+00,
            0.9501926400000000E+00,
            0.9750133887312993E+00,
            0.9900022327445249E+00,
            0.9949977837872073E+00,
            0.9989999621122122E+00,
            0.5687988496283079E+00,
            0.5351452100063650E+00,
            0.5143428032407864E+00,
            0.5000000000000000E+00
        };

        double[] x_vec =
        {
            1.00E+00,
            0.528E+00,
            1.89E+00,
            1.69E+00,
            1.60E+00,
            1.47E+00,
            4.06E+00,
            3.05E+00,
            2.09E+00,
            6.61E+00,
            3.71E+00,
            3.00E+00,
            10.01E+00,
            16.26E+00,
            22.78E+00,
            47.18E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
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
            a = 0;
            b = 0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            a = a_vec[n_data - 1];
            b = b_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static void f_noncentral_cdf_values(ref int n_data, ref int n1, ref int n2, ref double lambda,
            ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = NoncentralFRatioDistribution [ n1, n2, lambda ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 August 2004
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
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref int N1, ref int N2, the numerator and denominator
        //    degrees of freedom.
        //
        //    Output, ref double LAMBDA, the noncentrality parameter.
        //
        //    Output, ref double X, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 22;

        double[] fx_vec =
        {
            0.5000000000000000E+00,
            0.6367825323508774E+00,
            0.5840916116305482E+00,
            0.3234431872392788E+00,
            0.4501187879813550E+00,
            0.6078881441188312E+00,
            0.7059275551414605E+00,
            0.7721782003263727E+00,
            0.8191049017635072E+00,
            0.3170348430749965E+00,
            0.4327218008454471E+00,
            0.4502696915707327E+00,
            0.4261881186594096E+00,
            0.6753687206341544E+00,
            0.4229108778879005E+00,
            0.6927667261228938E+00,
            0.3632174676491226E+00,
            0.4210054012695865E+00,
            0.4266672258818927E+00,
            0.4464016600524644E+00,
            0.8445888579504827E+00,
            0.4339300273343604E+00
        };

        double[] lambda_vec =
        {
            0.00E+00,
            0.00E+00,
            0.25E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            2.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            2.00E+00,
            1.00E+00,
            1.00E+00,
            0.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00
        };

        int[] n1_vec =
        {
            1, 1, 1, 1,
            1, 1, 1, 1,
            1, 1, 2, 2,
            3, 3, 4, 4,
            5, 5, 6, 6,
            8, 16
        };

        int[] n2_vec =
        {
            1, 5, 5, 5,
            5, 5, 5, 5,
            5, 5, 5, 10,
            5, 5, 5, 5,
            1, 5, 6, 12,
            16, 8
        };

        double[] x_vec =
        {
            1.00E+00,
            1.00E+00,
            1.00E+00,
            0.50E+00,
            1.00E+00,
            2.00E+00,
            3.00E+00,
            4.00E+00,
            5.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            2.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            2.00E+00,
            2.00E+00
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
            n1 = 0;
            n2 = 0;
            lambda = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            n1 = n1_vec[n_data - 1];
            n2 = n2_vec[n_data - 1];
            lambda = lambda_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

}