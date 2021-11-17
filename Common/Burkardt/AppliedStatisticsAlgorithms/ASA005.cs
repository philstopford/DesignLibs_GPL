using System;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static double student_noncentral_cdf_values(ref int n_data, ref int df, ref double lambda,
            ref double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = NoncentralStudentTDistribution [ df, lambda ]
        //      CDF [ dist, x ]
        //
        //    Mathematica seems to have some difficulty computing this function
        //    to the desired number of digits.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 September 2004
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
        //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int *DF, double *LAMBDA, the parameters of the
        //    function.
        //
        //    Output, double *X, the argument of the function.
        //
        //    Output, double *FX, the value of the function.
        //
    {
        const int N_MAX = 30;

        int[] df_vec =
        {
            1, 2, 3,
            1, 2, 3,
            1, 2, 3,
            1, 2, 3,
            1, 2, 3,
            15, 20, 25,
            1, 2, 3,
            10, 10, 10,
            10, 10, 10,
            10, 10, 10
        };

        double[] fx_vec =
        {
            0.8975836176504333E+00,
            0.9522670169E+00,
            0.9711655571887813E+00,
            0.8231218864E+00,
            0.9049021510E+00,
            0.9363471834E+00,
            0.7301025986E+00,
            0.8335594263E+00,
            0.8774010255E+00,
            0.5248571617E+00,
            0.6293856597E+00,
            0.6800271741E+00,
            0.20590131975E+00,
            0.2112148916E+00,
            0.2074730718E+00,
            0.9981130072E+00,
            0.9994873850E+00,
            0.9998391562E+00,
            0.168610566972E+00,
            0.16967950985E+00,
            0.1701041003E+00,
            0.9247683363E+00,
            0.7483139269E+00,
            0.4659802096E+00,
            0.9761872541E+00,
            0.8979689357E+00,
            0.7181904627E+00,
            0.9923658945E+00,
            0.9610341649E+00,
            0.8688007350E+00
        };

        double[] lambda_vec =
        {
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            4.0E+00,
            4.0E+00,
            4.0E+00,
            7.0E+00,
            7.0E+00,
            7.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            2.0E+00,
            3.0E+00,
            4.0E+00,
            2.0E+00,
            3.0E+00,
            4.0E+00,
            2.0E+00,
            3.0E+00,
            4.0E+00
        };

        double[] x_vec =
        {
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            3.00E+00,
            15.00E+00,
            15.00E+00,
            15.00E+00,
            0.05E+00,
            0.05E+00,
            0.05E+00,
            4.00E+00,
            4.00E+00,
            4.00E+00,
            5.00E+00,
            5.00E+00,
            5.00E+00,
            6.00E+00,
            6.00E+00,
            6.00E+00
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        double fx;
        if (N_MAX < n_data)
        {
            n_data = 0;
            df = 0;
            lambda = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            df = df_vec[n_data - 1];
            lambda = lambda_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }

        return fx;
    }

}