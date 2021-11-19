using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class LogSeries
{
    public static double log_series_cdf(double x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_SERIES_CDF evaluates the Logarithmic Series CDF.
        //
        //  Discussion:
        //
        //    Simple summation is used, with a recursion to generate successive
        //    values of the PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Thanks:
        //
        //    Oscar van Vlijmen
        //
        //  Parameters:
        //
        //    Input, int X, the argument of the PDF.
        //    0 < X
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A < 1.0.
        //
        //    Output, double LOG_SERIES_CDF, the value of the CDF.
        //
    {
        int x2;

        double cdf = 0.0;
        double pdf = 0;

        for (x2 = 1; x2 <= x; x2++)
        {
            pdf = x2 switch
            {
                1 => -a / Math.Log(1.0 - a),
                _ => (x2 - 1) * a * pdf / x2
            };

            cdf += pdf;
        }

        return cdf;
    }

    public static int log_series_cdf_inv(double cdf, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_SERIES_CDF_INV inverts the Logarithmic Series CDF.
        //
        //  Discussion:
        //
        //    Simple summation is used.  The only protection against an
        //    infinite loop caused by roundoff is that X cannot be larger
        //    than 1000.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A < 1.0.
        //
        //    Output, double LOG_SERIES_CDF_INV, the argument X of the CDF for which
        //    CDF(X-1) <= CDF <= CDF(X).
        //
    {
        const int xmax = 1000;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("LOG_SERIES_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double cdf2 = 0.0;
        double pdf = 0;
        int x = 1;

        while (cdf2 < cdf && x < xmax)
        {
            pdf = x switch
            {
                1 => -a / Math.Log(1.0 - a),
                _ => (x - 1) * a * pdf / x
            };

            cdf2 += pdf;

            x += 1;
        }

        return x;
    }

    public static void log_series_cdf_values(ref int n_data, ref double t, ref int n, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_SERIES_CDF_VALUES returns some values of the log series CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`DiscreteDistributions`]
        //      dist = LogSeriesDistribution [ t ]
        //      CDF [ dist, n ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 August 2004
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
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &T, the parameter of the function.
        //
        //    Output, int &N, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 29;

        double[] fx_vec =
        {
            0.9491221581029903E+00,
            0.9433541128559735E+00,
            0.9361094611773272E+00,
            0.9267370278044118E+00,
            0.9141358246245129E+00,
            0.8962840235449100E+00,
            0.8690148741955517E+00,
            0.8221011541254772E+00,
            0.7213475204444817E+00,
            0.6068261510845583E+00,
            0.5410106403333613E+00,
            0.4970679476476894E+00,
            0.4650921887927060E+00,
            0.4404842934597863E+00,
            0.4207860535926143E+00,
            0.4045507673897055E+00,
            0.3908650337129266E+00,
            0.2149757685421097E+00,
            0.0000000000000000E+00,
            0.2149757685421097E+00,
            0.3213887739704539E+00,
            0.3916213575531612E+00,
            0.4437690508633213E+00,
            0.4850700239649681E+00,
            0.5191433267738267E+00,
            0.5480569580144867E+00,
            0.5731033910767085E+00,
            0.5951442521714636E+00,
            0.6147826594068904E+00
        };

        int[] n_vec =
        {
            1, 1, 1, 1, 1,
            1, 1, 1, 1, 1,
            1, 1, 1, 1, 1,
            1, 1, 1, 0, 1,
            2, 3, 4, 5, 6,
            7, 8, 9, 10
        };

        double[] t_vec =
        {
            0.1000000000000000E+00,
            0.1111111111111111E+00,
            0.1250000000000000E+00,
            0.1428571428571429E+00,
            0.1666666666666667E+00,
            0.2000000000000000E+00,
            0.2500000000000000E+00,
            0.3333333333333333E+00,
            0.5000000000000000E+00,
            0.6666666666666667E+00,
            0.7500000000000000E+00,
            0.8000000000000000E+00,
            0.8333333333333333E+00,
            0.8571485714857149E+00,
            0.8750000000000000E+00,
            0.8888888888888889E+00,
            0.9000000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00
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
            t = 0.0;
            n = 0;
            fx = 0.0;
        }
        else
        {
            t = t_vec[n_data - 1];
            n = n_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static bool log_series_check(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_SERIES_CHECK checks the parameter of the Logarithmic Series PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A < 1.0.
        //
        //    Output, bool LOG_SERIES_CHECK, is true if the parameters are legal.
        //
    {
        switch (a)
        {
            case <= 0.0:
            case >= 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("LOG_SERIES_CHECK - Warning!");
                Console.WriteLine("  A <= 0.0 or 1.0 <= A");
                return false;
            default:
                return true;
        }
    }

    public static double log_series_mean(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_SERIES_MEAN returns the mean of the Logarithmic Series PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A < 1.0.
        //
        //    Output, double LOG_SERIES_MEAN, the mean of the PDF.
        //
    {
        double mean = -a / ((1.0 - a) * Math.Log(1.0 - a));

        return mean;
    }

    public static double log_series_pdf(int x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_SERIES_PDF evaluates the Logarithmic Series PDF.
        //
        //  Discussion:
        //
        //    PDF(A;X) = - A^X / ( X * log ( 1 - A ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the argument of the PDF.
        //    0 < X
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A < 1.0.
        //
        //    Output, double LOG_SERIES_PDF, the value of the PDF.
        //
    {
        double pdf = x switch
        {
            <= 0 => 0.0,
            _ => -Math.Pow(a, x) / (x * Math.Log(1.0 - a))
        };

        return pdf;
    }

    public static int log_series_sample(double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_SERIES_SAMPLE samples the Logarithmic Series PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Luc Devroye,
        //    Non-Uniform Random Variate Generation,
        //    Springer-Verlag, New York, 1986, page 547.
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A < 1.0.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int LOG_SERIES_SAMPLE, a sample of the PDF.
        //
    {
        double u = UniformRNG.r8_uniform_01(ref seed);
        double v = UniformRNG.r8_uniform_01(ref seed);

        int x = (int) (1.0 + Math.Log(v) / Math.Log(1.0 - Math.Pow(1.0 - a, u)));

        return x;
    }

    public static double log_series_variance(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_SERIES_VARIANCE returns the variance of the Logarithmic Series PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A < 1.0.
        //
        //    Output, double LOG_SERIES_VARIANCE, the variance of the PDF.
        //
    {
        double alpha = -1.0 / Math.Log(1.0 - a);

        double variance = a * alpha * (1.0 - alpha * a) / Math.Pow(1.0 - a, 2);

        return variance;
    }
}