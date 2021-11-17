using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Extreme
{
    public static double extreme_values_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXTREME_VALUES_CDF evaluates the Extreme Values CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXTREME_VALUES_CDF, the value of the CDF.
        //
    {
        double cdf;
        double y;

        y = (x - a) / b;

        cdf = Math.Exp(-Math.Exp(-y));

        return cdf;
    }

    public static double extreme_values_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXTREME_VALUES_CDF_INV inverts the Extreme Values CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXTREME_VALUES_CDF_INV, the corresponding argument of the CDF.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("EXTREME_VALUES_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            default:
            {
                double x = a - b * Math.Log(-Math.Log(cdf));

                return x;
            }
        }
    }

    public static void extreme_values_cdf_values(ref int n_data, ref double alpha, ref double beta,
            ref double x, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXTREME_VALUES_CDF_VALUES returns some values of the Extreme Values CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = ExtremeValuesDistribution [ alpha, beta ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 September 2004
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
        //    Output, double &ALPHA, the first parameter of the distribution.
        //
        //    Output, double &BETA, the second parameter of the distribution.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 12;

        double[] alpha_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.5000000000000000E+01
            }
            ;

        double[] beta_vec =
            {
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.5000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01
            }
            ;

        double[] fx_vec =
            {
                0.3678794411714423E+00,
                0.8734230184931166E+00,
                0.9818510730616665E+00,
                0.9975243173927525E+00,
                0.5452392118926051E+00,
                0.4884435800065159E+00,
                0.4589560693076638E+00,
                0.4409910259429826E+00,
                0.5452392118926051E+00,
                0.3678794411714423E+00,
                0.1922956455479649E+00,
                0.6598803584531254E-01
            }
            ;

        double[] x_vec =
            {
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01
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
            alpha = 0.0;
            beta = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            alpha = alpha_vec[n_data - 1];
            beta = beta_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static bool extreme_values_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXTREME_VALUES_CHECK checks the parameters of the Extreme Values CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, bool EXTREME_VALUES_CHECK, is true if the parameters are legal.
        //
    {
        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("EXTREME_VALUES_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double extreme_values_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXTREME_VALUES_MEAN returns the mean of the Extreme Values PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXTREME_VALUES_MEAN, the mean of the PDF.
        //
    {
        double mean;

        mean = a + b * Misc.euler_constant();

        return mean;
    }

    public static double extreme_values_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXTREME_VALUES_PDF evaluates the Extreme Values PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) =
        //      ( 1 / B ) * exp ( ( A - X ) / B - exp ( ( A - X ) / B  ) ).
        //
        //    The Extreme Values PDF is also known as the Fisher-Tippet PDF
        //    and the Log-Weibull PDF.
        //
        //    The special case A = 0 and B = 1 is the Gumbel PDF.
        //
        //    The Extreme Values PDF is the limiting distribution for the
        //    smallest or largest value in a large sample drawn from
        //    any of a great variety of distributions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Eric Weisstein, editor,
        //    CRC Concise Encylopedia of Mathematics,
        //    CRC Press, 1998.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXTREME_VALUES_PDF, the value of the PDF.
        //
    {
        double pdf = 1.0 / b * Math.Exp((a - x) / b - Math.Exp((a - x) / b));

        return pdf;
    }

    public static double extreme_values_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXTREME_VALUES_SAMPLE samples the Extreme Values PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double EXTREME_VALUES_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = extreme_values_cdf_inv(cdf, a, b);

        return x;
    }

    public static double extreme_values_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXTREME_VALUES_VARIANCE returns the variance of the Extreme Values PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXTREME_VALUES_VARIANCE, the variance of the PDF.
        //
    {
            

        double variance = Math.PI * Math.PI * b * b / 6.0;

        return variance;
    }
}