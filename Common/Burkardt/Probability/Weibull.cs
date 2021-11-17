using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Weibull
{
    public static double weibull_cdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_CDF evaluates the Weibull CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //    A <= X.
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double CDF, the value of the CDF.
        //
    {
        double cdf;
        double y;

        if (x < a)
        {
            cdf = 0.0;
        }
        else
        {
            y = (x - a) / b;
            cdf = 1.0 - 1.0 / Math.Exp(Math.Pow(y, c));
        }

        return cdf;
    }

    public static double weibull_cdf_inv(double cdf, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_CDF_INV inverts the Weibull CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 < CDF < 1.0.
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double WEIBULL_CDF_INV, the corresponding argument of the CDF.
        //
    {
        double x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("WEIBULL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            default:
                x = a + b * Math.Pow(-Math.Log(1.0 - cdf), 1.0 / c);

                return x;
        }
    }

    public static void weibull_cdf_values(ref int n_data, ref double alpha, ref double beta,
            ref double x, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_CDF_VALUES returns some values of the Weibull CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = WeibullDistribution [ alpha, beta ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2004
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
        };

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
        };

        double[] fx_vec =
        {
            0.8646647167633873E+00,
            0.9816843611112658E+00,
            0.9975212478233336E+00,
            0.9996645373720975E+00,
            0.6321205588285577E+00,
            0.4865828809674080E+00,
            0.3934693402873666E+00,
            0.3296799539643607E+00,
            0.8946007754381357E+00,
            0.9657818816883340E+00,
            0.9936702845725143E+00,
            0.9994964109502630E+00
        };

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

    public static bool weibull_check(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_CHECK checks the parameters of the Weibull CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, bool WEIBULL_CHECK, is true if the parameters are legal.
        //
    {

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("WEIBULL_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
        }

        switch (c)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("WEIBULL_CHECK - Warning!");
                Console.WriteLine("  C <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double weibull_mean(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_MEAN returns the mean of the Weibull PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
    {
        double mean = b * Helpers.Gamma((c + 1.0) / c) + a;

        return mean;
    }

    public static double weibull_pdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_PDF evaluates the Weibull PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B,C;X) = ( C / B ) * ( ( X - A ) / B )^( C - 1 )
        //     * EXP ( - ( ( X - A ) / B )^C ).
        //
        //    The Weibull PDF is also known as the Frechet PDF.
        //
        //    WEIBULL_PDF(A,B,1;X) is the Exponential PDF.
        //
        //    WEIBULL_PDF(0,1,2;X) is the Rayleigh PDF.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    A <= X
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double PDF, the value of the PDF.
        //
    {
        double pdf;
        double y;

        if (x < a)
        {
            pdf = 0.0;
        }
        else
        {
            y = (x - a) / b;

            pdf = c / b * Math.Pow(y, c - 1.0) / Math.Exp(Math.Pow(y, c));
        }

        return pdf;
    }

    public static double weibull_sample(double a, double b, double c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_SAMPLE samples the Weibull PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double WEIBULL_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = weibull_cdf_inv(cdf, a, b, c);

        return x;
    }

    public static double weibull_variance(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_VARIANCE returns the variance of the Weibull PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
    {
        double g1 = Helpers.Gamma((c + 2.0) / c);
        double g2 = Helpers.Gamma((c + 1.0) / c);

        double variance = b * b * (g1 - g2 * g2);

        return variance;
    }

    public static double weibull_discrete_cdf(int x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_DISCRETE_CDF evaluates the Discrete Weibull CDF.
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
        //    Input, int X, the argument of the CDF.
        //    0 <= X.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= A <= 1.0,
        //    0.0 < B.
        //
        //    Output, double WEIBULL_DISCRETE_CDF, the value of the CDF.
        //
    {
        double cdf = x switch
        {
            < 0 => 0.0,
            _ => 1.0 - Math.Pow(1.0 - a, Math.Pow(x + 1, b))
        };

        return cdf;
    }

    public static int weibull_discrete_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_DISCRETE_CDF_INV inverts the Discrete Weibull CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= A <= 1.0,
        //    0.0 < B.
        //
        //    Output, int WEIBULL_DISCRETE_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("WEIBULL_DISCRETE_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            default:
            {
                int x = (int)Math.Ceiling(
                    Math.Pow(Math.Log(1.0 - cdf) / Math.Log(1.0 - a), 1.0 / b) - 1.0);

                return x;
            }
        }
    }

    public static bool weibull_discrete_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_DISCRETE_CHECK checks the parameters of the discrete Weibull CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= A <= 1.0,
        //    0.0 < B.
        //
        //    Output, bool WEIBULL_DISCRETE_CHECK, is true if the parameters
        //    are legal.
        //
    {
        switch (a)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("WEIBULL_DISCRETE_CHECK - Warning!");
                Console.WriteLine("  A < 0 or 1 < A.");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("WEIBULL_DISCRETE_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double weibull_discrete_pdf(int x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_DISCRETE_PDF evaluates the discrete Weibull PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = ( 1 - A )^X^B - ( 1 - A )^(X+1)^B.
        //
        //    WEIBULL_DISCRETE_PDF(A,1;X) = GEOMETRIC_PDF(A;X)
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
        //    Input, int X, the argument of the PDF.
        //    0 <= X
        //
        //    Input, double A, B, the parameters that define the PDF.
        //    0 <= A <= 1,
        //    0 < B.
        //
        //    Output, double WEIBULL_DISCRETE_PDF, the value of the PDF.
        //
    {
        double pdf = x switch
        {
            < 0 => 0.0,
            _ => Math.Pow(1.0 - a, Math.Pow(x, b)) - Math.Pow(1.0 - a, Math.Pow(x + 1, b))
        };

        return pdf;
    }

    public static int weibull_discrete_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_DISCRETE_SAMPLE samples the discrete Weibull PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= A <= 1.0,
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int WEIBULL_DISCRETE_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        int x = weibull_discrete_cdf_inv(cdf, a, b);

        return x;
    }
}