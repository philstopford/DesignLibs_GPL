using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Laplace
{
    public static void laplace_cdf_values(ref int n_data, ref double mu, ref double beta, ref double x,
            ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAPLACE_CDF_VALUES returns some values of the Laplace CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = LaplaceDistribution [ mu, beta ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 August 2004
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
        //    Output, double &MU, the mean of the distribution.
        //
        //    Output, double &BETA, the shape parameter.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 12;

        double[] beta_vec =
        {
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
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
            0.5000000000000000E+00,
            0.8160602794142788E+00,
            0.9323323583816937E+00,
            0.9751064658160680E+00,
            0.6967346701436833E+00,
            0.6417343447131054E+00,
            0.6105996084642976E+00,
            0.5906346234610091E+00,
            0.5000000000000000E+00,
            0.3032653298563167E+00,
            0.1839397205857212E+00,
            0.1115650800742149E+00
        };

        double[] mu_vec =
        {
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.4000000000000000E+01
        };

        double[] x_vec =
        {
            0.0000000000000000E+01,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01
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
            mu = 0.0;
            beta = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            mu = mu_vec[n_data - 1];
            beta = beta_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static double laplace_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAPLACE_CDF evaluates the Laplace CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 February 1999
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
        //    Output, double LAPLACE_CDF, the value of the PDF.
        //
    {
        double cdf;

        double y = (x - a) / b;

        if (x <= a)
        {
            cdf = 0.5 * Math.Exp(y);
        }
        else
        {
            cdf = 1.0 - 0.5 * Math.Exp(-y);
        }

        return cdf;
    }

    public static double laplace_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAPLACE_CDF_INV inverts the Laplace CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 1999
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
        //    Output, double LAPLACE_CDF_INV, the corresponding argument.
        //
    {
        double x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("LAPLACE_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case <= 0.5:
                x = a + b * Math.Log(2.0 * cdf);
                break;
            default:
                x = a - b * Math.Log(2.0 * (1.0 - cdf));
                break;
        }

        return x;
    }

    public static bool laplace_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAPLACE_CHECK checks the parameters of the Laplace PDF.
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
        //    Output, bool LAPLACE_CHECK, is true if the parameters are legal.
        //
    {
        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("LAPLACE_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double laplace_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAPLACE_MEAN returns the mean of the Laplace PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 1999
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
        //    Output, double LAPLACE_MEAN, the mean of the PDF.
        //
    {
        return a;
    }

    public static double laplace_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAPLACE_PDF evaluates the Laplace PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = exp ( - abs ( X - A ) / B ) / ( 2 * B )
        //
        //  Discussion:
        //
        //    The Laplace PDF is also known as the Double Exponential PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 February 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double LAPLACE_PDF, the value of the PDF.
        //
    {
        double pdf = Math.Exp(-Math.Abs(x - a) / b) / (2.0 * b);

        return pdf;
    }

    public static double laplace_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAPLACE_SAMPLE samples the Laplace PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 February 1999
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
        //    Output, double LAPLACE_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = laplace_cdf_inv(cdf, a, b);

        return x;
    }

    public static double laplace_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAPLACE_VARIANCE returns the variance of the Laplace PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 1999
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
        //    Output, double LAPLACE_VARIANCE, the variance of the PDF.
        //
    {
        double variance = 2.0 * b * b;

        return variance;
    }
}