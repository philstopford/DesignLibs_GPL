using System;
using Burkardt.CDFLib;
using Burkardt.Uniform;

namespace Burkardt.PDFLib;

public static partial class PDF
{
    public static bool log_normal_check(double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_CHECK checks the parameters of the Log Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, bool LOG_NORMAL_CHECK, is true if the parameters are legal.
        //
    {
        bool check = true;

        switch (sigma)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_CHECK - Fatal error!");
                Console.WriteLine("  SIGMA <= 0.");
                check = false;
                break;
        }

        return check;
    }

    public static double log_normal_mean(double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_MEAN returns the mean of the Log Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, double LOG_NORMAL_MEAN, the mean of the PDF.
        //
    {
        double mean = Math.Exp(mu + 0.5 * sigma * sigma);

        return mean;
    }

    public static double log_normal_pdf(double x, double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_PDF evaluates the Log Normal PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X)
        //      = Math.Exp ( - 0.5 * ( ( log ( X ) - MU ) / SIGMA )^2 )
        //        / ( SIGMA * X * sqrt ( 2 * PI ) )
        //
        //    The Log Normal PDF is also known as the Cobb-Douglas PDF,
        //    and as the Antilog_normal PDF.
        //
        //    The Log Normal PDF describes a variable X whose logarithm
        //    is normally distributed.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    0.0 < X
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, double LOG_NORMAL_PDF, the value of the PDF.
        //
    {
        double pdf;

        switch (x)
        {
            case <= 0.0:
                pdf = 0.0;
                break;
            default:
                double y = (Math.Log(x) - mu) / sigma;
                pdf = Math.Exp(-0.5 * y * y) / (sigma * x * Math.Sqrt(2.0 * Math.PI));
                break;
        }

        return pdf;
    }

    public static double log_normal_sample(double mu, double sigma, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_SAMPLE samples the Log Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double LOG_NORMAL_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = CDF.log_normal_cdf_inv(cdf, mu, sigma);

        return x;
    }

    public static double log_normal_variance(double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_VARIANCE returns the variance of the Log Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, double LOG_NORMAL_VARIANCE, the variance of the PDF.
        //
    {
        double variance = Math.Exp(2.0 * mu + sigma * sigma) * (Math.Exp(sigma * sigma) - 1.0);

        return variance;
    }
    //
}