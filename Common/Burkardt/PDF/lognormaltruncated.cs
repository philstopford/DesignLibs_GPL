using System;
using Burkardt.CDFLib;
using Burkardt.Uniform;

namespace Burkardt.PDFLib;

public static partial class PDF
{
    public static bool log_normal_truncated_ab_check(double mu, double sigma, double a,
            double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_TRUNCATED_AB_CHECK checks the Log Normal truncated AB PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 March 2016
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
        //    Input, double A, B, the lower and upper truncation limits.
        //    A < B.
        //
        //    Output, bool LOG_NORMAL_TRUNCATED_AB_CHECK, is true if the parameters 
        //    are legal.
        //
    {
        bool check = true;

        switch (sigma)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_CHECK - Fatal error!");
                Console.WriteLine("  SIGMA <= 0.");
                check = false;
                break;
        }

        return check;
    }

    public static double log_normal_truncated_ab_mean(double mu, double sigma, double a,
            double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_TRUNCATED_AB_MEAN: mean of the Log Normal truncated AB PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 March 2016
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
        //    Input, double A, B, the lower and upper truncation limits.
        //    A < B.
        //
        //    Output, double LOG_NORMAL_TRUNCATED_AB_MEAN, the mean of the PDF.
        //
    {
        bool check = log_normal_truncated_ab_check(mu, sigma, a, b);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_MEAN - Fatal error!");
                Console.WriteLine("  Parameters are not legal.");
                return 1;
        }

        double a0 = (Math.Log(a) - mu) / sigma;
        double b0 = (Math.Log(b) - mu) / sigma;

        double c1 = CDF.normal_01_cdf(sigma - a0);
        double c2 = CDF.normal_01_cdf(sigma - b0);
        double c3 = CDF.normal_01_cdf(+a0);
        double c4 = CDF.normal_01_cdf(+b0);

        double ln_mean = Math.Exp(mu + 0.5 * sigma * sigma);

        double mean = ln_mean * (c1 - c2) / (c4 - c3);

        return mean;
    }

    public static double log_normal_truncated_ab_pdf(double x, double mu, double sigma,
            double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_TRUNCATED_AB_PDF evaluates the Log Normal truncated AB PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 March 2016
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
        //    Input, double A, B, the lower and upper truncation limits.
        //    A < B.
        //
        //    Output, double LOG_NORMAL_TRUNCATED_AB_PDF, the value of the PDF.
        //
    {
        double pdf;

        bool check = log_normal_truncated_ab_check(mu, sigma, a, b);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_PDF - Fatal error!");
                Console.WriteLine("  Parameters are not legal.");
                return 1;
        }

        if (x <= a)
        {
            pdf = 0.0;
        }
        else if (b <= x)
        {
            pdf = 0.0;
        }
        else
        {
            double lncdf_a = CDF.log_normal_cdf(a, mu, sigma);
            double lncdf_b = CDF.log_normal_cdf(b, mu, sigma);
            double lnpdf_x = log_normal_pdf(x, mu, sigma);

            pdf = lnpdf_x / (lncdf_b - lncdf_a);
        }

        return pdf;
    }

    public static double log_normal_truncated_ab_sample(double mu, double sigma, double a,
            double b, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_TRUNCATED_AB_SAMPLE samples the Log Normal truncated AB PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 March 2016
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
        //    Input, double A, B, the lower and upper truncation limits.
        //    A < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double LOG_NORMAL_TRUNCATED_AB_SAMPLE, a sample of the PDF.
        //
    {
        bool check = log_normal_truncated_ab_check(mu, sigma, a, b);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_SAMPLE - Fatal error!");
                Console.WriteLine("  Parameters are not legal.");
                return 1;
        }

        double lncdf_a = CDF.log_normal_cdf(a, mu, sigma);
        double lncdf_b = CDF.log_normal_cdf(b, mu, sigma);

        double cdf = UniformRNG.r8_uniform_ab(lncdf_a, lncdf_b, ref seed);

        double x = CDF.log_normal_cdf_inv(cdf, mu, sigma);

        return x;
    }

    public static double log_normal_truncated_ab_variance(double mu, double sigma, double a,
            double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_TRUNCATED_AB_VARIANCE: variance of Log Normal truncated AB PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 March 2016
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
        //    Input, double A, B, the lower and upper truncation limits.
        //    A < B.
        //
        //    Output, double LOG_NORMAL_TRUNCATED_AB_VARIANCE, the variance of the PDF.
        //
    {
        bool check = log_normal_truncated_ab_check(mu, sigma, a, b);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_VARIANCE - Fatal error!");
                Console.WriteLine("  Parameters are not legal.");
                return 1;
        }

        double mean = log_normal_truncated_ab_mean(mu, sigma, a, b);

        double a0 = (Math.Log(a) - mu) / sigma;
        double b0 = (Math.Log(b) - mu) / sigma;

        double c1 = CDF.normal_01_cdf(2.0 * sigma - a0);
        double c2 = CDF.normal_01_cdf(2.0 * sigma - b0);
        double c3 = CDF.normal_01_cdf(+a0);
        double c4 = CDF.normal_01_cdf(+b0);

        double ln_xsquared = Math.Exp(2.0 * mu + 2.0 * sigma * sigma);

        double lntab_xsquared = ln_xsquared * (c1 - c2) / (c4 - c3);

        double variance = lntab_xsquared - mean * mean;

        return variance;
    }
}