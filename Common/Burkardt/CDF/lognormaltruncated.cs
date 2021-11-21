using System;
using Burkardt.PDFLib;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double log_normal_truncated_ab_cdf(double x, double mu, double sigma,
            double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_TRUNCATED_AB_CDF evaluates the Log Normal truncated AB CDF.
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
        //    0.0 < X.
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //    A < B.
        //
        //    Output, double LOG_NORMAL_TRUNCATED_AB_CDF, the value of the CDF.
        //
    {
        double cdf;

        bool check = PDF.log_normal_truncated_ab_check(mu, sigma, a, b);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_CDF - Fatal error!");
                Console.WriteLine("  Parameters are not legal.");
                return 1;
        }

        if (x <= a)
        {
            cdf = 0.0;
        }
        else if (b <= x)
        {
            cdf = 1.0;
        }
        else
        {
            double lncdf_a = log_normal_cdf(a, mu, sigma);
            double lncdf_b = log_normal_cdf(b, mu, sigma);
            double lncdf_x = log_normal_cdf(x, mu, sigma);

            cdf = (lncdf_x - lncdf_a) / (lncdf_b - lncdf_a);
        }

        return cdf;
    }

    public static double log_normal_truncated_ab_cdf_inv(double cdf, double mu, double sigma,
            double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_TRUNCATED_AB_CDF_INV inverts the Log Normal truncated AB CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //    A < B.
        //
        //    Input, double LOG_NORMAL_TRUNCATED_AB_CDF_INV, the corresponding argument.
        //
    {
        bool check = PDF.log_normal_truncated_ab_check(mu, sigma, a, b);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_CDF_INV - Fatal error!");
                Console.WriteLine("  Parameters are not legal.");
                return 1;
            default:
                double x;
                switch (cdf)
                {
                    case <= 0.0:
                        x = a;
                        break;
                    case >= 1.0:
                        x = b;
                        break;
                    default:
                        double lncdf_a = log_normal_cdf(a, mu, sigma);
                        double lncdf_b = log_normal_cdf(b, mu, sigma);

                        double lncdf_x = lncdf_a + cdf * (lncdf_b - lncdf_a);
                        x = log_normal_cdf_inv(lncdf_x, mu, sigma);
                        break;
                }

                return x;
        }
    }
}