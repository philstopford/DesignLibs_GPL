using System;
using Burkardt.CDFLib;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Truncated
    {
        public static double normal_truncated_ab_mean(double mu, double s, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_MEAN returns the mean of the truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviatione of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double NORMAL_TRUNCATED_AB_MEAN, the mean of the PDF.
        //
        {
            double alpha = (a - mu) / s;
            double beta = (b - mu) / s;

            double alpha_cdf = CDF.normal_01_cdf(alpha);
            double beta_cdf = CDF.normal_01_cdf(beta);

            double alpha_pdf = Normal.normal_01_pdf(alpha);
            double beta_pdf = Normal.normal_01_pdf(beta);

            double mean = mu + s * (alpha_pdf - beta_pdf) / (beta_cdf - alpha_cdf);

            return mean;
        }

        public static double normal_truncated_ab_pdf(double x, double mu, double s, double a,
            double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_PDF evaluates the truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double NORMAL_TRUNCATED_AB_PDF, the value of the PDF.
        //
        {
            double alpha = (a - mu) / s;
            double beta = (b - mu) / s;
            double xi = (x - mu) / s;

            double alpha_cdf = CDF.normal_01_cdf(alpha);
            double beta_cdf = CDF.normal_01_cdf(beta);
            double xi_pdf = Normal.normal_01_pdf(xi);

            double pdf = xi_pdf / (beta_cdf - alpha_cdf) / s;

            return pdf;
        }

        public static double normal_truncated_ab_sample(double mu, double s, double a, double b,
            ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_SAMPLE samples the truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Input/output, int &SEED, a seed for the random number
        //    generator.
        //
        //    Output, double NORMAL_TRUNCATED_AB_SAMPLE, a sample of the PDF.
        //
        {
            double alpha = (a - mu) / s;
            double beta = (b - mu) / s;

            double alpha_cdf = CDF.normal_01_cdf(alpha);
            double beta_cdf = CDF.normal_01_cdf(beta);

            double u = UniformRNG.r8_uniform_01(ref seed);
            double xi_cdf = alpha_cdf + u * (beta_cdf - alpha_cdf);
            double xi = CDF.normal_01_cdf_inv(xi_cdf);

            double x = mu + s * xi;

            return x;
        }

        public static double normal_truncated_ab_variance(double mu, double s, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_VARIANCE returns the variance of the truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double NORMAL_TRUNCATED_AB_VARIANCE, the variance of the PDF.
        //
        {
            double alpha = (a - mu) / s;
            double beta = (b - mu) / s;

            double alpha_pdf = Normal.normal_01_pdf(alpha);
            double beta_pdf = Normal.normal_01_pdf(beta);

            double alpha_cdf = CDF.normal_01_cdf(alpha);
            double beta_cdf = CDF.normal_01_cdf(beta);

            double variance = s * s * (1.0
                                       + (alpha * alpha_pdf - beta * beta_pdf) / (beta_cdf - alpha_cdf)
                                       - Math.Pow((alpha_pdf - beta_pdf) / (beta_cdf - alpha_cdf), 2));

            return variance;
        }

        public static double normal_truncated_a_mean(double mu, double s, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_MEAN returns the mean of the lower truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviatione of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_A_MEAN, the mean of the PDF.
        //
        {
            double alpha = (a - mu) / s;

            double alpha_cdf = CDF.normal_01_cdf(alpha);

            double alpha_pdf = Normal.normal_01_pdf(alpha);

            double mean = mu + s * alpha_pdf / (1.0 - alpha_cdf);

            return mean;
        }

        public static double normal_truncated_ab_moment(int order, double mu, double sigma,
                double a, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_AB_MOMENT: moments of the truncated Normal PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Phoebus Dhrymes,
            //    Moments of Truncated Normal Distributions,
            //    May 2005.
            //
            //  Parameters:
            //
            //    Input, int ORDER, the order of the moment.
            //    0 <= ORDER.
            //
            //    Input, double MU, SIGMA, the mean and standard deviation of the
            //    parent Normal distribution.
            //    0.0 < S.
            //
            //    Input, double A, B, the lower and upper truncation limits.
            //    A < B.
            //
            //    Output, double TRUNCATED_NORMAL_AB_MOMENT, the moment of the PDF.
            //
        {
            double a_cdf;
            double a_h;
            double a_pdf;
            double b_cdf;
            double b_h;
            double b_pdf;
            double ir;
            double irm1;
            double irm2;
            double moment;
            int r;

            if (order < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TRUNCATED_NORMAL_AB_MOMENT - Fatal error!");
                Console.WriteLine("  ORDER < 0.");
                return(1);
            }

            if (sigma <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("TRUNCATED_NORMAL_AB_MOMENT - Fatal error!");
                Console.WriteLine("  SIGMA <= 0.0.");
                return(1);
            }

            if (b <= a)
            {
                Console.WriteLine("");
                Console.WriteLine("TRUNCATED_NORMAL_AB_MOMENT - Fatal error!");
                Console.WriteLine("  B <= A.");
                return(1);
            }

            a_h = (a - mu) / sigma;
            a_pdf = Normal.normal_01_pdf(a_h);
            a_cdf = CDF.normal_01_cdf(a_h);

            if (a_cdf == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("TRUNCATED_NORMAL_AB_MOMENT - Fatal error!");
                Console.WriteLine("  PDF/CDF ratio fails, because A_CDF too small.");
                Console.WriteLine("  A_PDF = " + a_pdf + "");
                Console.WriteLine("  A_CDF = " + a_cdf + "");
                return (1);
            }

            b_h = (b - mu) / sigma;
            b_pdf = Normal.normal_01_pdf(b_h);
            b_cdf = CDF.normal_01_cdf(b_h);

            if (b_cdf == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("TRUNCATED_NORMAL_AB_MOMENT - Fatal error!");
                Console.WriteLine("  PDF/CDF ratio fails, because B_CDF too small.");
                Console.WriteLine("  B_PDF = " + b_pdf + "");
                Console.WriteLine("  B_CDF = " + b_cdf + "");
                return(1);
            }

            moment = 0.0;
            irm2 = 0.0;
            irm1 = 0.0;

            for (r = 0; r <= order; r++)
            {
                if (r == 0)
                {
                    ir = 1.0;
                }
                else if (r == 1)
                {
                    ir = -(b_pdf - a_pdf) / (b_cdf - a_cdf);
                }
                else
                {
                    ir = (double)(r - 1) * irm2
                         - (Math.Pow(b_h, r - 1) * b_pdf - Math.Pow(a_h, r - 1) * a_pdf)
                         / (b_cdf - a_cdf);
                }

                moment = moment + typeMethods.r8_choose(order, r) * Math.Pow(mu, order - r)
                                                      * Math.Pow(sigma, r) * ir;

                irm2 = irm1;
                irm1 = ir;
            }

            return moment;
        }
        
        public static double normal_truncated_a_moment ( int order, double mu, double sigma, 
                double a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_A_MOMENT: moments of the lower truncated Normal PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Phoebus Dhrymes,
            //    Moments of Truncated Normal Distributions,
            //    May 2005.
            //
            //  Parameters:
            //
            //    Input, int ORDER, the order of the moment.
            //    0 <= ORDER.
            //
            //    Input, double MU, SIGMA, the mean and standard deviation of the
            //    parent Normal distribution.
            //
            //    Input, double A, the lower truncation limit.
            //
            //    Output, double TRUNCATED_NORMAL_A_MOMENT, the moment of the PDF.
            //
        {
            double moment;

            moment = typeMethods.r8_mop ( order )
                     * truncated_normal_b_moment ( order, - mu, sigma, - a );

            return moment;
        }

        public static double normamomnormal_truncated_a_pdf(double x, double mu, double s, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_PDF evaluates the lower truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_A_PDF, the value of the PDF.
        //
        {
            double alpha = (a - mu) / s;
            double xi = (x - mu) / s;

            double alpha_cdf = CDF.normal_01_cdf(alpha);
            double xi_pdf = Normal.normal_01_pdf(xi);

            double pdf = xi_pdf / (1.0 - alpha_cdf) / s;

            return pdf;
        }

        public static double normal_truncated_a_pdf ( double x, double mu, double sigma, double a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_A_PDF evaluates the lower truncated Normal PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 January 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument of the PDF.
            //
            //    Input, double MU, SIGMA, the mean and standard deviation of the
            //    parent Normal distribution.
            //
            //    Input, double A, the lower truncation limit.
            //
            //    Output, double TRUNCATED_NORMAL_A_PDF, the value of the PDF.
            //
        {
            double alpha;
            double alpha_cdf;
            double pdf;
            double xi;
            double xi_pdf;

            if ( x < a )
            {
                pdf = 0.0;
            }
            else
            {
                alpha = ( a - mu ) / sigma;
                xi = ( x - mu ) / sigma;

                alpha_cdf = CDF.normal_01_cdf ( alpha );
                xi_pdf = Normal.normal_01_pdf ( xi );

                pdf = xi_pdf / ( 1.0 - alpha_cdf ) / sigma;
            }
  
            return pdf;
        }
        public static double normal_truncated_a_sample(double mu, double s, double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_SAMPLE samples the lower truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Input/output, int &SEED, a seed for the random number
        //    generator.
        //
        //    Output, double NORMAL_TRUNCATED_A_SAMPLE, a sample of the PDF.
        //
        {
            double alpha;
            double alpha_cdf;
            double u;
            double x;
            double xi;
            double xi_cdf;

            alpha = (a - mu) / s;

            alpha_cdf = CDF.normal_01_cdf(alpha);

            u = UniformRNG.r8_uniform_01(ref seed);
            xi_cdf = alpha_cdf + u * (1.0 - alpha_cdf);
            xi = CDF.normal_01_cdf_inv(xi_cdf);

            x = mu + s * xi;

            return x;
        }

        public static double normal_truncated_a_variance(double mu, double s, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_VARIANCE: variance of the lower truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_A_VARIANCE, the variance of the PDF.
        //
        {
            double alpha = (a - mu) / s;

            double alpha_pdf = Normal.normal_01_pdf(alpha);

            double alpha_cdf = CDF.normal_01_cdf(alpha);

            double variance = s * s * (1.0
                                       + alpha * alpha_pdf / (1.0 - alpha_cdf)
                                       - Math.Pow(alpha_pdf / (1.0 - alpha_cdf), 2));

            return variance;
        }

        public static double normal_truncated_b_mean(double mu, double s, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_MEAN returns the mean of the upper truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviatione of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_B_MEAN, the mean of the PDF.
        //
        {
            double beta = (b - mu) / s;

            double beta_cdf = CDF.normal_01_cdf(beta);

            double beta_pdf = Normal.normal_01_pdf(beta);

            double mean = mu - s * beta_pdf / beta_cdf;

            return mean;
        }

        public static double normal_truncated_b_moment ( int order, double mu, double sigma, double b )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_B_MOMENT: moments of the upper truncated Normal PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Phoebus Dhrymes,
            //    Moments of Truncated Normal Distributions,
            //    May 2005.
            //
            //  Parameters:
            //
            //    Input, int ORDER, the order of the moment.
            //    0 <= ORDER.
            //
            //    Input, double MU, SIGMA, the mean and standard deviation of the
            //    parent Normal distribution.
            //
            //    Input, double B, the upper truncation limit.
            //
            //    Output, double TRUNCATED_NORMAL_B_MOMENT, the moment of the PDF.
            //
        {
            double f;
            double h;
            double h_cdf;
            double h_pdf;
            double ir;
            double irm1;
            double irm2;
            double moment;
            int r;

            if ( order < 0 )
            {
                Console.WriteLine("");
                Console.WriteLine("TRUNCATED_NORMAL_B_MOMENT - Fatal error!");
                Console.WriteLine("  ORDER < 0.");
                return ( 1 );
            }

            h = ( b - mu ) / sigma;
            h_pdf = Normal.normal_01_pdf ( h );
            h_cdf = CDF.normal_01_cdf ( h );

            if ( h_cdf == 0.0 )
            {
                Console.WriteLine("");
                Console.WriteLine("TRUNCATED_NORMAL_B_MOMENT - Fatal error!");
                Console.WriteLine("  CDF((B-MU)/SIGMA) = 0.");
                return ( 1 );
            }

            f = h_pdf / h_cdf;

            moment = 0.0;
            irm2 = 0.0;
            irm1 = 0.0;

            for ( r = 0; r <= order; r++ )
            {
                if ( r == 0 )
                {
                    ir = 1.0;
                }
                else if ( r == 1 )
                {
                    ir = - f;
                }
                else
                {
                    ir = - Math.Pow ( h, r - 1 ) * f + ( double ) ( r - 1 ) * irm2;
                }

                moment = moment + typeMethods.r8_choose ( order, r ) * Math.Pow ( mu, order - r ) 
                                                         * Math.Pow ( sigma, r ) * ir;

                irm2 = irm1;
                irm1 = ir;
            }

            return moment;
        }
        
        public static double normal_truncated_b_pdf(double x, double mu, double s, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_PDF evaluates the upper truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_B_PDF, the value of the PDF.
        //
        {
            double beta = (b - mu) / s;
            double xi = (x - mu) / s;

            double beta_cdf = CDF.normal_01_cdf(beta);
            double xi_pdf = Normal.normal_01_pdf(xi);

            double pdf = xi_pdf / beta_cdf / s;

            return pdf;
        }

        public static double normal_truncated_b_sample(double mu, double s, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_SAMPLE samples the upper truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Input/output, int &SEED, a seed for the random number
        //    generator.
        //
        //    Output, double NORMAL_TRUNCATED_B_SAMPLE, a sample of the PDF.
        //
        {
            double beta = (b - mu) / s;

            double beta_cdf = CDF.normal_01_cdf(beta);

            double u = UniformRNG.r8_uniform_01(ref seed);
            double xi_cdf = u * beta_cdf;
            double xi = CDF.normal_01_cdf_inv(xi_cdf);

            double x = mu + s * xi;

            return x;
        }

        public static double normal_truncated_b_variance(double mu, double s, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_VARIANCE: variance of the upper truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_B_VARIANCE, the variance of the PDF.
        //
        {
            double beta = (b - mu) / s;

            double beta_pdf = Normal.normal_01_pdf(beta);

            double beta_cdf = CDF.normal_01_cdf(beta);

            double variance = s * s * (1.0
                                       - beta * beta_pdf / beta_cdf
                                       - Math.Pow(beta_pdf / beta_cdf, 2));

            return variance;
        }

        public static double truncated_normal_a_mean(double mu, double sigma, double a)
        {
            return normal_truncated_a_mean(mu, sigma, a);
        }

        public static double truncated_normal_a_moment(int order, double mu, double sigma,
            double a)
        {
            return normal_truncated_a_moment(order, mu, sigma, a);
        }
        
        public static double truncated_normal_ab_mean(double mu, double sigma, double a, double b)
        {
            return normal_truncated_ab_mean(mu, sigma, a, b);
        }

        public static double truncated_normal_ab_moment(int order, double mu, double sigma,
            double a, double b)
        {
            return normal_truncated_ab_moment(order, mu, sigma, a, b);
        }

        public static double truncated_normal_ab_pdf(double x, double mu, double sigma, double a,
            double b)
        {
            return normal_truncated_ab_pdf(x, mu, sigma, a, b);
        }

        public static double truncated_normal_ab_sample(double mu, double sigma, double a, double b,
            ref int seed)
        {
            return normal_truncated_ab_sample(mu, sigma, a, b, ref seed);
        }

        public static double truncated_normal_ab_variance(double mu, double sigma, double a,
            double b)
        {
            return normal_truncated_ab_variance(mu, sigma, a, b);
        }

        public static double truncated_normal_a_pdf(double x, double mu, double sigma, double a)
        {
            return normal_truncated_a_pdf(x, mu, sigma, a);
        }

        public static double truncated_normal_a_sample(double mu, double sigma, double a,
            ref int seed)
        {
            return normal_truncated_a_sample(mu, sigma, a, ref seed);
        }

        public static double truncated_normal_a_variance(double mu, double sigma, double a)
        {
            return normal_truncated_a_variance(mu, sigma, a);
        }

        public static double truncated_normal_b_mean(double mu, double sigma, double b)
        {
            return normal_truncated_b_mean(mu, sigma, b);
        }

        public static double truncated_normal_b_moment(int order, double mu, double sigma, double b)
        {
            return normal_truncated_b_moment(order, mu, sigma, b);
        }
        
        public static double truncated_normal_b_pdf(double x, double mu, double sigma, double b)
        {
            return normal_truncated_b_pdf(x, mu, sigma, b);
        }
        
        public static double truncated_normal_b_sample(double mu, double sigma, double b,
            ref int seed)
        {
            return normal_truncated_b_sample(mu, sigma, b, ref seed);
        }

        public static double truncated_normal_b_variance(double mu, double sigma, double b)
        {
            return normal_truncated_b_variance(mu, sigma, b);
        }

        
    }
}