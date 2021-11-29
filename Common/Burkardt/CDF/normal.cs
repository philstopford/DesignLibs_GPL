using System;
using Burkardt.Types;

namespace Burkardt.CDFLib;

public static partial class CDF
{

    public static double normal_01_cdf(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_CDF evaluates the Normal 01 CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 February 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    A G Adams,
        //    Areas Under the Normal Curve,
        //    Algorithm 39,
        //    Computer j.,
        //    Volume 12, pages 197-198, 1969.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Output, double CDF, the value of the CDF.
        //
    {
        const double a1 = 0.398942280444;
        const double a2 = 0.399903438504;
        const double a3 = 5.75885480458;
        const double a4 = 29.8213557808;
        const double a5 = 2.62433121679;
        const double a6 = 48.6959930692;
        const double a7 = 5.92885724438;
        const double b0 = 0.398942280385;
        const double b1 = 3.8052E-08;
        const double b2 = 1.00000615302;
        const double b3 = 3.98064794E-04;
        const double b4 = 1.98615381364;
        const double b5 = 0.151679116635;
        const double b6 = 5.29330324926;
        const double b7 = 4.8385912808;
        const double b8 = 15.1508972451;
        const double b9 = 0.742380924027;
        const double b10 = 30.789933034;
        const double b11 = 3.99019417011;
        double q;
        double y;
        switch (Math.Abs(x))
        {
            //
            //  |X| <= 1.28.
            //
            case <= 1.28:
                y = 0.5 * x * x;

                q = 0.5 - Math.Abs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5
                                                                         + a6 / (y + a7))));
                //
                //  1.28 < |X| <= 12.7
                //
                break;
            case <= 12.7:
                y = 0.5 * x * x;

                q = Math.Exp(-y) * b0 / (Math.Abs(x) - b1
                                         + b2 / (Math.Abs(x) + b3
                                                             + b4 / (Math.Abs(x) - b5
                                                                     + b6 / (Math.Abs(x) + b7
                                                                             - b8 / (Math.Abs(x) + b9
                                                                                 + b10 / (Math.Abs(x) + b11))))));
                //
                //  12.7 < |X|
                //
                break;
            default:
                q = 0.0;
                break;
        }

        double cdf = x switch
        {
            //
            //  Take account of negative X.
            //
            < 0.0 => q,
            _ => 1.0 - q
        };

        return cdf;
    }

    public static double normal_01_cdf_inv(double p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_CDF_INV inverts the standard normal CDF.
        //
        //  Discussion:
        //
        //    The result is accurate to about 1 part in 10**16.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 December 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Michael Wichura.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Michael Wichura,
        //    The Percentage Points of the Normal Distribution,
        //    Algorithm AS 241,
        //    Applied Statistics,
        //    Volume 37, Number 3, pages 477-484, 1988.
        //
        //  Parameters:
        //
        //    Input, double P, the value of the cumulative probability
        //    densitity function.  0 < P < 1.  If P is outside this range, an
        //    "infinite" value is returned.
        //
        //    Output, double NORMAL_01_CDF_INV, the normal deviate value
        //    with the property that the probability of a standard normal deviate being
        //    less than or equal to this value is P.
        //
    {
        double[] a =
        {
            3.3871328727963666080, 1.3314166789178437745e+2,
            1.9715909503065514427e+3, 1.3731693765509461125e+4,
            4.5921953931549871457e+4, 6.7265770927008700853e+4,
            3.3430575583588128105e+4, 2.5090809287301226727e+3
        };
        double[] b =
        {
            1.0, 4.2313330701600911252e+1,
            6.8718700749205790830e+2, 5.3941960214247511077e+3,
            2.1213794301586595867e+4, 3.9307895800092710610e+4,
            2.8729085735721942674e+4, 5.2264952788528545610e+3
        };
        double[] c =
        {
            1.42343711074968357734, 4.63033784615654529590,
            5.76949722146069140550, 3.64784832476320460504,
            1.27045825245236838258, 2.41780725177450611770e-1,
            2.27238449892691845833e-2, 7.74545014278341407640e-4
        };
        const double const1 = 0.180625;
        const double const2 = 1.6;
        double[] d =
        {
            1.0, 2.05319162663775882187,
            1.67638483018380384940, 6.89767334985100004550e-1,
            1.48103976427480074590e-1, 1.51986665636164571966e-2,
            5.47593808499534494600e-4, 1.05075007164441684324e-9
        };
        double[] e =
        {
            6.65790464350110377720, 5.46378491116411436990,
            1.78482653991729133580, 2.96560571828504891230e-1,
            2.65321895265761230930e-2, 1.24266094738807843860e-3,
            2.71155556874348757815e-5, 2.01033439929228813265e-7
        };
        double[] f =
        {
            1.0, 5.99832206555887937690e-1,
            1.36929880922735805310e-1, 1.48753612908506148525e-2,
            7.86869131145613259100e-4, 1.84631831751005468180e-5,
            1.42151175831644588870e-7, 2.04426310338993978564e-15
        };
        double r;
        const double r8_huge = 1.0E+30;
        const double split1 = 0.425;
        const double split2 = 5.0;
        double value;

        switch (p)
        {
            case <= 0.0:
                value = -r8_huge;
                return value;
            case >= 1.0:
                value = r8_huge;
                return value;
        }

        double q = p - 0.5;

        if (Math.Abs(q) <= split1)
        {
            r = const1 - q * q;
            value = q * typeMethods.r8poly_value_horner(7, a, r) / typeMethods.r8poly_value_horner(7, b, r);
        }
        else
        {
            r = q switch
            {
                < 0.0 => p,
                _ => 1.0 - p
            };

            switch (r)
            {
                case <= 0.0:
                    value = r8_huge;
                    break;
                default:
                {
                    r = Math.Sqrt(-Math.Log(r));

                    if (r <= split2)
                    {
                        r -= const2;
                        value = typeMethods.r8poly_value_horner(7, c, r) / typeMethods.r8poly_value_horner(7, d, r);
                    }
                    else
                    {
                        r -= split2;
                        value = typeMethods.r8poly_value_horner(7, e, r) / typeMethods.r8poly_value_horner(7, f, r);
                    }

                    break;
                }
            }

            value = q switch
            {
                < 0.0 => -value,
                _ => value
            };
        }

        return value;
    }

    public static double normal_cdf(double x, double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_CDF evaluates the Normal CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, double CDF, the value of the CDF.
        //
    {
        double y = (x - mu) / sigma;

        double cdf = normal_01_cdf(y);

        return cdf;
    }

    public static double normal_cdf_inv(double cdf, double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_CDF_INV inverts the Normal CDF.
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
        //  Reference:
        //
        //    Algorithm AS 111,
        //    Applied Statistics,
        //    Volume 26, pages 118-121, 1977.
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, double NORMAL_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("NORMAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double x2 = normal_01_cdf_inv(cdf);

        double x = mu + sigma * x2;

        return x;
    }
        
    public static double normal_ms_cdf ( double x, double mu, double sigma )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MS_CDF evaluates the Normal CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, double NORMAL_MS_CDF, the value of the CDF.
        //
    {
        double y = ( x - mu ) / sigma;

        double cdf = normal_01_cdf ( y );

        return cdf;
    }

    public static double normal_ms_cdf_inv ( double cdf, double mu, double sigma )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MS_CDF_INV inverts the Normal CDF.
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
        //  Reference:
        //
        //    Algorithm AS 111,
        //    Applied Statistics,
        //    Volume 26, pages 118-121, 1977.
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, double NORMAL_MS_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("NORMAL_MS_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double x2 = normal_01_cdf_inv ( cdf );

        double x = mu + sigma * x2;

        return x;
    }

    public static double truncated_normal_ab_cdf(double x, double mu, double sigma, double a,
            double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRUNCATED_NORMAL_AB_CDF evaluates the truncated Normal CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double TRUNCATED_NORMAL_AB_CDF, the value of the CDF.
        //
    {
        double cdf;

        if (x < a)
        {
            cdf = 0.0;
        }
        else if (x <= b)
        {
            double alpha = (a - mu) / sigma;
            double beta = (b - mu) / sigma;
            double xi = (x - mu) / sigma;

            double alpha_cdf = normal_01_cdf(alpha);
            double beta_cdf = normal_01_cdf(beta);
            double xi_cdf = normal_01_cdf(xi);

            cdf = (xi_cdf - alpha_cdf) / (beta_cdf - alpha_cdf);
        }
        else
        {
            cdf = 1.0;
        }

        return cdf;
    }

    public static double normal_truncated_b_cdf(double x, double mu, double s, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_CDF evaluates the upper truncated Normal CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_B_CDF, the value of the CDF.
        //
    {
        double beta = (b - mu) / s;
        double xi = (x - mu) / s;

        double beta_cdf = normal_01_cdf(beta);
        double xi_cdf = normal_01_cdf(xi);

        double cdf = xi_cdf / beta_cdf;

        return cdf;
    }

    public static double normal_truncated_b_cdf_inv(double cdf, double mu, double s, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_CDF_INV inverts the upper truncated Normal CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_B_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("NORMAL_TRUNCATED_B_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
        }

        double beta = (b - mu) / s;

        double beta_cdf = normal_01_cdf(beta);

        double xi_cdf = beta_cdf * cdf;
        double xi = normal_01_cdf_inv(xi_cdf);

        double x = mu + s * xi;

        return x;
    }

    public static double truncated_normal_ab_cdf_inv(double cdf, double mu, double sigma, double a,
            double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRUNCATED_NORMAL_AB_CDF_INV inverts the truncated Normal CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double TRUNCATED_NORMAL_AB_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("TRUNCATED_NORMAL_AB_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double alpha = (a - mu) / sigma;
        double beta = (b - mu) / sigma;

        double alpha_cdf = normal_01_cdf(alpha);
        double beta_cdf = normal_01_cdf(beta);

        double xi_cdf = (beta_cdf - alpha_cdf) * cdf + alpha_cdf;
        double xi = normal_01_cdf_inv(xi_cdf);

        double x = mu + sigma * xi;

        return x;
    }

    public static double normal_truncated_a_cdf(double x, double mu, double s, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_CDF evaluates the lower truncated Normal CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_A_CDF, the value of the CDF.
        //
    {
        double alpha = (a - mu) / s;
        double xi = (x - mu) / s;

        double alpha_cdf = normal_01_cdf(alpha);
        double xi_cdf = normal_01_cdf(xi);

        double cdf = (xi_cdf - alpha_cdf) / (1.0 - alpha_cdf);

        return cdf;
    }

    public static double normal_truncated_a_cdf_inv(double cdf, double mu, double s, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_CDF_INV inverts the lower truncated Normal CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_A_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("NORMAL_TRUNCATED_A_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
        }

        double alpha = (a - mu) / s;

        double alpha_cdf = normal_01_cdf(alpha);

        double xi_cdf = (1.0 - alpha_cdf) * cdf + alpha_cdf;
        double xi = normal_01_cdf_inv(xi_cdf);

        double x = mu + s * xi;

        return x;
    }

    public static double normal_truncated_ab_cdf(double x, double mu, double s, double a,
            double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_CDF evaluates the truncated Normal CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double NORMAL_TRUNCATED_AB_CDF, the value of the CDF.
        //
    {
        double alpha = (a - mu) / s;
        double beta = (b - mu) / s;
        double xi = (x - mu) / s;

        double alpha_cdf = normal_01_cdf(alpha);
        double beta_cdf = normal_01_cdf(beta);
        double xi_cdf = normal_01_cdf(xi);

        double cdf = (xi_cdf - alpha_cdf) / (beta_cdf - alpha_cdf);

        return cdf;
    }

    public static double normal_truncated_ab_cdf_inv(double cdf, double mu, double s, double a,
            double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_CDF_INV inverts the truncated Normal CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double NORMAL_TRUNCATED_AB_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("NORMAL_TRUNCATED_AB_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
        }

        double alpha = (a - mu) / s;
        double beta = (b - mu) / s;

        double alpha_cdf = normal_01_cdf(alpha);
        double beta_cdf = normal_01_cdf(beta);

        double xi_cdf = (beta_cdf - alpha_cdf) * cdf + alpha_cdf;
        double xi = normal_01_cdf_inv(xi_cdf);

        double x = mu + s * xi;

        return x;
    }

    public static double truncated_normal_a_cdf_inv(double x, double mu, double sigma, double a)
    {
        return normal_truncated_a_cdf_inv(x, mu, sigma, a);
    }

    public static double truncated_normal_a_cdf(double x, double mu, double sigma, double a)
    {
        return normal_truncated_a_cdf(x, mu, sigma, a);
    }
        
    public static double truncated_normal_b_cdf ( double x, double mu, double sigma, double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRUNCATED_NORMAL_B_CDF evaluates the upper truncated Normal CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double TRUNCATED_NORMAL_B_CDF, the value of the CDF.
        //
    {
        double cdf;

        if ( x <= b )
        {
            double beta = ( b - mu ) / sigma;
            double xi = ( x - mu ) / sigma;

            double beta_cdf = normal_01_cdf ( beta );
            double xi_cdf = normal_01_cdf ( xi );

            cdf = xi_cdf / beta_cdf;
        }
        else
        {
            cdf = 1.0;
        }
  
        return cdf;
    }
        
    public static double truncated_normal_b_cdf_inv ( double cdf, double mu, double sigma, 
            double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRUNCATED_NORMAL_B_CDF_INV inverts the upper truncated Normal CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double TRUNCATED_NORMAL_B_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("TRUNCATED_NORMAL_B_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double beta = ( b - mu ) / sigma;

        double beta_cdf = normal_01_cdf ( beta );

        double xi_cdf = beta_cdf * cdf;
        double xi = normal_01_cdf_inv ( xi_cdf );

        double x = mu + sigma * xi;

        return x;
    }
}