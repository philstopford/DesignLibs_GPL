using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Gamma
{
    public static double gamma_cdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_CDF evaluates the Gamma CDF.
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
        //    Input, double X, the argument of the PDF.
        //    A <= X
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double GAMMA_CDF, the value of the CDF.
        //
    {
        double x2 = (x - a) / b;
        double p2 = c;

        double cdf = typeMethods.r8_gamma_inc(p2, x2);

        return cdf;
    }

    public static void gamma_cdf_values(ref int n_data, ref double mu, ref double sigma, ref double x,
            ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_CDF_VALUES returns some values of the Gamma CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = GammaDistribution [ mu, sigma ]
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
        //    Output, double &MU, the mean of the distribution.
        //
        //    Output, double &SIGMA, the variance of the distribution.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 12;

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
            0.4421745996289254E+00,
            0.1911531694619419E+00,
            0.6564245437845009E-01,
            0.1857593622214067E-01
        };

        double[] mu_vec =
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

        double[] sigma_vec =
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
            mu = 0.0;
            sigma = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            mu = mu_vec[n_data - 1];
            sigma = sigma_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static bool gamma_check(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_CHECK checks the parameters of the Gamma PDF.
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
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, bool GAMMA_CHECK, is true if the parameters are legal.
        //
    {
        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("GAMMA_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                Console.WriteLine("  B = " + b + "");
                return false;
        }

        switch (c)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("GAMMA_CHECK - Warning!");
                Console.WriteLine("  C <= 0.");
                Console.WriteLine("  C = " + c + "");
                return false;
            default:
                return true;
        }
    }

    public static double gamma_mean(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_MEAN returns the mean of the Gamma PDF.
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
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double GAMMA_MEAN, the mean of the PDF.
        //
    {
        double mean = a + b * c;

        return mean;
    }

    public static double gamma_pdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_PDF evaluates the Gamma PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B,C;X) = EXP ( - ( X - A ) / B ) * ( ( X - A ) / B )^(C-1)
        //      / ( B * GAMMA ( C ) )
        //
        //    GAMMA_PDF(A,B,C;X), where C is an integer, is the Erlang PDF.
        //    GAMMA_PDF(A,B,1;X) is the Exponential PDF.
        //    GAMMA_PDF(0,2,C/2;X) is the Chi Squared PDF with C degrees of freedom.
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
        //    Input, double X, the argument of the PDF.
        //    A <= X
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    A controls the location of the peak;  A is often chosen to be 0.0.
        //    B is the "scale" parameter; 0.0 < B, and is often 1.0.
        //    C is the "shape" parameter; 0.0 < C, and is often 1.0.
        //
        //    Output, double GAMMA_PDF, the value of the PDF.
        //
    {
        double pdf;
        double y;

        if (x <= a)
        {
            pdf = 0.0;
        }
        else
        {
            y = (x - a) / b;

            pdf = Math.Pow(y, c - 1.0) / (b * Helpers.Gamma(c) * Math.Exp(y));
        }

        return pdf;
    }

    public static double gamma_sample(double a, double b, double c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_SAMPLE samples the Gamma PDF.
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
        //    Original FORTRAN77 version by Ahrens and U Dieter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Joachim Ahrens, Ulrich Dieter,
        //    Generating Gamma variates by a modified rejection technique,
        //    Communications of the ACM,
        //    Volume 25, Number 1, January 1982, pages 47-54.
        //
        //    Joachim Ahrens, Ulrich Dieter,
        //    Computer methods for sampling from Gamma, Beta, Poisson and
        //    binomial distributions,
        //    Computing,
        //    Volume 12, Number 3, September 1974, pages 223-246.
        //
        //    Joachim Ahrens, Klaus-Dieter Kohrt, Ulrich Dieter,<br>
        //    Algorithm 599:
        //    Sampling from Gamma and Poisson Distributions,<br>
        //    ACM Transactions on Mathematical Software,<br>
        //    Volume 9, Number 2, June 1983, pages 255-257.
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double GAMMA_SAMPLE, a sample of the PDF.
        //
    {
        double a1 = 0.3333333;
        double a2 = -0.2500030;
        double a3 = 0.2000062;
        double a4 = -0.1662921;
        double a5 = 0.1423657;
        double a6 = -0.1367177;
        double a7 = 0.1233795;
        double e1 = 1.0;
        double e2 = 0.4999897;
        double e3 = 0.1668290;
        double e4 = 0.0407753;
        double e5 = 0.0102930;
        double euler = 2.71828182845904;
        double q1 = 0.04166669;
        double q2 = 0.02083148;
        double q3 = 0.00801191;
        double q4 = 0.00144121;
        double q5 = -0.00007388;
        double q6 = 0.00024511;
        double q7 = 0.00024240;
        double s;
        double t;
        double u;
        double x;
        switch (c)
        {
            //
            //  Allow C = 0.
            //
            case 0.0:
                x = a;
                return x;
            //
            //  C < 1.
            //
            case < 1.0:
            {
                for (;;)
                {
                    u = UniformRNG.r8_uniform_01(ref seed);
                    t = 1.0 + c / euler;
                    double p = u * t;

                    s = Exponential.exponential_01_sample(ref seed);

                    if (p < 1.0)
                    {
                        x = Math.Exp(Math.Log(p) / c);
                        if (x <= s)
                        {
                            break;
                        }
                    }
                    else
                    {
                        x = -Math.Log((t - p) / c);
                        if ((1.0 - c) * Math.Log(x) <= s)
                        {
                            break;
                        }
                    }
                }

                x = a + b * x;
                return x;
            }
            //
            default:
            {
                double s2 = c - 0.5;
                s = Math.Sqrt(c - 0.5);
                double d = Math.Sqrt(32.0) - 12.0 * Math.Sqrt(c - 0.5);

                t = Normal.normal_01_sample(ref seed);
                x = Math.Pow(Math.Sqrt(c - 0.5) + 0.5 * t, 2);

                switch (t)
                {
                    case >= 0.0:
                        x = a + b * x;
                        return x;
                }

                u = UniformRNG.r8_uniform_01(ref seed);

                if (d * u <= t * t * t)
                {
                    x = a + b * x;
                    return x;
                }

                double r = 1.0 / c;

                double q0 = ((((((
                                     q7 * r
                                     + q6) * r
                                 + q5) * r
                                + q4) * r
                               + q3) * r
                              + q2) * r
                             + q1) * r;

                double bcoef;
                double si;
                double co;
                switch (c)
                {
                    case <= 3.686:
                        bcoef = 0.463 + s - 0.178 * s2;
                        si = 1.235;
                        co = 0.195 / s - 0.079 + 0.016 * s;
                        break;
                    case <= 13.022:
                        bcoef = 1.654 + 0.0076 * s2;
                        si = 1.68 / s + 0.275;
                        co = 0.062 / s + 0.024;
                        break;
                    default:
                        bcoef = 1.77;
                        si = 0.75;
                        co = 0.1515 / s;
                        break;
                }

                double v;
                double q;
                switch (Math.Sqrt(c - 0.5) + 0.5 * t)
                {
                    case > 0.0:
                    {
                        v = 0.5 * t / s;

                        q = Math.Abs(v) switch
                        {
                            > 0.25 => q0 - s * t + 0.25 * t * t + 2.0 * s2 * Math.Log(1.0 + v),
                            _ => q0 + 0.5 * t * t *
                                ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v
                        };

                        if (Math.Log(1.0 - u) <= q)
                        {
                            x = a + b * x;
                            return x;
                        }

                        break;
                    }
                }

                for (;;)
                {
                    double e = Exponential.exponential_01_sample(ref seed);

                    u = UniformRNG.r8_uniform_01(ref seed);

                    u = 2.0 * u - 1.0;
                    t = bcoef + Math.Abs(si * e) * typeMethods.r8_sign(u);

                    switch (t)
                    {
                        case >= -0.7187449:
                        {
                            v = 0.5 * t / s;

                            q = Math.Abs(v) switch
                            {
                                > 0.25 => q0 - s * t + 0.25 * t * t + 2.0 * s2 * Math.Log(1.0 + v),
                                _ => q0 + 0.5 * t * t *
                                    ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v
                            };

                            switch (q)
                            {
                                case > 0.0:
                                {
                                    double w = q switch
                                    {
                                        > 0.5 => Math.Exp(q) - 1.0,
                                        _ => ((((e5 * q + e4) * q + e3) * q + e2) * q + e1) * q
                                    };

                                    if (co * Math.Abs(u) <= w * Math.Exp(e - 0.5 * t * t))
                                    {
                                        x = a + b * Math.Pow(s + 0.5 * t, 2);
                                        return x;
                                    }

                                    break;
                                }
                            }

                            break;
                        }
                    }
                }

                break;
            }
        }
    }

    public static double gamma_variance(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_VARIANCE returns the variance of the Gamma PDF.
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
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
    {
        double variance = b * b * c;

        return variance;
    }

    public static void gamma_inc_values(ref int n_data, ref double a, ref double x, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
        //
        //  Discussion:
        //
        //    The incomplete Gamma function is defined as:
        //
        //      Integral ( X <= T < oo ) T^(A-1) * exp(-T) dT.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Gamma[A,X]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2010
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
        //    Output, double &A, the parameter of the function.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] a_vec =
        {
            0.10E+00,
            0.10E+00,
            0.10E+00,
            0.50E+00,
            0.50E+00,
            0.50E+00,
            0.10E+01,
            0.10E+01,
            0.10E+01,
            0.11E+01,
            0.11E+01,
            0.11E+01,
            0.20E+01,
            0.20E+01,
            0.20E+01,
            0.60E+01,
            0.60E+01,
            0.11E+02,
            0.26E+02,
            0.41E+02
        };

        double[] fx_vec =
        {
            2.490302836300570E+00,
            0.8718369702247978E+00,
            0.1079213896175866E+00,
            1.238121685818417E+00,
            0.3911298052193973E+00,
            0.01444722098952533E+00,
            0.9048374180359596E+00,
            0.3678794411714423E+00,
            0.006737946999085467E+00,
            0.8827966752611692E+00,
            0.3908330082003269E+00,
            0.008051456628620993E+00,
            0.9898141728888165E+00,
            0.5578254003710746E+00,
            0.007295055724436130E+00,
            114.9574754165633E+00,
            2.440923530031405E+00,
            280854.6620274718E+00,
            8.576480283455533E+24,
            2.085031346403364E+47
        };

        double[] x_vec =
        {
            0.30E-01,
            0.30E+00,
            0.15E+01,
            0.75E-01,
            0.75E+00,
            0.35E+01,
            0.10E+00,
            0.10E+01,
            0.50E+01,
            0.10E+00,
            0.10E+01,
            0.50E+01,
            0.15E+00,
            0.15E+01,
            0.70E+01,
            0.25E+01,
            0.12E+02,
            0.16E+02,
            0.25E+02,
            0.45E+02
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
            a = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            a = a_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }
}