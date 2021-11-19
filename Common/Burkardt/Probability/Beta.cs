using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Beta
{
    public static double beta_binomial_cdf(int x, double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_BINOMIAL_CDF evaluates the Beta Binomial CDF.
        //
        //  Discussion:
        //
        //    A simple summing approach is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the argument of the CDF.
        //
        //    Input, double A, B, parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input, int C, a parameter of the PDF.
        //    0 <= C.
        //
        //    Output, double BETA_BINOMIAL_CDF, the value of the CDF.
        //
    {
        double cdf = 0;

        switch (x)
        {
            case < 0:
                cdf = 0.0;
                break;
            default:
            {
                if (x < c)
                {
                    cdf = 0.0;
                    for (int y = 0; y <= x; y++)
                    {
                        double pdf = typeMethods.r8_beta(a + y, b + (c - y))
                                     / ((c + 1)
                                        * typeMethods.r8_beta(y + 1, c - y + 1)
                                        * typeMethods.r8_beta(a, b));
                        cdf += pdf;
                    }
                }
                else if (c <= x)
                {
                    cdf = 1.0;
                }

                break;
            }
        }

        return cdf;
    }

    public static int beta_binomial_cdf_inv(double cdf, double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_BINOMIAL_CDF_INV inverts the Beta Binomial CDF.
        //
        //  Discussion:
        //
        //    A simple discrete approach is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, double A, B, parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input, int C, a parameter of the PDF.
        //    0 <= C.
        //
        //    Output, int BETA_BINOMIAL_CDF_INV, the smallest X whose cumulative
        //    density function is greater than or equal to CDF.
        //
    {
        int x;
        int y;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("BETA_BINOMIAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double cum = 0.0;

        for (y = 0; y <= c; y++)
        {
            double pdf = typeMethods.r8_beta(a + y,
                b + (c - y)) / ((c + 1)
                                * typeMethods.r8_beta(y + 1,
                                    c - y + 1) * typeMethods.r8_beta(a, b));

            cum += pdf;

            if (!(cdf <= cum))
            {
                continue;
            }

            x = y;
            return x;

        }

        x = c;

        return x;
    }

    public static bool beta_binomial_check(double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_BINOMIAL_CHECK checks the parameters of the Beta Binomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input, int C, a parameter of the PDF.
        //    0 <= C.
        //
        //    Output, bool BETA_BINOMIAL_CHECK, is TRUE if the parameters are legal.
        //
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("BETA_BINOMIAL_CHECK - Warning!");
                Console.WriteLine("  A <= 0.");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("BETA_BINOMIAL_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
        }

        switch (c)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("BETA_BINOMIAL_CHECK - Warning!");
                Console.WriteLine("  C < 0.");
                return false;
            default:
                return true;
        }
    }

    public static double beta_binomial_mean(double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_BINOMIAL_MEAN returns the mean of the Beta Binomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input, int C, a parameter of the PDF.
        //    0 <= N.
        //
        //    Output, double BETA_BINOMIAL_MEAN, the mean of the PDF.
        //
    {
        double mean = c * a / (a + b);

        return mean;
    }

    public static double beta_binomial_pdf(int x, double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_BINOMIAL_PDF evaluates the Beta Binomial PDF.
        //
        //  Discussion:
        //
        //    PDF(X)(A,B,C) = Beta(A+X,B+C-X)
        //      / ( (C+1) * Beta(X+1,C-X+1) * Beta(A,B) )  for 0 <= X <= C.
        //
        //    This PDF can be reformulated as:
        //
        //      The beta binomial probability density function for X successes
        //      out of N trials is
        //
        //      PDF2(X)( N, MU, THETA ) =
        //        C(N,X) * Product ( 0 <= R <= X - 1 ) ( MU + R * THETA )
        //               * Product ( 0 <= R <= N - X - 1 ) ( 1 - MU + R * THETA )
        //               / Product ( 0 <= R <= N - 1 )  ( 1 + R * THETA )
        //
        //      where
        //
        //        C(N,X) is the combinatorial coefficient;
        //        MU is the expectation of the underlying Beta distribution;
        //        THETA is a shape parameter.
        //
        //      A THETA value of 0 ( or A+B --> Infinity ) results in the binomial
        //      distribution:
        //
        //        PDF2(X) ( N, MU, 0 ) = C(N,X) * MU^X * ( 1 - MU )^(N-X)
        //
        //    Given A, B, C for PDF, then the equivalent PDF2 has:
        //
        //      N     = C
        //      MU    = A / ( A + B )
        //      THETA = 1 / ( A + B )
        //
        //    Given N, MU, THETA for PDF2, the equivalent PDF has:
        //
        //      A = MU / THETA
        //      B = ( 1 - MU ) / THETA
        //      C = N
        //
        //    BETA_BINOMIAL_PDF(X)(1,1,C) = UNIFORM_DISCRETE_PDF(X)(0,C-1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the argument of the PDF.
        //
        //    Input, double A, B, parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input, int C, a parameter of the PDF.
        //    0 <= C.
        //
        //    Output, double BETA_BINOMIAL_PDF, the value of the PDF.
        //
    {
        double pdf = 0;

        switch (x)
        {
            case < 0:
                pdf = 0.0;
                break;
            default:
            {
                if (x <= c)
                {
                    pdf = typeMethods.r8_beta(a + x, b + (c - x))
                          / ((c + 1)
                             * typeMethods.r8_beta(x + 1,
                                 c - x + 1) * typeMethods.r8_beta(a, b));
                }
                else if (c < x)
                {
                    pdf = 0.0;
                }

                break;
            }
        }

        return pdf;
    }

    public static int beta_binomial_sample(double a, double b, int c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_BINOMIAL_SAMPLE samples the Beta Binomial CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input, int C, a parameter of the PDF.
        //    0 <= C.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int BETA_BINOMIAL_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        int x = beta_binomial_cdf_inv(cdf, a, b, c);

        return x;
    }

    public static double beta_binomial_variance(double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_BINOMIAL_VARIANCE returns the variance of the Beta Binomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input, int C, a parameter of the PDF.
        //    0 <= C.
        //
        //    Output, double BETA_BINOMIAL_VARIANCE, the variance of the PDF.
        //
    {
        double variance = c * a * b
                          * (a + b + c)
                          / ((a + b) * (a + b) * (a + b + 1.0));

        return variance;
    }

    public static double beta_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_CDF evaluates the Beta CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2004
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
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double BETA_CDF, the value of the CDF.
        //
    {
        double cdf = x switch
        {
            <= 0.0 => 0.0,
            <= 1.0 => beta_inc(a, b, x),
            _ => 1.0
        };

        return cdf;
    }

    public static double beta_cdf_inv(double cdf, double p, double q)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_CDF_INV inverts the Beta CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 April 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    GW Cran, KJ Martin, GE Thomas,
        //    Remark AS R19 and Algorithm AS 109:
        //    A Remark on Algorithms AS 63: The Incomplete Beta Integral
        //    and AS 64: Inverse of the Incomplete Beta Integeral,
        //    Applied Statistics,
        //    Volume 26, Number 1, 1977, pages 111-114.
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the incomplete Beta
        //    function.  0 <= CDF <= 1.
        //
        //    Input, double P, Q, the parameters of the incomplete
        //    Beta function.
        //
        //    Output, double BETA_CDF_INV, the argument of the Beta CDF which 
        //    produces the value CDF.
        //
        //  Local Parameters:
        //
        //    Local, double SAE, the most negative decimal exponent
        //    which does not cause an underflow.
        //
    {
        double a;
        const double beta_log = 0.0;
        bool indx;
        double pp;
        double qq;
        const double sae = -37.0;
        double t;
        double value;

        double fpu = Math.Pow(10.0, sae);
        switch (p)
        {
            //
            //  Test for admissibility of parameters.
            //
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("BETA_CDF_INV - Fatal error!");
                Console.WriteLine("  P <= 0.0");
                value = -1.0;
                return value;
        }

        switch (q)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("BETA_CDF_INV - Fatal error!");
                Console.WriteLine("  Q <= 0.0");
                value = -1.0;
                return value;
        }

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("BETA_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0.0 or 1.0 < CDF");
                value = -1.0;
                return value;
            //
            //  Return immediately if the answer is easy.
            //
            case 0.0:
                value = 0.0;
                return value;
            case 1.0:
                value = 1.0;
                return value;
            //
            //  Change tail if necessary.
            //
            case > 0.5:
                a = 1.0 - cdf;
                pp = q;
                qq = p;
                indx = true;
                break;
            default:
                a = cdf;
                pp = p;
                qq = q;
                indx = false;
                break;
        }

        //
        //  Calculate the initial approximation.
        //
        double r = Math.Sqrt(-Math.Log(a * a));

        double y = r - (2.30753 + 0.27061 * r)
            / (1.0 + (0.99229 + 0.04481 * r) * r);

        switch (pp)
        {
            case > 1.0 when 1.0 < qq:
            {
                r = (y * y - 3.0) / 6.0;
                double s = 1.0 / (pp + pp - 1.0);
                t = 1.0 / (qq + qq - 1.0);
                double h = 2.0 / (s + t);
                double w = y * Math.Sqrt(h + r) / h - (t - s)
                    * (r + 5.0 / 6.0 - 2.0 / (3.0 * h));
                value = pp / (pp + qq * Math.Exp(w + w));
                break;
            }
            default:
            {
                r = qq + qq;
                t = 1.0 / (9.0 * qq);
                t = r * Math.Pow(1.0 - t + y * Math.Sqrt(t), 3);

                switch (t)
                {
                    case <= 0.0:
                        value = 1.0 - Math.Exp((Math.Log((1.0 - a) * qq) + beta_log) / qq);
                        break;
                    default:
                    {
                        t = (4.0 * pp + r - 2.0) / t;

                        value = t switch
                        {
                            <= 1.0 => Math.Exp((Math.Log(a * pp) + beta_log) / pp),
                            _ => 1.0 - 2.0 / (t + 1.0)
                        };

                        break;
                    }
                }

                break;
            }
        }

        //
        //  Solve for X by a modified Newton-Raphson method.
        //
        r = 1.0 - pp;
        t = 1.0 - qq;
        double yprev = 0.0;
        double sq = 1.0;
        double prev = 1.0;

        value = value switch
        {
            > 0.9999 => 0.9999,
            _ => value switch
            {
                < 0.0001 => 0.0001,
                _ => value
            }
        };

        int iex = (int)Math.Max(-5.0 / pp / pp - 1.0 / Math.Pow(a, 0.2) - 13.0, sae);

        double acu = Math.Pow(10.0, iex);

        for (;;)
        {
            y = beta_inc(pp, qq, value);

            double xin = value;
            y = (y - a) * Math.Exp(beta_log + r * Math.Log(xin) + t * Math.Log(1.0 - xin));

            prev = (y * yprev) switch
            {
                <= 0.0 => Math.Max(sq, fpu),
                _ => prev
            };

            double g = 1.0;

            double tx;
            for (;;)
            {
                for (;;)
                {
                    double adj = g * y;
                    sq = adj * adj;

                    if (sq < prev)
                    {
                        tx = value - adj;

                        if (tx is >= 0.0 and <= 1.0)
                        {
                            break;
                        }
                    }

                    g /= 3.0;
                }

                if (prev <= acu)
                {
                    value = indx switch
                    {
                        true => 1.0 - value,
                        _ => value
                    };

                    return value;
                }

                if (y * y <= acu)
                {
                    value = indx switch
                    {
                        true => 1.0 - value,
                        _ => value
                    };

                    return value;
                }

                if (tx != 0.0 && Math.Abs(tx - 1.0) > double.Epsilon)
                {
                    break;
                }

                g /= 3.0;
            }

            if (Math.Abs(tx - value) <= double.Epsilon)
            {
                break;
            }

            value = tx;
            yprev = y;
        }

        value = indx switch
        {
            true => 1.0 - value,
            _ => value
        };

        return value;
    }

    public static double beta_cdf_inv_old(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_CDF_INV_OLD inverts the Beta CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Abernathy and Smith.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Abernathy, Smith,
        //    Algorithm 724,
        //    ACM Transactions on Mathematical Software,
        //    Volume 19, Number 4, December 1993, pages 481-483.
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double BETA_CDF_INV, the argument of the CDF.
        //
    {
        const int MAXK = 20;

        double[] d = new double[MAXK * (MAXK - 1)];
        const double error = 0.0001;
        const double errapp = 0.01;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("BETA_CDF_INV_OLD - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
        }

        //
        //  Estimate the solution.
        //
        double x = a / (a + b);

        double xold = 0.0;
        int loopct = 2;

        while (errapp <= Math.Abs((x - xold) / x) && loopct != 0)
        {
            xold = x;
            loopct -= 1;
            //
            //  CDF_X = PROB { BETA(A,B) <= X }.
            //  Q = ( CDF - CDF_X ) / PDF_X.
            //
            double cdf_x = beta_cdf(x, a, b);

            double pdf_x = beta_pdf(x, a, b);

            double q = (cdf - cdf_x) / pdf_x;
            //
            //  D(N,K) = C(N,K) * Q**(N+K-1) / (N-1)!
            //
            double t = 1.0 - x;
            double s1 = q * (b - 1.0) / t;
            double s2 = q * (1.0 - a) / x;
            d[2 - 1 + 0 * MAXK] = s1 + s2;
            double tail = d[2 - 1 + 0 * MAXK] * q / 2.0;
            x = x + q + tail;

            int k = 3;

            while (error < Math.Abs(tail / x) && k <= MAXK)
            {
                //
                //  Find D(2,K-2).
                //
                s1 = q * (k - 2.0) * s1 / t;
                s2 = q * (2.0 - k) * s2 / x;
                d[2 - 1 + (k - 2) * MAXK] = s1 + s2;
                //
                //  Find D(3,K-3), D(4,K-4), D(5,K-5), ... , D(K-1,1).
                //
                int i;
                for (i = 3; i <= k - 1; i++)
                {
                    double sum2 = d[2 - 1 + 0 * MAXK] * d[i - 2 + (k - i) * MAXK];
                    double bcoeff = 1.0;

                    int j;
                    for (j = 1; j <= k - i; j++)
                    {
                        bcoeff = bcoeff * (k - i - j + 1)
                                 / j;
                        sum2 += bcoeff * d[2 - 1 + j * MAXK] * d[i - 2 + (k - i - j) * MAXK];
                    }

                    d[i - 1 + (k - i) * MAXK] = sum2 + d[i - 2 + (k - i + 1) * MAXK] / (i - 1);
                }

                //
                //  Compute D(K,0) and use it to expand the series.
                //
                d[k - 1 + 0 * MAXK] = d[2 - 1 + 0 * MAXK] * d[k - 2 + 0 * MAXK] + d[k - 2 + 1 * MAXK]
                    / (k - 1);
                tail = d[k - 1 + 0 * MAXK] * q / k;
                x += tail;
                switch (x)
                {
                    //
                    //  Check for divergence.
                    //
                    case <= 0.0:
                    case >= 1.0:
                        Console.WriteLine(" ");
                        Console.WriteLine("BETA_CDF_INV_OLD - Fatal error!");
                        Console.WriteLine("  The series has diverged.");
                        x = -1.0;
                        return x;
                    default:
                        k += 1;
                        break;
                }
            }
        }

        return x;
    }

    public static bool beta_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_CHECK checks the parameters of the Beta PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, bool BETA_CHECK, is true if the parameters are legal.
        //
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("BETA_CHECK - Warning!");
                Console.WriteLine("  A <= 0.");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("BETA_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static void beta_cdf_values(ref int n_data, ref double a, ref double b, ref double x,
            ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_CDF_VALUES returns some values of the Beta CDF.
        //
        //  Discussion:
        //
        //    The incomplete Beta function may be written
        //
        //      BETA_INC(A,B,X) = Integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
        //                      / Integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
        //
        //    Thus,
        //
        //      BETA_INC(A,B,0.0) = 0.0;
        //      BETA_INC(A,B,1.0) = 1.0
        //
        //    The incomplete Beta function is also sometimes called the
        //    "modified" Beta function, or the "normalized" Beta function
        //    or the Beta CDF (cumulative density function.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      BETA[X,A,B] / BETA[A,B]
        //
        //    The function can also be evaluated by using the Statistics package:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = BetaDistribution [ a, b ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 April 2013
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
        //    Karl Pearson,
        //    Tables of the Incomplete Beta Function,
        //    Cambridge University Press, 1968.
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
        //    Output, double &A, &B, the parameters of the function.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 45;

        double[] a_vec =
        {
            0.5E+00,
            0.5E+00,
            0.5E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            5.5E+00,
            10.0E+00,
            10.0E+00,
            10.0E+00,
            10.0E+00,
            20.0E+00,
            20.0E+00,
            20.0E+00,
            20.0E+00,
            20.0E+00,
            30.0E+00,
            30.0E+00,
            40.0E+00,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.2E+01,
            0.3E+01,
            0.4E+01,
            0.5E+01,
            1.30625,
            1.30625,
            1.30625
        };

        double[] b_vec =
        {
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            1.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            5.0E+00,
            0.5E+00,
            5.0E+00,
            5.0E+00,
            10.0E+00,
            5.0E+00,
            10.0E+00,
            10.0E+00,
            20.0E+00,
            20.0E+00,
            10.0E+00,
            10.0E+00,
            20.0E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.2E+01,
            0.3E+01,
            0.4E+01,
            0.5E+01,
            0.2E+01,
            0.2E+01,
            0.2E+01,
            0.2E+01,
            11.7562,
            11.7562,
            11.7562
        };

        double[] fx_vec =
        {
            0.6376856085851985E-01,
            0.2048327646991335E+00,
            0.1000000000000000E+01,
            0.0000000000000000E+00,
            0.5012562893380045E-02,
            0.5131670194948620E-01,
            0.2928932188134525E+00,
            0.5000000000000000E+00,
            0.2800000000000000E-01,
            0.1040000000000000E+00,
            0.2160000000000000E+00,
            0.3520000000000000E+00,
            0.5000000000000000E+00,
            0.6480000000000000E+00,
            0.7840000000000000E+00,
            0.8960000000000000E+00,
            0.9720000000000000E+00,
            0.4361908850559777E+00,
            0.1516409096347099E+00,
            0.8978271484375000E-01,
            0.1000000000000000E+01,
            0.5000000000000000E+00,
            0.4598773297575791E+00,
            0.2146816102371739E+00,
            0.9507364826957875E+00,
            0.5000000000000000E+00,
            0.8979413687105918E+00,
            0.2241297491808366E+00,
            0.7586405487192086E+00,
            0.7001783247477069E+00,
            0.5131670194948620E-01,
            0.1055728090000841E+00,
            0.1633399734659245E+00,
            0.2254033307585166E+00,
            0.3600000000000000E+00,
            0.4880000000000000E+00,
            0.5904000000000000E+00,
            0.6723200000000000E+00,
            0.2160000000000000E+00,
            0.8370000000000000E-01,
            0.3078000000000000E-01,
            0.1093500000000000E-01,
            0.918884684620518,
            0.21052977489419,
            0.1824130512500673
        };

        double[] x_vec =
        {
            0.01E+00,
            0.10E+00,
            1.00E+00,
            0.00E+00,
            0.01E+00,
            0.10E+00,
            0.50E+00,
            0.50E+00,
            0.10E+00,
            0.20E+00,
            0.30E+00,
            0.40E+00,
            0.50E+00,
            0.60E+00,
            0.70E+00,
            0.80E+00,
            0.90E+00,
            0.50E+00,
            0.90E+00,
            0.50E+00,
            1.00E+00,
            0.50E+00,
            0.80E+00,
            0.60E+00,
            0.80E+00,
            0.50E+00,
            0.60E+00,
            0.70E+00,
            0.80E+00,
            0.70E+00,
            0.10E+00,
            0.20E+00,
            0.30E+00,
            0.40E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.30E+00,
            0.30E+00,
            0.30E+00,
            0.30E+00,
            0.225609,
            0.0335568,
            0.0295222
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
            b = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            a = a_vec[n_data - 1];
            b = b_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static double beta_inc(double a, double b, double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_INC returns the value of the incomplete Beta function.
        //
        //  Discussion:
        //
        //    This calculation requires an iteration.  In some cases, the iteration
        //    may not converge rapidly, or may become inaccurate.
        //
        //    BETA_INC(A,B,X)
        //
        //      =   Integral ( 0 <= T <= X ) T^(A-1) (1-T)^(B-1) dT
        //        / Integral ( 0 <= T <= 1 ) T^(A-1) (1-T)^(B-1) dT
        //
        //      =   Integral ( 0 <= T <= X ) T^(A-1) (1-T)^(B-1) dT
        //        / BETA(A,B)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Majumder, Bhattacharjee.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Majumder, Bhattacharjee,
        //    Algorithm AS63,
        //    Applied Statistics,
        //    1973, volume 22, number 3.
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the function.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input, double X, the argument of the function.
        //    Normally, 0.0 <= X <= 1.0.
        //
        //    Output, double BETA_INC, the value of the function.
        //
    {
        double cx;
        const int it_max = 1000;
        bool indx;
        double pp;
        double qq;
        const double tol = 1.0E-07;
        double value;
        double xx;

        switch (a)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("BETA_INC - Fatal error!");
                Console.WriteLine("  A <= 0.");
                return 1.0;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("BETA_INC - Fatal error!");
                Console.WriteLine("  B <= 0.");
                return 1.0;
        }

        switch (x)
        {
            case <= 0.0:
                value = 0.0;
                return value;
            case >= 1.0:
                value = 1.0;
                return value;
        }

        //
        //  Change tail if necessary and determine S.
        //
        double psq = a + b;

        if (a < (a + b) * x)
        {
            xx = 1.0 - x;
            cx = x;
            pp = b;
            qq = a;
            indx = true;
        }
        else
        {
            xx = x;
            cx = 1.0 - x;
            pp = a;
            qq = b;
            indx = false;
        }

        double term = 1.0;
        int i = 1;
        value = 1.0;

        int ns = (int) (qq + cx * (a + b));
        //
        //  Use Soper's reduction formulas.
        //
        double rx = xx / cx;

        double temp = qq - i;
        rx = ns switch
        {
            0 => xx,
            _ => rx
        };

        int it = 0;

        for (;;)
        {
            it += 1;

            if (it_max < it)
            {
                Console.WriteLine("");
                Console.WriteLine("BETA_INC - Fatal error!");
                Console.WriteLine("  Maximum number of iterations exceeded!");
                Console.WriteLine("  IT_MAX = " + it_max + "");
                return 1.0;
            }

            term = term * temp * rx / (pp + i);
            value += term;
            temp = Math.Abs(term);

            if (temp <= tol && temp <= tol * value)
            {
                break;
            }

            i += 1;
            ns -= 1;

            switch (ns)
            {
                case >= 0:
                {
                    temp = qq - i;
                    rx = ns switch
                    {
                        0 => xx,
                        _ => rx
                    };

                    break;
                }
                default:
                    temp = psq;
                    psq += 1.0;
                    break;
            }
        }

        //
        //  Finish calculation.
        //
        value = value * Math.Exp(pp * Math.Log(xx)
                                 + (qq - 1.0) * Math.Log(cx)) / (typeMethods.r8_beta(a, b) * pp);

        switch (indx)
        {
            case true:
                value = 1.0 - value;
                break;
        }

        return value;
    }

    public static void beta_inc_values(ref int n_data, ref double a, ref double b, ref double x,
            ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_INC_VALUES returns some values of the incomplete Beta function.
        //
        //  Discussion:
        //
        //    The incomplete Beta function may be written
        //
        //      BETA_INC(A,B,X) = Integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
        //                      / Integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
        //
        //    Thus,
        //
        //      BETA_INC(A,B,0.0) = 0.0;
        //      BETA_INC(A,B,1.0) = 1.0
        //
        //    The incomplete Beta function is also sometimes called the
        //    "modified" Beta function, or the "normalized" Beta function
        //    or the Beta CDF (cumulative density function.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      BETA[X,A,B] / BETA[A,B]
        //
        //    The function can also be evaluated by using the Statistics package:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = BetaDistribution [ a, b ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 April 2013
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
        //    Karl Pearson,
        //    Tables of the Incomplete Beta Function,
        //    Cambridge University Press, 1968.
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
        //    Output, double &A, &B, the parameters of the function.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 45;

        double[] a_vec =
        {
            0.5E+00,
            0.5E+00,
            0.5E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            5.5E+00,
            10.0E+00,
            10.0E+00,
            10.0E+00,
            10.0E+00,
            20.0E+00,
            20.0E+00,
            20.0E+00,
            20.0E+00,
            20.0E+00,
            30.0E+00,
            30.0E+00,
            40.0E+00,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.2E+01,
            0.3E+01,
            0.4E+01,
            0.5E+01,
            1.30625,
            1.30625,
            1.30625
        };

        double[] b_vec =
        {
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            1.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            2.0E+00,
            5.0E+00,
            0.5E+00,
            5.0E+00,
            5.0E+00,
            10.0E+00,
            5.0E+00,
            10.0E+00,
            10.0E+00,
            20.0E+00,
            20.0E+00,
            10.0E+00,
            10.0E+00,
            20.0E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.2E+01,
            0.3E+01,
            0.4E+01,
            0.5E+01,
            0.2E+01,
            0.2E+01,
            0.2E+01,
            0.2E+01,
            11.7562,
            11.7562,
            11.7562
        };

        double[] fx_vec =
        {
            0.6376856085851985E-01,
            0.2048327646991335E+00,
            0.1000000000000000E+01,
            0.0000000000000000E+00,
            0.5012562893380045E-02,
            0.5131670194948620E-01,
            0.2928932188134525E+00,
            0.5000000000000000E+00,
            0.2800000000000000E-01,
            0.1040000000000000E+00,
            0.2160000000000000E+00,
            0.3520000000000000E+00,
            0.5000000000000000E+00,
            0.6480000000000000E+00,
            0.7840000000000000E+00,
            0.8960000000000000E+00,
            0.9720000000000000E+00,
            0.4361908850559777E+00,
            0.1516409096347099E+00,
            0.8978271484375000E-01,
            0.1000000000000000E+01,
            0.5000000000000000E+00,
            0.4598773297575791E+00,
            0.2146816102371739E+00,
            0.9507364826957875E+00,
            0.5000000000000000E+00,
            0.8979413687105918E+00,
            0.2241297491808366E+00,
            0.7586405487192086E+00,
            0.7001783247477069E+00,
            0.5131670194948620E-01,
            0.1055728090000841E+00,
            0.1633399734659245E+00,
            0.2254033307585166E+00,
            0.3600000000000000E+00,
            0.4880000000000000E+00,
            0.5904000000000000E+00,
            0.6723200000000000E+00,
            0.2160000000000000E+00,
            0.8370000000000000E-01,
            0.3078000000000000E-01,
            0.1093500000000000E-01,
            0.918884684620518,
            0.21052977489419,
            0.1824130512500673
        };

        double[] x_vec =
        {
            0.01E+00,
            0.10E+00,
            1.00E+00,
            0.00E+00,
            0.01E+00,
            0.10E+00,
            0.50E+00,
            0.50E+00,
            0.10E+00,
            0.20E+00,
            0.30E+00,
            0.40E+00,
            0.50E+00,
            0.60E+00,
            0.70E+00,
            0.80E+00,
            0.90E+00,
            0.50E+00,
            0.90E+00,
            0.50E+00,
            1.00E+00,
            0.50E+00,
            0.80E+00,
            0.60E+00,
            0.80E+00,
            0.50E+00,
            0.60E+00,
            0.70E+00,
            0.80E+00,
            0.70E+00,
            0.10E+00,
            0.20E+00,
            0.30E+00,
            0.40E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.30E+00,
            0.30E+00,
            0.30E+00,
            0.30E+00,
            0.225609,
            0.0335568,
            0.0295222
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
            b = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            a = a_vec[n_data - 1];
            b = b_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static double beta_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_MEAN returns the mean of the Beta PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double BETA_MEAN, the mean of the PDF.
        //
    {
        double mean = a / (a + b);

        return mean;
    }

    public static double beta_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_PDF evaluates the Beta PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = X^(A-1) * (1-X)^(B-1) / BETA(A,B).
        //
        //    A = B = 1 yields the Uniform distribution on [0,1].
        //    A = B = 1/2 yields the Arcsin distribution.
        //        B = 1 yields the power function distribution.
        //    A = B -> Infinity tends to the Normal distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    0.0 <= X <= 1.0.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double BETA_PDF, the value of the PDF.
        //
    {
        double pdf;

        switch (x)
        {
            case < 0.0:
            case > 1.0:
                pdf = 0.0;
                break;
            default:
                pdf = Math.Pow(x, a - 1.0) * Math.Pow(1.0 - x, b - 1.0)
                      / typeMethods.r8_beta(a, b);
                break;
        }

        return pdf;
    }

    public static double beta_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_SAMPLE samples the Beta PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Kennedy, James Gentle,
        //    Algorithm BN,
        //    Statistical Computing,
        //    Dekker, 1980.
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double BETA_SAMPLE, a sample of the PDF.
        //
    {
        double x;

        double mu = (a - 1.0) / (a + b - 2.0);
        double stdev = 0.5 / Math.Sqrt(a + b - 2.0);

        for (;;)
        {
            double y = Normal.normal_01_sample(ref seed);

            x = mu + stdev * y;

            switch (x)
            {
                case < 0.0:
                case > 1.0:
                    continue;
            }

            double u = UniformRNG.r8_uniform_01(ref seed);

            double test = (a - 1.0) * Math.Log(x / (a - 1.0))
                          + (b - 1.0) * Math.Log((1.0 - x) / (b - 1.0))
                          + (a + b - 2.0) * Math.Log(a + b - 2.0) + 0.5 * y * y;

            if (Math.Log(u) <= test)
            {
                break;
            }

        }

        return x;
    }

    public static void beta_values(ref int n_data, ref double x, ref double y, ref double fxy )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_VALUES returns some values of the Beta function.
        //
        //  Discussion:
        //
        //    Beta(X,Y) = ( Gamma(X) * Gamma(Y) ) / Gamma(X+Y)
        //
        //    Both X and Y must be greater than 0.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Beta[X,Y]
        //
        //  Properties:
        //
        //    Beta(X,Y) = Beta(Y,X).
        //    Beta(X,Y) = Integral ( 0 <= T <= 1 ) T^(X-1) (1-T)^(Y-1) dT.
        //    Beta(X,Y) = Gamma(X) * Gamma(Y) / Gamma(X+Y)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2004
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
        //    Output, double &X, &Y, the arguments of the function.
        //
        //    Output, double &FXY, the value of the function.
        //
    {
        const int N_MAX = 17;

        double[] b_vec =
        {
            0.5000000000000000E+01,
            0.2500000000000000E+01,
            0.1666666666666667E+01,
            0.1250000000000000E+01,
            0.5000000000000000E+01,
            0.2500000000000000E+01,
            0.1000000000000000E+01,
            0.1666666666666667E+00,
            0.3333333333333333E-01,
            0.7142857142857143E-02,
            0.1587301587301587E-02,
            0.2380952380952381E-01,
            0.5952380952380952E-02,
            0.1984126984126984E-02,
            0.7936507936507937E-03,
            0.3607503607503608E-03,
            0.8325008325008325E-04
        };

        double[] x_vec =
        {
            0.2E+00,
            0.4E+00,
            0.6E+00,
            0.8E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            2.0E+00,
            3.0E+00,
            4.0E+00,
            5.0E+00,
            6.0E+00,
            6.0E+00,
            6.0E+00,
            6.0E+00,
            6.0E+00,
            7.0E+00
        };

        double[] y_vec =
        {
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            0.2E+00,
            0.4E+00,
            1.0E+00,
            2.0E+00,
            3.0E+00,
            4.0E+00,
            5.0E+00,
            2.0E+00,
            3.0E+00,
            4.0E+00,
            5.0E+00,
            6.0E+00,
            7.0E+00
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
            x = 0.0;
            y = 0.0;
            fxy = 0.0;
        }
        else
        {
            x = x_vec[n_data - 1];
            y = y_vec[n_data - 1];
            fxy = b_vec[n_data - 1];
        }
    }

    public static double beta_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_VARIANCE returns the variance of the Beta PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double BETA_VARIANCE, the variance of the PDF.
        //
    {
        double variance = a * b / ((a + b) * (a + b) * (1.0 + a + b));

        return variance;
    }
}