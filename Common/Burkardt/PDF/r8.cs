﻿using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.PDFLib
{
    public static partial class PDF
    {
        public static double r8_beta_pdf(double alpha, double beta, double rval)
//****************************************************************************80
//
//  Purpose:
//
//    R8_BETA_PDF evaluates the PDF of a beta distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2015
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double ALPHA, BETA, shape parameters.
//    0.0 < ALPHA, BETA.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_BETA_PDF, the value of the PDF at RVAL.
//
        {
            double temp;
            double value;

            if (alpha <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_BETA_PDF - Fatal error!");
                Console.WriteLine("  Parameter ALPHA is not positive.");
                return (1);
            }

            if (beta <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_BETA_PDF- Fatal error!");
                Console.WriteLine("  Parameter BETA is not positive.");
                return (1);
            }

            if (rval <= 0.0 || 1.0 <= rval)
            {
                value = 0.0;
            }
            else
            {
                temp = typeMethods.r8_gamma_log(alpha + beta) - typeMethods.r8_gamma_log(alpha)
                                                              - typeMethods.r8_gamma_log(beta);

                value = Math.Exp(temp) * Math.Pow(rval, alpha - 1.0)
                                       * Math.Pow(1.0 - rval, beta - 1.0);
            }

            return value;
        }

        public static double r8_beta_sample(double aa, double bb)
//****************************************************************************80
//
//  Purpose:
//
//    R8_BETA_SAMPLE generates a beta random deviate.
//
//  Discussion:
//
//    This procedure returns a single random deviate from the beta distribution
//    with parameters A and B.  The density is
//
//      x^(a-1) * (1-x)^(b-1) / Beta(a,b) for 0 < x < 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Russell Cheng,
//    Generating Beta Variates with Nonintegral Shape Parameters,
//    Communications of the ACM,
//    Volume 21, Number 4, April 1978, pages 317-322.
//
//  Parameters:
//
//    Input, double AA, the first parameter of the beta distribution.
//    0.0 < AA.
//
//    Input, double BB, the second parameter of the beta distribution.
//    0.0 < BB.
//
//    Output, double R8_BETA_SAMPLE, a beta random variate.
//
        {
            double a;
            double alpha;
            double b;
            double beta;
            double delta;
            double gamma;
            double k1;
            double k2;
            const double log4 = 1.3862943611198906188;
            const double log5 = 1.6094379124341003746;
            double r;
            double s;
            double t;
            double u1;
            double u2;
            double v;
            double value;
            double w;
            double y;
            double z;

            if (aa <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_BETA_SAMPLE - Fatal error!");
                Console.WriteLine("  AA <= 0.0");
                return (1);
            }

            if (bb <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_BETA_SAMPLE - Fatal error!");
                Console.WriteLine("  BB <= 0.0");
                return (1);
            }

//
//  Algorithm BB
//
            if (1.0 < aa && 1.0 < bb)
            {
                if (aa < bb)
                {
                    a = aa;
                    b = bb;
                }
                else
                {
                    a = bb;
                    b = aa;
                }

                alpha = a + b;
                beta = Math.Sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
                gamma = a + 1.0 / beta;

                for (;;)
                {
                    u1 = r8_uniform_01_sample();
                    u2 = r8_uniform_01_sample();
                    v = beta * Math.Log(u1 / (1.0 - u1));
                    w = a * Math.Exp(v);

                    z = u1 * u1 * u2;
                    r = gamma * v - log4;
                    s = a + r - w;

                    if (5.0 * z <= s + 1.0 + log5)
                    {
                        break;
                    }

                    t = Math.Log(z);
                    if (t <= s)
                    {
                        break;
                    }

                    if (t <= (r + alpha * Math.Log(alpha / (b + w))))
                    {
                        break;
                    }
                }
            }
//
//  Algorithm BC
//
            else
            {
                if (aa < bb)
                {
                    a = bb;
                    b = aa;
                }
                else
                {
                    a = aa;
                    b = bb;
                }

                alpha = a + b;
                beta = 1.0 / b;
                delta = 1.0 + a - b;
                k1 = delta * (1.0 / 72.0 + b / 24.0)
                     / (a / b - 7.0 / 9.0);
                k2 = 0.25 + (0.5 + 0.25 / delta) * b;

                for (;;)
                {
                    u1 = r8_uniform_01_sample();
                    u2 = r8_uniform_01_sample();

                    if (u1 < 0.5)
                    {
                        y = u1 * u2;
                        z = u1 * y;

                        if (k1 <= 0.25 * u2 + z - y)
                        {
                            continue;
                        }
                    }
                    else
                    {
                        z = u1 * u1 * u2;

                        if (z <= 0.25)
                        {
                            v = beta * Math.Log(u1 / (1.0 - u1));
                            w = a * Math.Exp(v);

                            if (aa == a)
                            {
                                value = w / (b + w);
                            }
                            else
                            {
                                value = b / (b + w);
                            }

                            return value;
                        }

                        if (k2 < z)
                        {
                            continue;
                        }
                    }

                    v = beta * Math.Log(u1 / (1.0 - u1));
                    w = a * Math.Exp(v);

                    if (Math.Log(z) <= alpha * (Math.Log(alpha / (b + w)) + v) - log4)
                    {
                        break;
                    }
                }
            }

            if (aa == a)
            {
                value = w / (b + w);
            }
            else
            {
                value = b / (b + w);
            }

            return value;
        }

        public static double r8_chi_pdf(double df, double rval)
//****************************************************************************80
//
//  Purpose:
//
//    R8_CHI_PDF evaluates the PDF of a chi-squared distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C  by John Burkardt.
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_CHI_PDF, the value of the PDF at RVAL.
//
        {
            double temp1;
            double temp2;
            double value;

            if (df <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_CHI_PDF - Fatal error!");
                Console.WriteLine("  Degrees of freedom must be positive.");
                return (1);
            }

            if (rval <= 0.0)
            {
                value = 0.0;
            }
            else
            {
                temp2 = df * 0.5;

                temp1 = (temp2 - 1.0) * Math.Log(rval) - 0.5 * rval
                                                       - temp2 * Math.Log(2.0) - typeMethods.r8_gamma_log(temp2);

                value = Math.Exp(temp1);
            }

            return value;
        }

        public static double r8_chi_sample(double df)
//****************************************************************************80
//
//  Purpose:
//
//    R8_CHI_SAMPLE generates a Chi-Square random deviate.
//
//  Discussion:
//
//    This procedure generates a random deviate from the chi square distribution
//    with DF degrees of freedom random variable.
//
//    The algorithm exploits the relation between chisquare and gamma.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Output, double R8_CHI_SAMPLE, a random deviate from the distribution.
//
        {
            double arg1;
            double arg2;
            double value;

            if (df <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_CHI_SAMPLE - Fatal error!");
                Console.WriteLine("  DF <= 0.");
                Console.WriteLine("  Value of DF: " + df + "");
                return (1);
            }

            arg1 = 1.0;
            arg2 = df / 2.0;

            value = 2.0 * r8_gamma_sample(arg1, arg2);

            return value;
        }

        public static double r8_choose(int n, int k)
//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, double R8_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
        {
            int i;
            int mn;
            int mx;
            double value;

            if (k < n - k)
            {
                mn = k;
                mx = n - k;
            }
            else
            {
                mn = n - k;
                mx = k;
            }

            if (mn < 0)
            {
                value = 0.0;
            }
            else if (mn == 0)
            {
                value = 1.0;
            }
            else
            {
                value = (double) (mx + 1);

                for (i = 2; i <= mn; i++)
                {
                    value = (value * (double) (mx + i)) / (double) i;
                }
            }

            return value;
        }

        public static double r8_exponential_pdf(double beta, double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_EXPONENTIAL_PDF evaluates the PDF of an exponential distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double BETA, the scale value.
//    0.0 < BETA.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_EXPONENTIAL_PDF, the value of the PDF at RVAL.
//
        {
            double value;

            if (beta <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_EXPONENTIAL_PDF - Fatal error!");
                Console.WriteLine("  BETA parameter must be positive.");
                return (1);
            }

            if (rval < 0.0)
            {
                value = 0.0;
            }
            else
            {
                value = Math.Exp(-rval / beta) / beta;
            }

            return value;
        }

        public static double r8_exponential_sample(double lambda)

//****************************************************************************80
//
//  Purpose:
//
//    R8_EXPONENTIAL_SAMPLE samples the exponential PDF.
//
//  Discussion:
//
//    Note that the parameter LAMBDA is a multiplier.  In some formulations,
//    it is used as a divisor instead.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double LAMBDA, the parameter of the PDF.
//
//    Output, double R8_EXPONENTIAL_SAMPLE, a sample of the PDF.
//
        {
            double r;
            double value;

            r = r8_uniform_01_sample();

            value = -Math.Log(r) * lambda;

            return value;
        }

        public static double r8_exponential_01_pdf(double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_EXPONENTIAL_01_PDF: PDF of the standard exponential distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_EXPONENTIAL_01_PDF, the value of the PDF at RVAL.
//
        {
            double value;

            if (rval < 0.0)
            {
                value = 0.0;
            }
            else
            {
                value = Math.Exp(-rval);
            }

            return value;
        }

        public static double r8_exponential_01_sample()

//****************************************************************************80
//
//  Purpose:
//
//    R8_EXPONENTIAL_01_SAMPLE samples the standard exponential PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EXPONENTIAL_01_SAMPLE, a sample of the PDF.
//
        {
            double r;
            double value;

            r = r8_uniform_01_sample();

            value = -Math.Log(r);

            return value;
        }

        public static double r8_gamma_pdf(double beta, double alpha, double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_PDF evaluates the PDF of a gamma distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double BETA, the rate parameter.
//    0.0 < BETA.
//
//    Input, double ALPHA, the shape parameter.
//    0.0 < ALPHA.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_GAMMA_PDF, the value of the PDF at RVAL.
//
        {
            double temp;
            double value;

            if (alpha <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GAMMA_PDF - Fatal error!");
                Console.WriteLine("  Parameter ALPHA is not positive.");
                return (1);
            }

            if (beta <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GAMMA_PDF - Fatal error!");
                Console.WriteLine("  Parameter BETA is not positive.");
                return (1);
            }

            if (rval <= 0.0)
            {
                value = 0.0;
            }
            else
            {
                temp = alpha * Math.Log(beta) + (alpha - 1.0) * Math.Log(rval)
                       - beta * rval - typeMethods.r8_gamma_log(alpha);

                value = Math.Exp(temp);
            }

            return value;
        }

        public static double r8_gamma_sample(double a, double r)

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_SAMPLE generates a Gamma random deviate.
//
//  Discussion:
//
//    This procedure generates random deviates from the gamma distribution whose
//    density is (A^R)/Gamma(R) * X^(R-1) * Exp(-A*X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joachim Ahrens, Ulrich Dieter,
//    Generating Gamma Variates by a Modified Rejection Technique,
//    Communications of the ACM,
//    Volume 25, Number 1, January 1982, pages 47-54.
//
//    Joachim Ahrens, Ulrich Dieter,
//    Computer Methods for Sampling from Gamma, Beta, Poisson and
//    Binomial Distributions,
//    Computing,
//    Volume 12, Number 3, September 1974, pages 223-246.
//
//  Parameters:
//
//    Input, double A, the location parameter.
//    A nonzero.
//
//    Input, double R, the shape parameter.
//    0.0 < R.
//
//    Output, double R8_GAMMA_SAMPLE, a random deviate from the distribution.
//
        {
            double value;

            value = r8_gamma_01_sample(r) / a;

            return value;
        }

        public static double r8_gamma_01_pdf(double alpha, double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_01_PDF evaluates the PDF of a standard gamma distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double ALPHA, the shape parameter.
//    0.0 < ALPHA.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_GAMMA_01_PDF, the value of the PDF at RVAL.
//
        {
            double temp;
            double value;

            if (alpha <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GAMMA_01_PDF - Fatal error!");
                Console.WriteLine("  Parameter ALPHA is not positive.");
                return (1);
            }

            if (rval <= 0.0)
            {
                value = 0.0;
            }
            else
            {
                temp = (alpha - 1.0) * Math.Log(rval) - rval - typeMethods.r8_gamma_log(alpha);

                value = Math.Exp(temp);
            }

            return value;
        }

        public static double r8_gamma_01_sample(double a)

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_01_SAMPLE samples the standard Gamma distribution.
//
//  Discussion:
//
//    This procedure corresponds to algorithm GD in the reference.
//
//    pdf ( a; x ) = 1/gamma(a) * x^(a-1) * exp ( - x )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joachim Ahrens, Ulrich Dieter,
//    Generating Gamma Variates by a Modified Rejection Technique,
//    Communications of the ACM,
//    Volume 25, Number 1, January 1982, pages 47-54.
//
//  Parameters:
//
//    Input, double A, the shape parameter.
//    0.0 < A.
//
//    Output, double R8_GAMMA_01_SAMPLE, a random deviate from the distribution.
//
        {
            double a1 = 0.3333333;
            double a2 = -0.2500030;
            double a3 = 0.2000062;
            double a4 = -0.1662921;
            double a5 = 0.1423657;
            double a6 = -0.1367177;
            double a7 = 0.1233795;
            double b;
            double c;
            double d;
            double e;
            double e1 = 1.0;
            double e2 = 0.4999897;
            double e3 = 0.1668290;
            double e4 = 0.0407753;
            double e5 = 0.0102930;
            double p;
            double q;
            double q0;
            double q1 = 0.04166669;
            double q2 = 0.02083148;
            double q3 = 0.00801191;
            double q4 = 0.00144121;
            double q5 = -0.00007388;
            double q6 = 0.00024511;
            double q7 = 0.00024240;
            double r;
            double s;
            double s2;
            double si;
            double sqrt32 = 5.6568542494923801952;
            double t;
            double u;
            double v;
            double value = 0;
            double w;
            double x;

            if (1.0 <= a)
            {
                s2 = a - 0.5;
                s = Math.Sqrt(s2);
                d = sqrt32 - 12.0 * s;
//
//  Immediate acceptance.
//
                t = r8_normal_01_sample();
                x = s + 0.5 * t;
                value = x * x;

                if (0.0 <= t)
                {
                    return value;
                }

//
//  Squeeze acceptance.
//
                u = r8_uniform_01_sample();
                if (d * u <= t * t * t)
                {
                    return value;
                }

                r = 1.0 / a;
                q0 = ((((((
                              q7 * r
                              + q6) * r
                          + q5) * r
                         + q4) * r
                        + q3) * r
                       + q2) * r
                      + q1) * r;
//
//  Approximation depending on size of parameter A.
//
                if (13.022 < a)
                {
                    b = 1.77;
                    si = 0.75;
                    c = 0.1515 / s;
                }
                else if (3.686 < a)
                {
                    b = 1.654 + 0.0076 * s2;
                    si = 1.68 / s + 0.275;
                    c = 0.062 / s + 0.024;
                }
                else
                {
                    b = 0.463 + s + 0.178 * s2;
                    si = 1.235;
                    c = 0.195 / s - 0.079 + 0.16 * s;
                }

//
//  Quotient test.
//
                if (0.0 < x)
                {
                    v = 0.5 * t / s;

                    if (0.25 < Math.Abs(v))
                    {
                        q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * Math.Log(1.0 + v);
                    }
                    else
                    {
                        q = q0 + 0.5 * t * t * ((((((
                                                        a7 * v
                                                        + a6) * v
                                                    + a5) * v
                                                   + a4) * v
                                                  + a3) * v
                                                 + a2) * v
                                                + a1) * v;
                    }

                    if (Math.Log(1.0 - u) <= q)
                    {
                        return value;
                    }
                }

                for (;;)
                {
                    e = r8_exponential_01_sample();
                    u = 2.0 * r8_uniform_01_sample() - 1.0;

                    if (0.0 <= u)
                    {
                        t = b + Math.Abs(si * e);
                    }
                    else
                    {
                        t = b - Math.Abs(si * e);
                    }

//
//  Possible rejection.
//
                    if (t < -0.7187449)
                    {
                        continue;
                    }

//
//  Calculate V and quotient Q.
//
                    v = 0.5 * t / s;

                    if (0.25 < Math.Abs(v))
                    {
                        q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * Math.Log(1.0 + v);
                    }
                    else
                    {
                        q = q0 + 0.5 * t * t * ((((((
                                                        a7 * v
                                                        + a6) * v
                                                    + a5) * v
                                                   + a4) * v
                                                  + a3) * v
                                                 + a2) * v
                                                + a1) * v;
                    }

//
//  Hat acceptance.
//
                    if (q <= 0.0)
                    {
                        continue;
                    }

                    if (0.5 < q)
                    {
                        w = Math.Exp(q) - 1.0;
                    }
                    else
                    {
                        w = ((((
                                   e5 * q
                                   + e4) * q
                               + e3) * q
                              + e2) * q
                             + e1) * q;
                    }

//
//  May have to sample again.
//
                    if (c * Math.Abs(u) <= w * Math.Exp(e - 0.5 * t * t))
                    {
                        break;
                    }
                }

                x = s + 0.5 * t;
                value = x * x;
            }
//
//  Method for A < 1.
//
            else if (a < 1.0)
            {
                b = 1.0 + 0.3678794 * a;

                for (;;)
                {
                    p = b * r8_uniform_01_sample();

                    if (p < 1.0)
                    {
                        value = Math.Exp(Math.Log(p) / a);
                        if (value <= r8_exponential_01_sample())
                        {
                            break;
                        }
                    }
                    else
                    {
                        value = -Math.Log((b - p) / a);
                        if ((1.0 - a) * Math.Log(value) <= r8_exponential_01_sample())
                        {
                            break;
                        }
                    }
                }
            }

            return value;
        }

        public static double r8_invchi_pdf(double df, double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_INVCHI_PDF evaluates the PDF of an inverse chi-squared distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_INVCHI_PDF, the value of the PDF at RVAL.
//
        {
            double temp1;
            double temp2;
            double value;

            if (df <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_INVCHI_PDF - Fatal error!");
                Console.WriteLine("  Degrees of freedom must be positive.");
                return (1);
            }

            if (rval <= 0.0)
            {
                value = 0.0;
            }
            else
            {
                temp2 = df * 0.5;

                temp1 = -temp2 * Math.Log(2.0) - (temp2 + 1.0) * Math.Log(rval)
                                               - 0.5 / rval - typeMethods.r8_gamma_log(temp2);

                value = Math.Exp(temp1);
            }

            return value;
        }

        public static double r8_invchi_sample(double df)

//****************************************************************************80
//
//  Purpose:
//
//    R8_INVCHI_SAMPLE samples the inverse chi-squared distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Output, double R8_INVCHI_SAMPLE, the sample of the PDF.
//
        {
            double a;
            double b;
            double value;

            a = 0.5;
            b = 0.5 * df;
            value = r8_gamma_sample(a, b);

            if (value != 0.0)
            {
                value = 1.0 / value;
            }

            return value;
        }

        public static double r8_invgam_pdf(double beta, double alpha, double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_INVGAM_PDF evaluates the PDF of an inverse gamma distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double BETA, the rate parameter.
//    0.0 < BETA.
//
//    Input, double ALPHA, the shape parameter.
//    0.0 < ALPHA.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_INVGAM_PDF, the value of the PDF at RVAL.
//
        {
            double temp;
            double value;

            if (alpha <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_INVGAM_PDF - Fatal error!");
                Console.WriteLine("  Parameter ALPHA is not positive.");
                return (1);
            }

            if (beta <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_INVGAM_PDF - Fatal error!");
                Console.WriteLine("  Parameter BETA is not positive.");
                return (1);
            }

            if (rval <= 0.0)
            {
                value = 0.0;
            }
            else
            {
                temp = alpha * Math.Log(beta) - (alpha + 1.0) * Math.Log(rval)
                                              - beta / rval - typeMethods.r8_gamma_log(alpha);

                value = Math.Exp(temp);
            }

            return value;
        }

        public static double r8_invgam_sample(double beta, double alpha)

//****************************************************************************80
//
//  Purpose:
//
//    R8_INVGAM_SAMPLE samples an inverse gamma distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double BETA, the rate parameter.
//    0.0 < BETA.
//
//    Input, double ALPHA, the shape parameter.
//    0.0 < ALPHA.
//
//    Output, double R8_INVGAM_SAMPLE, a sample of the PDF.
//
        {
            double value;

            value = r8_gamma_sample(beta, alpha);
            if (value != 0.0)
            {
                value = 1.0 / value;
            }


            return value;
        }


        public static double r8_normal_pdf(double av, double sd, double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_PDF evaluates the PDF of a normal distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double AV, the mean value.
//
//    Input, double SD, the standard deviation.
//    0.0 < SD.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_NORMAL_PDF, the value of the PDF at RVAL.
//
        {
            double pi = 3.141592653589793;
            double rtemp;
            double value;

            if (sd <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_NORMAL_PDF - Fatal error!");
                Console.WriteLine("  Standard deviation must be positive.");
                return (1);
            }

            rtemp = (rval - av) * (rval - av) * 0.5 / (sd * sd);

            value = Math.Exp(-rtemp) / sd / Math.Sqrt(2.0 * pi);

            return value;
        }

        public static double r8_normal_sample(double av, double sd)

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_SAMPLE generates a normal random deviate.
//
//  Discussion:
//
//    This procedure generates a single random deviate from a normal 
//    distribution with mean AV, and standard deviation SD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joachim Ahrens, Ulrich Dieter,
//    Extensions of Forsythe's Method for Random
//    Sampling from the Normal Distribution,
//    Mathematics of Computation,
//    Volume 27, Number 124, October 1973, page 927-937.
//
//  Parameters:
//
//    Input, double AV, the mean.
//
//    Input, double SD, the standard deviation.
//
//    Output, double R8_NORMAL_SAMPLE, a random deviate from the distribution.
//
        {
            double value;

            value = sd * r8_normal_01_sample() + av;

            return value;
        }

        public static double r8_normal_01_pdf(double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01_PDF evaluates the PDF of a standard normal distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_NORMAL_01_PDF, the value of the PDF at RVAL.
//
        {
            double pi = 3.141592653589793;
            double value;

            value = Math.Exp(-0.5 * rval * rval) / Math.Sqrt(2.0 * pi);

            return value;
        }

        public static double r8_normal_01_sample()

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01_SAMPLE samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    The Box-Muller method is used, which is efficient, but
//    generates two values at a time.
//
//    Typically, we would use one value and save the other for the next call.
//    However, the fact that this function has saved memory makes it difficult
//    to correctly handle cases where we want to re-initialize the code,
//    or to run in parallel.  Therefore, we will instead use the first value
//    and DISCARD the second.
//
//    EFFICIENCY must defer to SIMPLICITY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_NORMAL_01_SAMPLE, a normally distributed random value.
//
        {
            const double pi = 3.14159265358979323;
            double r1;
            double r2;
            double x;

            r1 = r8_uniform_01_sample();
            r2 = r8_uniform_01_sample();

            x = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Cos(2.0 * pi * r2);

            return x;
        }

        public static double r8_scinvchi_pdf(double df, double s, double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_SCINVCHI_PDF: PDF for a scaled inverse chi-squared distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Input, double S, the scale factor.
//    0.0 < S.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_SCINVCHI_PDF, the value of the PDF at RVAL.
//    inverse-chi-square distribution.
//
        {
            double temp1;
            double temp2;
            double value;

            if (df <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_SCINVCHI_PDF - Fatal error!");
                Console.WriteLine("  Degrees of freedom must be positive.");
                return (1);
            }

            if (s <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_SCINVCHI_PDF - Fatal error!");
                Console.WriteLine("  Scale parameter must be positive.");
                return (1);
            }

            if (rval <= 0.0)
            {
                value = 0.0;
            }
            else
            {
                temp2 = df * 0.5;

                temp1 = temp2 * Math.Log(temp2) + temp2 * Math.Log(s)
                        - (temp2 * s / rval)
                        - (temp2 + 1.0) * Math.Log(rval) - typeMethods.r8_gamma_log(temp2);

                value = Math.Exp(temp1);
            }

            return value;
        }

        public static double r8_scinvchi_sample(double df, double s)

//****************************************************************************80
//
//  Purpose:
//
//    R8_SCINVCHI_SAMPLE: sample a scaled inverse chi-squared distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Input, double S, the scale factor.
//    0.0 < S.
//
//    Output, double R8_SCINVCHI_SAMPLE, a sample of the distribution
//
        {
            double a;
            double b;
            double value;

            a = 0.5 * df * s;
            b = 0.5 * df;
            value = r8_gamma_sample(a, b);
            if (value != 0.0)
            {
                value = 1.0 / value;
            }

            return value;
        }

        public static double r8_uniform_pdf(double lower, double upper, double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_PDF evaluates the PDF of a uniform distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double LOWER, UPPER, the lower and upper range limits.
//    LOWER < UPPER.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_UNIFORM_PDF, the value of the PDF at RVAL.
//
        {
            double value;

            if (upper <= lower)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_UNIFORM_PDF - Fatal error!");
                Console.WriteLine("  For uniform PDF, the lower limit must be");
                Console.WriteLine("  less than the upper limit");
                return (1);
            }

            if (rval < lower)
            {
                value = 0.0;
            }
            else if (rval <= upper)
            {
                value = 1.0 / (upper - lower);
            }
            else
            {
                value = 0.0;
            }

            return value;
        }

        public static double r8_uniform_sample(double low, double high)

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_SAMPLE generates a uniform random deviate.
//
//  Discussion:
//
//    This procedure generates a real deviate uniformly distributed between
//    LOW and HIGH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double LOW, HIGH, the lower and upper bounds.
//
//    Output, double R8_UNIFORM_SAMPLE, a random deviate from the distribution.
//
        {
            double value;

            value = low + (high - low) * r8_uniform_01_sample();

            return value;
        }

        public static double r8_uniform_01_pdf(double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01_PDF evaluates the PDF of a standard uniform distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_UNIFORM_01_PDF, the value of the PDF at RVAL.
//
        {
            double value;

            if (rval < 0.0)
            {
                value = 0.0;
            }
            else if (rval <= 1.0)
            {
                value = 1.0;
            }
            else
            {
                value = 0.0;
            }

            return value;
        }

        public static double r8_uniform_01_sample()

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01_SAMPLE generates a uniform random deviate from [0,1].
//
//  Discussion:
//
//    This function should be the only way that the package accesses random
//    numbers.
//
//    Setting OPTION to 0 accesses the R8_UNI_01() function in RNGLIB,
//    for which there are versions in various languages, which should result
//    in the same values being returned.  This should be the only place
//    in this library that accesses a function from RNGLIB.
//
//    Setting OPTION to 1 in the C++ version calls the system
//    random number generator "rand()".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Output, double R8_UNIFORM_01_SAMPLE, a random deviate.
//
        {
            const int option = 1;
            double value;

            value = entropyRNG.RNG.nextdouble();

/*  
if ( option == 0 )
{
value = r8_uni_01 ( );
}
else
{
value = ( double ) rand ( ) / ( double ) RAND_MAX;
}
*/
            return value;
        }

        public static double r8vec_multinormal_pdf(int n, double[] mu, double[] r, double c_det,
            double[] x)

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MULTINORMAL_PDF evaluates a multivariate normal PDF.
//
//  Discussion:
//
//    PDF ( MU(1:N), C(1:N,1:N); X(1:N) ) = 
//      1 / ( 2 * pi ) ^ ( N / 2 ) * 1 / sqrt ( det ( C ) )
//      * exp ( - ( X - MU )' * inverse ( C ) * ( X - MU ) / 2 )
//
//    Here,
//
//      X is the argument vector of length N,
//      MU is the mean vector of length N,
//      C is an N by N positive definite symmetric covariance matrix.
//
//    The properties of C guarantee that it has an upper triangular
//    matrix R, the Cholesky factor, such that C = R' * R.  It is the
//    matrix R that is required by this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double MU[N], the mean vector.
//
//    Input, double R[N*N], the upper triangular Cholesky
//    factor of the covariance matrix C.
//
//    Input, double C_DET, the determinant of the
//    covariance matrix C.
//
//    Input, double X[N], a sample of the distribution.
//
//    Output, double R8VEC_MULTINORMAL_PDF, the PDF evaluated
//    at X.
//
        {
            double[] b;
            int i;
            double pdf;
            double pi = 3.141592653589793;
            double xcx;
            double[] y;
//
//  Compute:
//    inverse(R')*(x-mu) = y
//  by solving:
//    R'*y = x-mu
//
            b = new double[n];

            for (i = 0; i < n; i++)
            {
                b[i] = x[i] - mu[i];
            }

            y = typeMethods.r8mat_utsol(n, r, b);
//
//  Compute:
//    (x-mu)' * inv(C)          * (x-mu)
//  = (x-mu)' * inv(R'*R)       * (x-mu)
//  = (x-mu)' * inv(R) * inv(R) * (x-mu)
//  = y' * y.
//
            xcx = typeMethods.r8vec_dot_product(n, y, y);

            pdf = 1.0 / Math.Sqrt(Math.Pow(2.0 * pi, n))
                  * 1.0 / Math.Sqrt(c_det)
                  * Math.Exp(-0.5 * xcx);

            return pdf;
        }

        public static double[] r8vec_multinormal_sample(int n, double[] mu, double[] r)

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MULTINORMAL_SAMPLE samples a multivariate normal PDF.
//
//  Discussion:
//
//    PDF ( MU(1:N), C(1:N,1:N); X(1:N) ) = 
//      1 / ( 2 * pi ) ^ ( N / 2 ) * 1 / det ( C )
//      * exp ( - ( X - MU )' * inverse ( C ) * ( X - MU ) / 2 )
//
//    Here,
//
//      X is the argument vector of length N,
//      MU is the mean vector of length N,
//      C is an N by N positive definite symmetric covariance matrix.
//
//    The properties of C guarantee that it has an upper triangular
//    matrix R, the Cholesky factor, such that C = R' * R.  It is the
//    matrix R that is required by this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double MU[N], the mean vector.
//
//    Input, double R[N*N], the upper triangular Cholesky
//    factor of the covariance matrix C.
//
//    Output, double R8VEC_MULTINORMAL_SAMPLE[N], a sample of the distribution.
//
        {
            int i;
            int j;
            double[] x;
            double[] z;
//
//  Compute X = MU + R' * Z
//  where Z is a vector of standard normal variates.
//
            z = new double[n];
            for (j = 0; j < n; j++)
            {
                z[j] = r8_normal_01_sample();
            }

            x = new double[n];
            for (i = 0; i < n; i++)
            {
                x[i] = mu[i];
                for (j = 0; j <= i; j++)
                {
                    x[i] = x[i] + r[j + i * n] * z[j];
                }
            }

            return x;
        }


    }
}