using System;
using Burkardt.Types;

namespace Burkardt.PDFLib;

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
        double value;

        switch (alpha)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_BETA_PDF - Fatal error!");
                Console.WriteLine("  Parameter ALPHA is not positive.");
                return 1;
        }

        switch (beta)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_BETA_PDF- Fatal error!");
                Console.WriteLine("  Parameter BETA is not positive.");
                return 1;
        }

        switch (rval)
        {
            case <= 0.0:
            case >= 1.0:
                value = 0.0;
                break;
            default:
                double temp = typeMethods.r8_gamma_log(alpha + beta) - typeMethods.r8_gamma_log(alpha)
                                                                     - typeMethods.r8_gamma_log(beta);

                value = Math.Exp(temp) * Math.Pow(rval, alpha - 1.0)
                                       * Math.Pow(1.0 - rval, beta - 1.0);
                break;
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
        const double log4 = 1.3862943611198906188;
        const double log5 = 1.6094379124341003746;
        double u1;
        double u2;
        double v;
        double value;
        double w;
        double z;

        switch (aa)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_BETA_SAMPLE - Fatal error!");
                Console.WriteLine("  AA <= 0.0");
                return 1;
        }

        switch (bb)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_BETA_SAMPLE - Fatal error!");
                Console.WriteLine("  BB <= 0.0");
                return 1;
        }

        switch (aa)
        {
            //
            //  Algorithm BB
            //
            case > 1.0 when 1.0 < bb:
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
                double gamma = a + 1.0 / beta;

                for (;;)
                {
                    u1 = r8_uniform_01_sample();
                    u2 = r8_uniform_01_sample();
                    v = beta * Math.Log(u1 / (1.0 - u1));
                    w = a * Math.Exp(v);

                    z = u1 * u1 * u2;
                    double r = gamma * v - log4;
                    double s = a + r - w;

                    if (5.0 * z <= s + 1.0 + log5)
                    {
                        break;
                    }

                    double t = Math.Log(z);
                    if (t <= s)
                    {
                        break;
                    }

                    if (t <= r + alpha * Math.Log(alpha / (b + w)))
                    {
                        break;
                    }
                }

                break;
            }
            //
            default:
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
                double delta = 1.0 + a - b;
                double k1 = delta * (1.0 / 72.0 + b / 24.0)
                            / (a / b - 7.0 / 9.0);
                double k2 = 0.25 + (0.5 + 0.25 / delta) * b;

                for (;;)
                {
                    u1 = r8_uniform_01_sample();
                    u2 = r8_uniform_01_sample();

                    switch (u1)
                    {
                        case < 0.5:
                        {
                            double y = u1 * u2;
                            z = u1 * y;

                            if (k1 <= 0.25 * u2 + z - y)
                            {
                                continue;
                            }

                            break;
                        }
                        default:
                        {
                            z = u1 * u1 * u2;

                            switch (z)
                            {
                                case <= 0.25:
                                {
                                    v = beta * Math.Log(u1 / (1.0 - u1));
                                    w = a * Math.Exp(v);

                                    if (Math.Abs(aa - a) <= typeMethods.r8_epsilon())
                                    {
                                        value = w / (b + w);
                                    }
                                    else
                                    {
                                        value = b / (b + w);
                                    }

                                    return value;
                                }
                            }

                            if (k2 < z)
                            {
                                continue;
                            }

                            break;
                        }
                    }

                    v = beta * Math.Log(u1 / (1.0 - u1));
                    w = a * Math.Exp(v);

                    if (Math.Log(z) <= alpha * (Math.Log(alpha / (b + w)) + v) - log4)
                    {
                        break;
                    }
                }

                break;
            }
        }

        if (Math.Abs(aa - a) <= typeMethods.r8_epsilon())
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
        double value;

        switch (df)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_CHI_PDF - Fatal error!");
                Console.WriteLine("  Degrees of freedom must be positive.");
                return 1;
        }

        switch (rval)
        {
            case <= 0.0:
                value = 0.0;
                break;
            default:
                double temp2 = df * 0.5;

                double temp1 = (temp2 - 1.0) * Math.Log(rval) - 0.5 * rval
                                                              - temp2 * Math.Log(2.0) - typeMethods.r8_gamma_log(temp2);

                value = Math.Exp(temp1);
                break;
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
        switch (df)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_CHI_SAMPLE - Fatal error!");
                Console.WriteLine("  DF <= 0.");
                Console.WriteLine("  Value of DF: " + df + "");
                return 1;
        }

        const double arg1 = 1.0;
        double arg2 = df / 2.0;

        double value = 2.0 * r8_gamma_sample(arg1, arg2);

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
        int mn;
        int mx;
        double value = 0;

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

        switch (mn)
        {
            case < 0:
                value = 0.0;
                break;
            case 0:
                value = 1.0;
                break;
            default:
            {
                value = mx + 1;

                int i;
                for (i = 2; i <= mn; i++)
                {
                    value = value * (mx + i) / i;
                }

                break;
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
        double value = 0;

        switch (beta)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_EXPONENTIAL_PDF - Fatal error!");
                Console.WriteLine("  BETA parameter must be positive.");
                return 1;
        }

        value = rval switch
        {
            < 0.0 => 0.0,
            _ => Math.Exp(-rval / beta) / beta
        };

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
        double r = r8_uniform_01_sample();

        double value = -Math.Log(r) * lambda;

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
        double value = rval switch
        {
            < 0.0 => 0.0,
            _ => Math.Exp(-rval)
        };

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
        double r = r8_uniform_01_sample();

        double value = -Math.Log(r);

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
        double value;

        switch (alpha)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_GAMMA_PDF - Fatal error!");
                Console.WriteLine("  Parameter ALPHA is not positive.");
                return 1;
        }

        switch (beta)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_GAMMA_PDF - Fatal error!");
                Console.WriteLine("  Parameter BETA is not positive.");
                return 1;
        }

        switch (rval)
        {
            case <= 0.0:
                value = 0.0;
                break;
            default:
                double temp = alpha * Math.Log(beta) + (alpha - 1.0) * Math.Log(rval)
                              - beta * rval - typeMethods.r8_gamma_log(alpha);

                value = Math.Exp(temp);
                break;
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
        double value = 0;

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
        double value;

        switch (alpha)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_GAMMA_01_PDF - Fatal error!");
                Console.WriteLine("  Parameter ALPHA is not positive.");
                return 1;
        }

        switch (rval)
        {
            case <= 0.0:
                value = 0.0;
                break;
            default:
                double temp = (alpha - 1.0) * Math.Log(rval) - rval - typeMethods.r8_gamma_log(alpha);

                value = Math.Exp(temp);
                break;
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
        const double a1 = 0.3333333;
        const double a2 = -0.2500030;
        const double a3 = 0.2000062;
        const double a4 = -0.1662921;
        const double a5 = 0.1423657;
        const double a6 = -0.1367177;
        const double a7 = 0.1233795;
        double b;
        const double e1 = 1.0;
        const double e2 = 0.4999897;
        const double e3 = 0.1668290;
        const double e4 = 0.0407753;
        const double e5 = 0.0102930;
        const double q1 = 0.04166669;
        const double q2 = 0.02083148;
        const double q3 = 0.00801191;
        const double q4 = 0.00144121;
        const double q5 = -0.00007388;
        const double q6 = 0.00024511;
        const double q7 = 0.00024240;
        const double sqrt32 = 5.6568542494923801952;
        double value = 0;

        switch (a)
        {
            case >= 1.0:
            {
                double s2 = a - 0.5;
                double s = Math.Sqrt(s2);
                double d = sqrt32 - 12.0 * s;
                //
                //  Immediate acceptance.
                //
                double t = r8_normal_01_sample();
                double x = s + 0.5 * t;
                value = x * x;

                switch (t)
                {
                    case >= 0.0:
                        return value;
                }

                //
                //  Squeeze acceptance.
                //
                double u = r8_uniform_01_sample();
                if (d * u <= t * t * t)
                {
                    return value;
                }

                double r = 1.0 / a;
                double q0 = ((((((
                                     q7 * r
                                     + q6) * r
                                 + q5) * r
                                + q4) * r
                               + q3) * r
                              + q2) * r
                             + q1) * r;
                double si;
                double c;
                switch (a)
                {
                    //
                    //  Approximation depending on size of parameter A.
                    //
                    case > 13.022:
                        b = 1.77;
                        si = 0.75;
                        c = 0.1515 / s;
                        break;
                    case > 3.686:
                        b = 1.654 + 0.0076 * s2;
                        si = 1.68 / s + 0.275;
                        c = 0.062 / s + 0.024;
                        break;
                    default:
                        b = 0.463 + s + 0.178 * s2;
                        si = 1.235;
                        c = 0.195 / s - 0.079 + 0.16 * s;
                        break;
                }

                double q;
                double v;
                switch (x)
                {
                    //
                    //  Quotient test.
                    //
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
                            return value;
                        }

                        break;
                    }
                }

                for (;;)
                {
                    double e = r8_exponential_01_sample();
                    u = 2.0 * r8_uniform_01_sample() - 1.0;

                    t = u switch
                    {
                        >= 0.0 => b + Math.Abs(si * e),
                        _ => b - Math.Abs(si * e)
                    };

                    switch (t)
                    {
                        //
                        //  Possible rejection.
                        //
                        case < -0.7187449:
                            continue;
                    }

                    //
                    //  Calculate V and quotient Q.
                    //
                    v = 0.5 * t / s;

                    q = Math.Abs(v) switch
                    {
                        > 0.25 => q0 - s * t + 0.25 * t * t + 2.0 * s2 * Math.Log(1.0 + v),
                        _ => q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v
                    };

                    double w;
                    switch (q)
                    {
                        //
                        //  Hat acceptance.
                        //
                        case <= 0.0:
                            continue;
                        case > 0.5:
                            w = Math.Exp(q) - 1.0;
                            break;
                        default:
                            w = ((((
                                       e5 * q
                                       + e4) * q
                                   + e3) * q
                                  + e2) * q
                                 + e1) * q;
                            break;
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
                break;
            }
            //
            //  Method for A < 1.
            //
            case < 1.0:
            {
                b = 1.0 + 0.3678794 * a;

                for (;;)
                {
                    double p = b * r8_uniform_01_sample();

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

                break;
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
        double value;

        switch (df)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_INVCHI_PDF - Fatal error!");
                Console.WriteLine("  Degrees of freedom must be positive.");
                return 1;
        }

        switch (rval)
        {
            case <= 0.0:
                value = 0.0;
                break;
            default:
                double temp2 = df * 0.5;

                double temp1 = -temp2 * Math.Log(2.0) - (temp2 + 1.0) * Math.Log(rval)
                                                      - 0.5 / rval - typeMethods.r8_gamma_log(temp2);

                value = Math.Exp(temp1);
                break;
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
        const double a = 0.5;
        double b = 0.5 * df;
        double value = r8_gamma_sample(a, b);

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
        double value;

        switch (alpha)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_INVGAM_PDF - Fatal error!");
                Console.WriteLine("  Parameter ALPHA is not positive.");
                return 1;
        }

        switch (beta)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_INVGAM_PDF - Fatal error!");
                Console.WriteLine("  Parameter BETA is not positive.");
                return 1;
        }

        switch (rval)
        {
            case <= 0.0:
                value = 0.0;
                break;
            default:
                double temp = alpha * Math.Log(beta) - (alpha + 1.0) * Math.Log(rval)
                                                     - beta / rval - typeMethods.r8_gamma_log(alpha);

                value = Math.Exp(temp);
                break;
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
        double value = 0;

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
        switch (sd)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_NORMAL_PDF - Fatal error!");
                Console.WriteLine("  Standard deviation must be positive.");
                return 1;
        }

        double rtemp = (rval - av) * (rval - av) * 0.5 / (sd * sd);

        double value = Math.Exp(-rtemp) / sd / Math.Sqrt(2.0 * Math.PI);

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
        double value = 0;

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
            
        double value = 0;

        value = Math.Exp(-0.5 * rval * rval) / Math.Sqrt(2.0 * Math.PI);

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
        double r1 = r8_uniform_01_sample();
        double r2 = r8_uniform_01_sample();

        double x = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Cos(2.0 * Math.PI * r2);

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
        double value;

        switch (df)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_SCINVCHI_PDF - Fatal error!");
                Console.WriteLine("  Degrees of freedom must be positive.");
                return 1;
        }

        switch (s)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_SCINVCHI_PDF - Fatal error!");
                Console.WriteLine("  Scale parameter must be positive.");
                return 1;
        }

        switch (rval)
        {
            case <= 0.0:
                value = 0.0;
                break;
            default:
                double temp2 = df * 0.5;

                double temp1 = temp2 * Math.Log(temp2) + temp2 * Math.Log(s)
                               - temp2 * s / rval
                               - (temp2 + 1.0) * Math.Log(rval) - typeMethods.r8_gamma_log(temp2);

                value = Math.Exp(temp1);
                break;
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
        double a = 0.5 * df * s;
        double b = 0.5 * df;
        double value = r8_gamma_sample(a, b);
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
        double value = 0;

        if (upper <= lower)
        {
            Console.WriteLine("");
            Console.WriteLine("R8_UNIFORM_PDF - Fatal error!");
            Console.WriteLine("  For uniform PDF, the lower limit must be");
            Console.WriteLine("  less than the upper limit");
            return 1;
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
        double value = 0;

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
        double value = rval switch
        {
            < 0.0 => 0.0,
            <= 1.0 => 1.0,
            _ => 0.0
        };

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

        double value = entropyRNG.RNG.nextdouble();

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
            double[] x, int xIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MULTINORMAL_PDF evaluates a multivariate normal PDF.
        //
        //  Discussion:
        //
        //    PDF ( MU(1:N), C(1:N,1:N); X(1:N) ) = 
        //      1 / ( 2 * Math.PI ) ^ ( N / 2 ) * 1 / sqrt ( det ( C ) )
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
        int i;

        //
        //  Compute:
        //    inverse(R')*(x-mu) = y
        //  by solving:
        //    R'*y = x-mu
        //
        double[] b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = x[xIndex + i] - mu[i];
        }

        double[] y = typeMethods.r8mat_utsol(n, r, b);
        //
        //  Compute:
        //    (x-mu)' * inv(C)          * (x-mu)
        //  = (x-mu)' * inv(R'*R)       * (x-mu)
        //  = (x-mu)' * inv(R) * inv(R) * (x-mu)
        //  = y' * y.
        //
        double xcx = typeMethods.r8vec_dot_product(n, y, y);

        double pdf = 1.0 / Math.Sqrt(Math.Pow(2.0 * Math.PI, n))
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
        //      1 / ( 2 * Math.PI ) ^ ( N / 2 ) * 1 / det ( C )
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
        //
        //  Compute X = MU + R' * Z
        //  where Z is a vector of standard normal variates.
        //
        double[] z = new double[n];
        for (j = 0; j < n; j++)
        {
            z[j] = r8_normal_01_sample();
        }

        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = mu[i];
            for (j = 0; j <= i; j++)
            {
                x[i] += r[j + i * n] * z[j];
            }
        }

        return x;
    }


}