using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Binomial
    {
        static double binomial_cdf(double x, int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOMIAL_CDF evaluates the Binomial CDF.
        //
        //  Discussion:
        //
        //    CDF(X)(A,B) is the probability of at most X successes in A trials,
        //    given that the probability of success on a single trial is B.
        //
        //    A sequence of trials with fixed probability of success on
        //    any trial is known as a sequence of Bernoulli trials.
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
        //    Input, int X, the desired number of successes.
        //    0 <= X <= A.
        //
        //    Input, int A, the number of trials.
        //    1 <= A.
        //
        //    Input, double B, the probability of success on one trial.
        //    0.0 <= B <= 1.0.
        //
        //    Output, double CDF, the value of the CDF.
        //
        {
            int cnk;
            double cdf;
            double pr;

            if (x < 0)
            {
                cdf = 0.0;
            }
            else if (a <= x)
            {
                cdf = 1.0;
            }
            else if (b == 0.0)
            {
                cdf = 1.0;
            }
            else if (b == 1.0)
            {
                cdf = 0.0;
            }
            else
            {
                cdf = 0.0;

                for (int j = 0; j <= x; j++)
                {
                    cnk = typeMethods.i4_choose(a, j);

                    pr = (double) (cnk) * Math.Pow(b, j) * Math.Pow((1.0 - b), (a - j));

                    cdf = cdf + pr;

                }

            }

            return cdf;
        }

        static int binomial_cdf_inv(double cdf, int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOMIAL_CDF_INV inverts the Binomial CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, int A, the number of trials.
        //    1 <= A.
        //
        //    Input, double B, the probability of success on one trial.
        //    0.0 <= B <= 1.0.
        //
        //    Output, int BINOMIAL_CDF_INV, the corresponding argument.
        //
        {
            int x2;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("BINOMIAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            }

            int x = 0;
            double cdf2 = 0.0;

            for (x2 = 0; x2 <= a; x2++)
            {
                double pdf = binomial_pdf(x2, a, b);

                cdf2 = cdf2 + pdf;

                if (cdf <= cdf2)
                {
                    x = x2;
                    return x;
                }

            }

            return x;
        }

        static bool binomial_check(int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOMIAL_CHECK checks the parameter of the Binomial PDF.
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
        //    Input, int A, the number of trials.
        //    1 <= A.
        //
        //    Input, double B, the probability of success on one trial.
        //    0.0 <= B <= 1.0.
        //
        //    Output, bool BINOMIAL_CHECK, is true if the parameters are legal.
        //
        {
            if (a < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("BINOMIAL_CHECK - Warning!");
                Console.WriteLine("  A < 1.");
                return false;
            }

            if (b < 0.0 || 1.0 < b)
            {
                Console.WriteLine("");
                Console.WriteLine("BINOMIAL_CHECK - Warning!");
                Console.WriteLine("  B < 0 or 1 < B.");
                return false;
            }

            return true;
        }

        static double binomial_mean(int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOMIAL_MEAN returns the mean of the Binomial PDF.
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
        //    Input, int A, the number of trials.
        //    1 <= A.
        //
        //    Input, double B, the probability of success on one trial.
        //    0.0 <= B <= 1.0.
        //
        //    Output, double BINOMIAL_MEAN, the expected value of the number of
        //    successes in A trials.
        //
        {
            double mean = (double) (a) * b;

            return mean;
        }

        static double binomial_pdf(int x, int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOMIAL_PDF evaluates the Binomial PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) is the probability of exactly X successes in A trials,
        //    given that the probability of success on a single trial is B.
        //
        //    PDF(A,B;X) = C(N,X) * B^X * ( 1.0 - B )^( A - X )
        //
        //    Binomial_PDF(1,B;X) = Bernoulli_PDF(B;X).
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
        //    Input, int X, the desired number of successes.
        //    0 <= X <= A.
        //
        //    Input, int A, the number of trials.
        //    1 <= A.
        //
        //    Input, double B, the probability of success on one trial.
        //    0.0 <= B <= 1.0.
        //
        //    Output, double BINOMIAL_PDF, the value of the PDF.
        //
        {
            int cnk;
            double pdf;

            if (a < 1)
            {
                pdf = 0.0;
            }
            else if (x < 0 || a < x)
            {
                pdf = 0.0;
            }
            else if (b == 0.0)
            {
                if (x == 0)
                {
                    pdf = 1.0;
                }
                else
                {
                    pdf = 0.0;
                }
            }
            else if (b == 1.0)
            {
                if (x == a)
                {
                    pdf = 1.0;
                }
                else
                {
                    pdf = 0.0;
                }
            }
            else
            {
                cnk = typeMethods.i4_choose(a, x);

                pdf = (double) (cnk) * Math.Pow(b, x) * Math.Pow((1.0 - b), (a - x));
            }

            return pdf;
        }

        static int binomial_sample(int a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOMIAL_SAMPLE samples the Binomial PDF.
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
        //    Algorithm BU,
        //    Statistical Computing,
        //    Dekker, 1980.
        //
        //  Parameters:
        //
        //    Input, int A, the number of trials.
        //    1 <= A.
        //
        //    Input, double B, the probability of success on one trial.
        //    0.0 <= B <= 1.0.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int BINOMIAL_SAMPLE, a sample of the PDF.
        //
        {
            int i;
            double u;
            int x;

            x = 0;

            for (i = 1; i <= a; i++)
            {
                u = UniformRNG.r8_uniform_01(ref seed);

                if (u <= b)
                {
                    x = x + 1;
                }

            }

            return x;
        }

        static double binomial_variance(int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOMIAL_VARIANCE returns the variance of the Binomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the number of trials.
        //    1 <= A.
        //
        //    Input, double B, the probability of success on one trial.
        //    0.0 <= B <= 1.0.
        //
        //    Output, double BINOMIAL_VARIANCE, the variance of the PDF.
        //
        {
            double variance = (double) (a) * b * (1.0 - b);

            return variance;
        }
    }
}