using System;
using Burkardt.Types;

namespace Burkardt.Probability
{
    public static class Multinomial
    {


        public static bool multinomial_check(int a, int b, double[] c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTINOMIAL_CHECK checks the parameters of the Multinomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the number of trials.
        //
        //    Input, int B, the number of outcomes possible on one trial.
        //    1 <= B.
        //
        //    Input, double C(B).  C(I) is the probability of outcome I on
        //    any trial.
        //    0.0 <= C(I) <= 1.0,
        //    Sum ( 1 <= I <= B ) C(I) = 1.0.
        //
        //    Output, bool MULTINOMIAL_CHECK, is true if the parameters are legal.
        //
        {
            if (b < 1)
            {
                Console.WriteLine(" ");
                Console.WriteLine("MULTINOMIAL_CHECK - Warning!");
                Console.WriteLine("  B < 1.");
                return false;
            }

            for (int i = 0; i < b; i++)
            {
                if (c[i] < 0.0 || 1.0 < c[i])
                {
                    Console.WriteLine(" ");
                    Console.WriteLine("MULTINOMIAL_CHECK - Warning!");
                    Console.WriteLine("  Input C(I) is out of range.");
                    return false;
                }
            }

            double c_sum = typeMethods.r8vec_sum(b, c);

            if (0.0001 < Math.Abs(1.0 - c_sum))
            {
                Console.WriteLine(" ");
                Console.WriteLine("MULTINOMIAL_CHECK - Warning!");
                Console.WriteLine("  The probabilities do not sum to 1.");
                return false;
            }

            return true;
        }

        public static double[] multinomial_covariance(int a, int b, double[] c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTINOMIAL_COVARIANCE returns the covariances of the Multinomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 February 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the number of trials.
        //
        //    Input, int B, the number of outcomes possible on one trial.
        //    1 <= B.
        //
        //    Input, double C(B).  C(I) is the probability of outcome I on
        //    any trial.
        //    0.0 <= C(I) <= 1.0,
        //    SUM ( 1 <= I <= B) C(I) = 1.0.
        //
        //    Output, double MULTINOMIAL_COVARIANCE[B*B], the covariance matrix.
        //
        {
            double[] covariance = new double[b * b];

            for (int i = 0; i < b; i++)
            {
                for (int j = 0; j < b; j++)
                {
                    if (i == j)
                    {
                        covariance[i + j * b] = (double) (a) * c[i] * (1.0 - c[i]);
                    }
                    else
                    {
                        covariance[i + j * b] = -(double) (a) * c[i] * c[j];
                    }
                }
            }

            return covariance;
        }

        public static double[] multinomial_mean(int a, int b, double[] c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTINOMIAL_MEAN returns the means of the Multinomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the number of trials.
        //
        //    Input, int B, the number of outcomes possible on one trial.
        //    1 <= B.
        //
        //    Input, double C(B).  C(I) is the probability of outcome I on
        //    any trial.
        //    0.0 <= C(I) <= 1.0,
        //    SUM ( 1 <= I <= B) C(I) = 1.0.
        //
        //    Output, double MEAN(B), MEAN(I) is the expected value of the
        //    number of outcome I in N trials.
        //
        {
            double[] mean = new double[b];

            for (int i = 0; i < b; i++)
            {
                mean[i] = (double) (a) * c[i];
            }

            return mean;
        }

        public static double multinomial_pdf(int[] x, int a, int b, double[] c )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTINOMIAL_PDF computes a Multinomial PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B,C;X) = Comb(A,B,X) * Product ( 1 <= I <= B ) C(I)^X(I)
        //
        //    where Comb(A,B,X) is the multinomial coefficient
        //      C( A; X(1), X(2), ..., X(B) )
        //
        //    PDF(A,B,C;X) is the probability that in A trials there
        //    will be exactly X(I) occurrences of event I, whose probability
        //    on one trial is C(I), for I from 1 to B.
        //
        //    As soon as A or B gets large, the number of possible X's explodes,
        //    and the probability of any particular X can become extremely small.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X[B]; X(I) counts the number of occurrences of
        //    outcome I, out of the total of A trials.
        //
        //    Input, int A, the total number of trials.
        //
        //    Input, int B, the number of different possible outcomes on
        //    one trial.
        //
        //    Input, double C[B]; C(I) is the probability of outcome I on
        //    any one trial.
        //
        //    Output, double MULTINOMIAL_PDF, the value of the multinomial PDF.
        //
        {
            //
            //  To try to avoid overflow, do the calculation in terms of logarithms.
            //  Note that Gamma(A+1) = A factorial.
            //
            double pdf_log = Helpers.LogGamma((double) (a + 1));

            for (int i = 0; i < b; i++)
            {
                pdf_log = pdf_log + x[i] * Math.Log(c[i])
                          - Helpers.LogGamma((double) (x[i] + 1));
            }

            double pdf = Math.Exp(pdf_log);

            return pdf;
        }

        public static int[] multinomial_sample(int a, int b, double[] c, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTINOMIAL_SAMPLE samples the Multinomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Luc Devroye,
        //    Non-Uniform Random Variate Generation,
        //    Springer-Verlag, New York, 1986, page 559.
        //
        //  Parameters:
        //
        //    Input, int A, the total number of trials.
        //    0 <= A.
        //
        //    Input, int B, the number of outcomes possible on one trial.
        //    1 <= B.
        //
        //    Input, double C[B].  C(I) is the probability of outcome I on
        //    any trial.
        //    0.0 <= C(I) <= 1.0,
        //    SUM ( 1 <= I <= B) C(I) = 1.0.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int MULTINOMIAL_SAMPLE[B]; The Ith entry is the number of
        //    occurrences of event I during the N trials.
        //
        {
            int[] x = new int[b];
            int ntot = a;

            double sum2 = 1.0;

            for (int i = 0; i < b; i++)
            {
                x[i] = 0;
            }

            for (int ifactor = 0; ifactor < b - 1; ifactor++)
            {
                double prob = c[ifactor] / sum2;
                //
                //  Generate a binomial random deviate for NTOT trials with
                //  single trial success probability PROB.
                //
                x[ifactor] = Binomial.binomial_sample(ntot, prob, ref seed);

                ntot = ntot - x[ifactor];
                if (ntot <= 0)
                {
                    return x;
                }

                sum2 = sum2 - c[ifactor];
            }

            //
            //  The last factor gets what's left.
            //
            x[b - 1] = ntot;

            return x;
        }

        public static double[] multinomial_variance(int a, int b, double[] c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTINOMIAL_VARIANCE returns the variances of the Multinomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the number of trials.
        //
        //    Input, int B, the number of outcomes possible on one trial.
        //    1 <= B.
        //
        //    Input, double C[B].  C(I) is the probability of outcome I on
        //    any trial.
        //    0.0 <= C(I) <= 1.0,
        //    SUM ( 1 <= I <= B ) C(I) = 1.0.
        //
        //    Output, double VARIANCE(B), VARIANCE(I) is the variance of the
        //    total number of events of type I.
        //
        {
            double[] variance = new double[b];

            for (int i = 0; i < b; i++)
            {
                variance[i] = (double) (a) * c[i] * (1.0 - c[i]);
            }

            return variance;
        }
    }
}