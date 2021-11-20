using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Bernoulli
{
    public static double bernoulli_cdf(int x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_CDF evaluates the Bernoulli CDF.
        //
        //
        //  Modified:
        //
        //    22 May 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the number of successes on a single trial.
        //    X = 0 or 1.
        //
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Output, double BERNOULLI_CDF, the value of the CDF.
        //
    {
        double cdf = x switch
        {
            < 0 => 0.0,
            0 => 1.0 - a,
            _ => 1.0
        };

        return cdf;
    }

    public static int bernoulli_cdf_inv(double cdf, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_CDF_INV inverts the Bernoulli CDF.
        //
        //  Modified:
        //
        //    22 May 2004
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
        //    Input, double A, the parameter of the PDF.
        //    0.0 <= A <= 1.0.
        //
        //    Output, int BERNOULLI_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("BERNOULLI_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        int x = cdf <= 1.0 - a ? 0 : 1;

        return x;
    }

    public static bool bernoulli_check(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_CHECK checks the parameter of the Bernoulli CDF.
        //
        //  Modified:
        //
        //    22 May 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 <= A <= 1.0.
        //
        //    Output, bool BERNOULLI_CHECK, is TRUE if the data is acceptable.
    {
        switch (a)
        {
            case < 0.0:
            case > 1.0:
                return false;
            default:
                return true;
        }
    }

    public static double bernoulli_mean(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_MEAN returns the mean of the Bernoulli PDF.
        //
        //  Modified:
        //
        //    22 May 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the probability of success.
        //    0.0 <= A <= 1.0.
        //
        //    Output, double BERNOULLI_MEAN, the mean of the PDF.
        //
    {
        return a;
    }

    public static double bernoulli_pdf(int x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_PDF evaluates the Bernoulli PDF.
        //
        //  Formula:
        //
        //    PDF(X)(A) = A^X * ( 1.0 - A )^( X - 1 )
        //
        //    X = 0 or 1.
        //
        //  Discussion:
        //
        //    The Bernoulli PDF describes the simple case in which a single trial
        //    is carried out, with two possible outcomes, called "success" and
        //    "failure"; the probability of success is A.
        //
        //  Modified:
        //
        //    22 May 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the number of successes on a single trial.
        //    X = 0 or 1.
        //
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Output, double BERNOULLI_PDF, the value of the PDF.
        //
    {
        double pdf = x switch
        {
            < 0 => 0.0,
            0 => 1.0 - a,
            1 => a,
            _ => 0.0
        };

        return pdf;
    }

    public static int bernoulli_sample(double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_SAMPLE samples the Bernoulli PDF.
        //
        //  Modified:
        //
        //    22 May 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int BERNOULLI_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        int x = bernoulli_cdf_inv(cdf, a);

        return x;
    }

    public static double bernoulli_variance(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNOULLI_VARIANCE returns the variance of the Bernoulli PDF.
        //
        //  Modified:
        //
        //    22 May 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Output, double BERNOULLI_VARIANCE, the variance of the PDF.
        //
    {
        double variance = a * (1.0 - a);

        return variance;
    }
}