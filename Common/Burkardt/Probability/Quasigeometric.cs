using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Quasigeometric
{
    public static double quasigeometric_cdf(int x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUASIGEOMETRIC_CDF evaluates the Quasigeometric CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the maximum number of trials.
        //
        //    Input, double A, the probability of 0 successes.
        //    0.0 <= A <= 1.0.
        //
        //    Input, double B, the depreciation constant.
        //    0.0 <= B < 1.0.
        //
        //    Output, double QUASIGEOMETRIC_CDF, the value of the CDF.
        //
    {
        double cdf;

        switch (x)
        {
            case < 0:
                cdf = 0.0;
                break;
            case 0:
                cdf = a;
                break;
            default:
            {
                cdf = b switch
                {
                    0.0 => 1.0,
                    _ => a + (1.0 - a) * (1.0 - Math.Pow(b, x))
                };

                break;
            }
        }

        return cdf;
    }

    public static int quasigeometric_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUASIGEOMETRIC_CDF_INV inverts the Quasigeometric CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0
        //
        //    Input, double A, the probability of 0 successes.
        //    0.0 <= A <= 1.0.
        //
        //    Input, double B, the depreciation constant.
        //    0.0 <= B < 1.0.
        //
        //    Output, int QUASIGEOMETRIC_CDF_INV, the corresponding value of X.
        //
    {
        int x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("QUASIGEOMETRIC_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        if (cdf < a)
        {
            x = 0;
        }
        else
        {
            x = b switch
            {
                0.0 => 1,
                _ => 1 + (int) ((Math.Log(1.0 - cdf) - Math.Log(1.0 - a)) / Math.Log(b))
            };
        }

        return x;
    }

    public static bool quasigeometric_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUASIGEOMETRIC_CHECK checks the parameters of the Quasigeometric CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the probability of 0 successes.
        //    0.0 <= A <= 1.0.
        //
        //    Input, double B, the depreciation constant.
        //    0.0 <= B < 1.0.
        //
        //    Output, bool QUASIGEOMETRIC_CHECK, is true if the parameters are legal.
        //
    {
        bool check = true;

        switch (a)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("QUASIGEOMETRIC_CHECK - Warning!");
                Console.WriteLine("  A < 0 or 1 < A.");
                check = false;
                break;
        }

        switch (b)
        {
            case < 0.0:
            case >= 1.0:
                Console.WriteLine("");
                Console.WriteLine("QUASIGEOMETRIC_CHECK - Warning!");
                Console.WriteLine("  B < 0 or 1 <= B.");
                check = false;
                break;
        }

        return check;
    }


    public static double quasigeometric_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUASIGEOMETRIC_MEAN returns the mean of the Quasigeometric PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the probability of 0 successes.
        //    0.0 <= A <= 1.0.
        //
        //    Input, double B, the depreciation constant.
        //    0.0 <= B < 1.0.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
    {
        double mean = (1.0 - a) / (1.0 - b);

        return mean;
    }

    public static double quasigeometric_pdf(int x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUASIGEOMETRIC_PDF evaluates the Quasigeometric PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) =    A                     if 0  = X;
        //               = (1-A) * (1-B) * B^(X-1)  if 1 <= X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Darren Glass, Philip Lowry,
        //    Quasiquasigeometric Distributions and Extra Inning Baseball Games,
        //    Mathematics Magazine,
        //    Volume 81, Number 2, April 2008, pages 127-137.
        //
        //    Paul Nahin,
        //    Digital Dice: Computational Solutions to Practical Probability Problems,
        //    Princeton University Press, 2008,
        //    ISBN13: 978-0-691-12698-2,
        //    LC: QA273.25.N34.
        //
        //  Parameters:
        //
        //    Input, int X, the independent variable.
        //    0 <= X
        //
        //    Input, double A, the probability of 0 successes.
        //    0.0 <= A <= 1.0.
        //
        //    Input, double B, the depreciation constant.
        //    0.0 <= B < 1.0.
        //
        //    Output, double PDF, the value of the PDF.
        //
    {
        double pdf;

        switch (x)
        {
            case < 0:
                pdf = 0.0;
                break;
            case 0:
                pdf = a;
                break;
            default:
            {
                pdf = b switch
                {
                    0.0 when x == 1 => 1.0,
                    0.0 => 0.0,
                    _ => (1.0 - a) * (1.0 - b) * Math.Pow(b, x - 1)
                };

                break;
            }
        }

        return pdf;
    }

    public static int quasigeometric_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUASIGEOMETRIC_SAMPLE samples the Quasigeometric PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the probability of 0 successes.
        //    0.0 <= A <= 1.0.
        //
        //    Input, double B, the depreciation constant.
        //    0.0 <= B < 1.0.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, int QUASIGEOMETRIC_SAMPLE, a sample of the PDF.
        //
    {
        double cdf;
        int x;

        cdf = UniformRNG.r8_uniform_01(ref seed);

        x = quasigeometric_cdf_inv(cdf, a, b);

        return x;
    }

    public static double quasigeometric_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUASIGEOMETRIC_VARIANCE returns the variance of the Quasigeometric PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the probability of 0 successes.
        //    0.0 <= A <= 1.0.
        //
        //    Input, double B, the depreciation constant.
        //    0.0 <= B < 1.0.
        //
        //    Output, double QUASIGEOMETRIC_VARIANCE, the variance of the PDF.
        //
    {
        double variance = (1.0 - a) * (a + b) / (1.0 - b) / (1.0 - b);

        return variance;
    }
}