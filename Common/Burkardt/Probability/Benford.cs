using System;
using Burkardt.Types;

namespace Burkardt.Probability;

public static class Benford
{
    public static double benford_cdf(int x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BENFORD_CDF returns the Benford CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the string of significant digits to be 
        //    checked.  If X is 1, then we are asking for the Benford probability that
        //    a value will have first digit 1.  If X is 123, we are asking for
        //    the probability that the first three digits will be 123, and so on.
        //
        //    Output, double BENFORD_CDF, the Benford probability that an item taken
        //    from a real world distribution will have the initial digit X or less.
        //
    {
        double value = 0;

        switch (x)
        {
            case <= 0:
                value = 0.0;
                break;
            default:
            {
                if (typeMethods.i4_is_power_of_10(x + 1))
                {
                    value = 1.0;
                }
                else
                {
                    value = Math.Log10(x + 1);
                    value %= 1.0;
                }

                break;
            }
        }

        return value;
    }

    public static double benford_pdf(int x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BENFORD_PDF returns the Benford probability of one or more significant digits.
        //
        //  Discussion:
        //
        //    Benford's law is an empirical formula explaining the observed
        //    distribution of initial digits in lists culled from newspapers,
        //    tax forms, stock market prices, and so on.  It predicts the observed
        //    high frequency of the initial digit 1, for instance.
        //
        //    Note that the probabilities of digits 1 through 9 are guaranteed
        //    to add up to 1, since
        //      LOG10 ( 2/1 ) + LOG10 ( 3/2) + LOG10 ( 4/3 ) + ... + LOG10 ( 10/9 )
        //      = LOG10 ( 2/1 * 3/2 * 4/3 * ... * 10/9 ) = LOG10 ( 10 ) = 1.
        //
        //  Formula:
        //
        //    PDF(X) = LOG10 ( ( X + 1 ) / X ).
        //
        //  Reference:
        //
        //    Frank Benford,
        //    The Law of Anomalous Numbers,
        //    Proceedings of the American Philosophical Society,
        //    Volume 78, pages 551-572, 1938.
        //
        //    Ted Hill,
        //    The First Digit Phenomenon,
        //    American Scientist,
        //    Volume 86, July/August 1998, pages 358 - 363.
        //
        //    R Raimi,
        //    The Peculiar Distribution of First Digits,
        //    Scientific American,
        //    December 1969, pages 109-119.
        //
        //  Modified:
        //
        //    13 August 1998
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the string of significant digits to be checked.
        //    If X is 1, then we are asking for the Benford probability that
        //    a value will have first digit 1.  If X is 123, we are asking for
        //    the probability that the first three digits will be 123, and so on.
        //
        //    Output, double BENFORD_PDF, the Benford probability that an item taken
        //    from a real world distribution will have the initial digits X.
        //
    {
        double pdf = x switch
        {
            <= 0 => 0.0,
            _ => Math.Log10((x + 1) / (double) x)
        };

        return pdf;
    }
}