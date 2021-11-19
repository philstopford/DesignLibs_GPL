using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Birthday
{
    public static double birthday_cdf(int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIRTHDAY_CDF returns the Birthday Concurrence CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of people whose birthdays have been
        //    disclosed.
        //
        //    Output, double BIRTHDAY_CDF, the probability that at least
        //    two of the N people have matching birthays..
        //
    {
        double cdf;

        switch (n)
        {
            case < 1:
                cdf = 0.0;
                return cdf;
            case > 365:
                cdf = 1.0;
                return cdf;
        }

        //
        //  Compute the probability that N people have distinct birthdays.
        //
        cdf = 1.0;
        for (int i = 1; i <= n; i++)
        {
            cdf = cdf * (365 + 1 - i) / 365.0;
        }

        //
        //  Compute the probability that it is NOT the case that N people
        //  have distinct birthdays.  This is the cumulative probability
        //  that person 2 matches person 1, or person 3 matches 1 or 2,
        //  etc.
        //
        cdf = 1.0 - cdf;

        return cdf;
    }

    public static int birthday_cdf_inv(double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIRTHDAY_CDF_INV inverts the Birthday Concurrence CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the probability that at least
        //    two of the N people have matching birthays..
        //
        //    Output, int BIRTHDAY_CDF_INV, the corresponding number of people whose
        //    birthdays need to be disclosed.
        //
    {
        int n;

        switch (cdf)
        {
            case <= 0.0:
                n = 0;
                return n;
            case >= 1.0:
                n = 365;
                return n;
        }

        //
        //  Compute the probability that N people have distinct birthdays.
        //
        double cdf_not = 1.0;

        for (int i = 1; i <= 365; i++)
        {
            cdf_not = cdf_not * (365 + 1 - i) / 365.0;
            if (!(cdf <= 1.0 - cdf_not))
            {
                continue;
            }

            n = i;
            return n;
        }

        n = 365;
        return n;
    }

    public static double birthday_pdf(int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIRTHDAY_PDF returns the Birthday Concurrence PDF.
        //
        //  Discussion:
        //
        //    The probability is the probability that the N-th person is the
        //    first one to match a birthday with someone earlier.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of people whose birthdays have been
        //    disclosed.
        //
        //    Output, double BIRTHDAY_PDF, the probability that the N-th person
        //    is the first to match a birthday with someone earlier.
        //
    {
        double pdf;

        switch (n)
        {
            case < 1:
            case > 365:
                pdf = 0.0;
                return pdf;
        }

        pdf = 1.0;
        //
        //  Compute the probability that N-1 people have distinct birthdays.
        //
        for (int i = 1; i <= n - 1; i++)
        {
            pdf = pdf * (365 + 1 - i) / 365.0;
        }

        //
        //  Compute the probability that person N has one of those N-1 birthdays.
        //
        pdf = pdf * (n - 1) / 365.0;

        return pdf;
    }

    public static int birthday_sample(int n, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIRTHDAY_SAMPLE samples the Birthday Concurrence PDF.
        //
        //  Discussion:
        //
        //    The probability is the probability that the N-th person is the
        //    first one to match a birthday with someone earlier.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of people whose birthdays have been
        //    disclosed.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int BIRTHDAY_SAMPLE,
        //    * 1 if the first N-1 people had distinct
        //      birthdays, but person N had a birthday in common with a previous person,
        //    * 0 otherwise.
        //
    {
        int value;

        switch (n)
        {
            case < 1:
                value = 0;
                return value;
        }

        //
        //  Choose N birthdays at random.
        //
        int[] b = UniformRNG.i4vec_uniform_ab_new(n, 1, 365, ref seed);
        //
        //  Are the first N-1 birthdays unique?
        //
        int u1 = typeMethods.i4vec_unique_count(n - 1, b);

        if (u1 < n - 1)
        {
            value = 0;
            return value;
        }

        //
        //  Does the N-th birthday match an earlier one?
        //
        int u2 = typeMethods.i4vec_unique_count(n, b);

        value = u2 == n - 1 ? 1 : 0;

        return value;
    }

}