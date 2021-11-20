using System;
using Burkardt.Types;

namespace Burkardt.PolynomialNS;

public static class Multinomial
{
    public static int commul(int n, int nfact, int[] iarray)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMMUL computes a multinomial combinatorial coefficient.
        //
        //  Discussion:
        //
        //    The multinomial coefficient is a generalization of the binomial
        //    coefficient.  It may be interpreted as the number of combinations of
        //    N objects, where IARRAY(1) objects are indistinguishable of type 1,
        //    ... and IARRAY(K) are indistinguishable of type NFACT.
        //
        //    The formula is:
        //
        //      COMMUL = N! / ( IARRAY(1)! IARRAY(2)! ... IARRAY(NFACT)! )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, determines the numerator.
        //
        //    Input, int NFACT, the number of factors in the numerator.
        //
        //    Input, int IARRAY(NFACT).
        //    IARRAY contains the NFACT values used in the denominator.
        //    Note that the sum of these entries should be N,
        //    and that all entries should be nonnegative.
        //
        //    Output, int COMMUL, the value of the multinomial coefficient.
        //
    {
        int i;

        for (i = 0; i < nfact; i++)
        {
            switch (iarray[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("COMMUL - Fatal error");
                    Console.WriteLine("  Entry " + i + " of IARRAY = " + iarray[i] + "");
                    Console.WriteLine("  But this value must be nonnegative.");
                    return 1;
            }
        }

        int isum = 0;
        for (i = 0; i < nfact; i++)
        {
            isum += iarray[i];
        }

        if (isum != n)
        {
            Console.WriteLine("");
            Console.WriteLine("COMMUL - Fatal error!");
            Console.WriteLine("  The sum of the IARRAY entries is " + isum + "");
            Console.WriteLine("  But it must equal N = " + n + "");
            return 1;
        }

        double arg = n + 1;
        double facn = Helpers.LogGamma(arg);

        for (i = 0; i < nfact; i++)
        {
            arg = iarray[i] + 1;
            double fack = Helpers.LogGamma(arg);
            facn -= fack;
        }

        return (int)typeMethods.r8_nint(Math.Exp(facn));
    }

    public static bool multicoef_check(int nfactor, int[] factor)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTICOEF_CHECK checks the parameters of the multinomial coefficient.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NFACTOR, the number of factors.
        //    1 <= NFACTOR.
        //
        //    Input, int FACTOR(NFACTOR), contains the factors.
        //    0.0 <= FACTOR(I).
        //
        //    Output, bool MULTICOEF_CHECK, is true if the parameters are legal.
        //
    {
        switch (nfactor)
        {
            case < 1:
                Console.WriteLine(" ");
                Console.WriteLine("MULTICOEF_CHECK - Warning!");
                Console.WriteLine("  NFACTOR < 1.");
                return false;
        }

        for (int i = 0; i < nfactor; i++)
        {
            switch (factor[i])
            {
                case < 0:
                    Console.WriteLine(" ");
                    Console.WriteLine("MULTICOEF_CHECK - Warning");
                    Console.WriteLine("  Factor[" + i + "] = " + factor[i] + "");
                    Console.WriteLine("  But this value must be nonnegative.");
                    return false;
            }
        }

        return true;
    }

    public static int multinomial_coef1(int nfactor, int[] factor)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTINOMIAL_COEF1 computes a multinomial coefficient.
        //
        //  Discussion:
        //
        //    The multinomial coefficient is a generalization of the binomial
        //    coefficient.  It may be interpreted as the number of combinations of
        //    N objects, where FACTOR(1) objects are indistinguishable of type 1,
        //    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
        //    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
        //
        //    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
        //
        //    The log of the gamma function is used, to avoid overflow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NFACTOR, the number of factors.
        //
        //    Input, int FACTOR[NFACTOR], contains the factors.
        //    0 <= FACTOR(I)
        //
        //    Output, int MULTINOMIAL_COEF1, the value of the multinomial coefficient.
        //
    {
        int i;
        //
        //  Each factor must be nonnegative.
        //
        for (i = 0; i < nfactor; i++)
        {
            switch (factor[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("MULTINOMIAL_COEF1 - Fatal error");
                    Console.WriteLine("  Factor " + i + " = " + factor[i] + "");
                    Console.WriteLine("  But this value must be nonnegative.");
                    return 1;
            }
        }

        //
        //  The factors sum to N.
        //
        int n = typeMethods.i4vec_sum(nfactor, factor);

        double arg = n + 1;
        double facn = Helpers.LogGamma(arg);

        for (i = 0; i < nfactor; i++)
        {
            arg = factor[i] + 1;
            double fack = Helpers.LogGamma(arg);
            facn -= fack;
        }

        int value = (int)typeMethods.r8_nint(Math.Exp(facn));

        return value;
    }

    public static int multinomial_coef2(int nfactor, int[] factor)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTINOMIAL_COEF2 computes a multinomial coefficient.
        //
        //  Discussion:
        //
        //    The multinomial coefficient is a generalization of the binomial
        //    coefficient.  It may be interpreted as the number of combinations of
        //    N objects, where FACTOR(1) objects are indistinguishable of type 1,
        //    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
        //    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
        //
        //    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
        //
        //    A direct method is used, which should be exact.  However, there
        //    is a possibility of intermediate overflow of the result.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NFACTOR, the number of factors.
        //
        //    Input, int FACTOR[NFACTOR], contains the factors.
        //    0 <= FACTOR(I)
        //
        //    Output, int MULTINOMIAL_COEF2, the value of the multinomial coefficient.
        //
    {
        int i;
        //
        //  Each factor must be nonnegative.
        //
        for (i = 0; i < nfactor; i++)
        {
            switch (factor[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("MULTINOMIAL_COEF2 - Fatal error!");
                    Console.WriteLine("  Factor " + i + " = " + factor[i] + "");
                    Console.WriteLine("  But this value must be nonnegative.");
                    return 1;
            }
        }

        int value = 1;
        int k = 0;

        for (i = 0; i < nfactor; i++)
        {
            int j;
            for (j = 1; j <= factor[i]; j++)
            {
                k += 1;
                value = value * k / j;
            }
        }

        return value;
    }
}