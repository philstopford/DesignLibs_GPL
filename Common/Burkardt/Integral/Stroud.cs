using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS;

public static class Stroud
{
    public static double en_her_monomial_integral ( int n, int[] alpha )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_HER_MONOMIAL_INTEGRAL evaluates monomial integrals in EN_HER.
        //
        //  Discussion:
        //
        //    ALPHA is the set of polynomial exponents.
        //
        //    EN_HER is the entire N-dimensional space with weight function
        //
        //      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
        //
        //    The integral to be evaluated is
        //
        //      value = integral ( EN ) x(1)^alpha(1) * x(2)^alpha(2) * ... 
        //        * x(n)^alpha(n) * w(x) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int ALPHA[N], the polynomial exponents.
        //    0 <= ALPHA[*].
        //
        //    Output, double EN_HER_MONOMIAL_INTEGRAL, the value of the integral.
        //
    {
        int i;

        for ( i = 0; i < n; i++ )
        {
            switch (alpha[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("EN_HER_MONOMIAL_INTEGRAL - Fatal error//");
                    Console.WriteLine("  ALPHA[" + i + "] < 0.");
                    return 1;
            }
        }

        double value = 1.0;
        for ( i = 0; i < n; i++ )
        {
            if ( alpha[i] % 2 == 1 )
            {
                value = 0.0;
                break;
            }

            double arg = (alpha[i] + 1) / 2.0;
            value *= typeMethods.r8_gamma ( arg );
        }

        return value;
    }
    public static double en_r2_monomial_integral(int n, int[] alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EN_R2_MONOMIAL_INTEGRAL evaluates monomial integrals in EN_R2.
        //
        //  Discussion:
        //
        //    ALPHA is the set of polynomial exponents.
        //
        //    EN_R2 is the entire N-dimensional space with weight function
        //
        //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
        //
        //    The integral to be evaluated is
        //
        //      value = integral ( EN ) x(1)^alpha(1) * x(2)^alpha(2) * ... 
        //        * x(n)^alpha(n) * w(x) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int ALPHA[N], the polynomial exponents.
        //    0 <= ALPHA[*].
        //
        //    Output, double EN_R2_MONOMIAL_INTEGRAL, the value of the integral.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            switch (alpha[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  ALPHA[" + i + "] < 0.");
                    return 1;
            }
        }

        double value = 1.0;
        for (i = 0; i < n; i++)
        {
            if (alpha[i] % 2 == 1)
            {
                value = 0.0;
                break;
            }

            double arg = (alpha[i] + 1) / 2.0;
            value *= typeMethods.r8_gamma(arg);
        }

        return value;
    }

}