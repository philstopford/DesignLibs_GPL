using Burkardt.Types;

namespace Burkardt.IntegralNS;

public static class Epn_Lag
{

    public static double ep1_glg_monomial_integral(int expon, double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EP1_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EP1.
        //
        //  Discussion:
        //
        //    EP1_GLG is the interval [0,+oo) with generalized Laguerre weight function:
        //
        //      w(alpha;x) = x^alpha exp ( - x )
        //
        //    value = integral ( 0 <= x < +oo ) x^expon x^alpha exp ( - x ) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int EXPON, the exponent.
        //    0 <= EXPON.
        //
        //    Input, double ALPHA, the exponent of X in the weight function.
        //    -1.0 < ALPHA.
        //
        //    Output, double EP1_GLG_MONOMIAL_INTEGRAL, the value of the integral.
        //
    {
        double arg;
        double exact;

        arg = alpha + (expon + 1);

        exact = typeMethods.r8_gamma(arg);

        return exact;
    }

    public static double ep1_lag_monomial_integral(int expon)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EP1_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EP1.
        //
        //  Discussion:
        //
        //    EP1 is the interval [0,+oo) with exponential or Laguerre weight function:
        //
        //      w(x) = exp ( - x )
        //
        //    value = integral ( 0 <= x < oo ) x^expon exp ( - x ) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int EXPON, the exponent.
        //    0 <= EXPON.
        //
        //    Output, double EP1_LAG_MONOMIAL_INTEGRAL, the value of the integral.
        //
    {
        double value = 0;

        value = typeMethods.r8_factorial(expon);

        return value;
    }

    public static double epn_glg_monomial_integral(int n, int[] expon, double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EPN.
        //
        //  Discussion:
        //
        //    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
        //    Laguerre weight function:
        //
        //      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
        //
        //    value = integral ( EPN ) 
        //      product ( 1 <= i <= n ) x(I)^expon(i) x(i)^alpha exp ( - x(i) ) dx(i)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int EXPON[N], the exponents.
        //
        //    Input, double ALPHA, the exponent of X in the weight function.
        //    -1.0 < ALPHA.
        //
        //    Output, double EPN_GLG_MONOMIAL_INTEGRAL, the value of the integral.
        //
    {
        int i;
        double value = 0;
        double value2;

        value = 1.0;
        for (i = 0; i < n; i++)
        {
            value2 = ep1_glg_monomial_integral(expon[i], alpha);
            value *= value2;
        }

        return value;
    }

    public static double epn_lag_monomial_integral(int n, int[] expon)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EPN_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EPN.
        //
        //  Discussion:
        //
        //    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
        //    or Laguerre weight function:
        //
        //      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
        //
        //    value = integral ( EPN ) 
        //      product ( 1 <= i <= n ) x(I)^expon(i) exp ( -x(i) ) dx(i)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int EXPON(N), the exponents.
        //
        //    Output, double EPN_LAG_MONOMIAL_VALUE, the value of the integral.
        //
    {
        int i;
        double value = 0;
        double value2;

        value = 1.0;
        for (i = 0; i < n; i++)
        {
            value2 = ep1_lag_monomial_integral(expon[i]);
            value *= value2;
        }

        return value;
    }

}