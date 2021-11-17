using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8_factorial(int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FACTORIAL computes the factorial of N.
        //
        //  Discussion:
        //
        //    factorial ( N ) = product ( 1 <= I <= N ) I
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 January 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument of the factorial function.
        //    If N is less than 1, the function value is returned as 1.
        //
        //    Output, double R8_FACTORIAL, the factorial of N.
        //
    {
        double value = 1.0;

        for (int i = 1; i <= n; i++)
        {
            value *= i;
        }

        return value;
    }

    public static double r8_factorial_stirling ( int n )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FACTORIAL_STIRLING computes Stirling's approximation to N!.
        //
        //  Discussion:
        //
        //    N! = Product ( 1 <= I <= N ) I
        //
        //    Stirling ( N ) = sqrt ( 2 * PI * N ) * ( N / E )^N * E^(1/(12*N) )
        //
        //    This routine returns the raw approximation for all nonnegative
        //    values of N.  If N is less than 0, the value is returned as 0,
        //    and if N is 0, the value of 1 is returned.  In all other cases,
        //    Stirling's formula is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument of the function.
        //
        //    Output, double R8_FACTORIAL_STIRLING, an approximation to N!.
        //
    {
        double value = n switch
        {
            < 0 => 0.0,
            0 => 1.0,
            _ => Math.Sqrt(2.0 * Math.PI * n) * Math.Pow(n / Math.E, n) * Math.Exp(1.0 / (12 * n))
        };

        return value;
    }

    public static double r8_factorial_log ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FACTORIAL_LOG computes the natural logarithm of the factorial function.
        //
        //  Discussion:
        //
        //    LOG ( FACTORIAL ( N ) )
        //      = LOG ( product ( 1 <= I <= N ) I )
        //      = sum ( ( 1 <= I <= N ) LOG ( I ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument of the factorial function.
        //    If N is less than 1, R8_FACTORIAL_LOG is returned as 0.
        //
        //    Output, double R8_FACTORIAL_LOG, the logarithm of the factorial of N.
        //
    {
        int i;

        double value = 0.0;

        for ( i = 1; i <= n; i++ )
        {
            value += Math.Log ( i );
        }

        return value;
    }

    public static double r8_factorial2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FACTORIAL2 computes the double factorial function.
        //
        //  Discussion:
        //
        //    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
        //                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
        //
        //  Example:
        //
        //     N Value
        //
        //     0     1
        //     1     1
        //     2     2
        //     3     3
        //     4     8
        //     5    15
        //     6    48
        //     7   105
        //     8   384
        //     9   945
        //    10  3840
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument of the double factorial
        //    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
        //
        //    Output, double R8_FACTORIAL2, the value of Factorial2(N).
        //
    {
        double value = 1.0;

        switch (n)
        {
            case < 1:
                return value;
        }

        int n_copy = n;

        while (1 < n_copy)
        {
            value *= n_copy;
            n_copy -= 2;
        }

        return value;
    }


}