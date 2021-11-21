using System;

namespace Burkardt.IntegralNS;

public static class SinPower
{
    public static double sin_power_int ( double a, double b, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIN_POWER_INT evaluates the Math.Sine Math.Power integral.
        //
        //  Discussion:
        //
        //    The function is defined by
        //
        //      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( Math.Sin ( t ))^n dt
        //
        //    The algorithm uses the following fact:
        //
        //      Integral Math.Sin^n ( t ) = (1/n) * (
        //        Math.Sin^(n-1)(t) * Math.Cos(t) + ( n-1 ) * Integral Math.Sin^(n-2) ( t ) dt )
        //
        //  LicenMath.Sing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, double A, B, the limits of integration.
        //
        //    Input, integer N, the Math.Power of the Math.Sine function.
        //
        //    Output, double SIN_POWER_INT, the value of the integral.
        //
    {
        int m;
        int mlo;
        double value;

        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("SIN_POWER_INT - Fatal error!");
                Console.WriteLine("  Power N < 0.");
                return 1;
        }

        double sa = Math.Sin ( a );
        double sb = Math.Sin ( b );
        double ca = Math.Cos ( a );
        double cb = Math.Cos ( b );

        switch (n % 2)
        {
            case 0:
                value = b - a;
                mlo = 2;
                break;
            default:
                value = ca - cb;
                mlo = 3;
                break;
        }

        for ( m = mlo; m <= n; m += 2 )
        {
            value = ( (m - 1) * value 
                        + Math.Pow ( sa, m-1 ) * ca - Math.Pow ( sb, m-1 ) * cb ) 
                    / m;
        }

        return value;
    }
}