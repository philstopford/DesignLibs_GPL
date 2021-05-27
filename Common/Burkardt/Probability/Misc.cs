using System;

namespace Burkardt.Probability
{
    public static class Misc
    {
        public static double sin_power_int ( double a, double b, int n )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIN_POWER_INT evaluates the sine power integral.
        //
        //  Discussion:
        //
        //    The function is defined by
        //
        //      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
        //
        //    The algorithm uses the following fact:
        //
        //      Integral sin^n ( t ) = (1/n) * (
        //        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
        //
        //  Licensing:
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
        //    Input, int N, the power of the sine function.
        //
        //    Output, double SIN_POWER_INT, the value of the integral.
        //
        {
            int mlo;
            double value;
            if ( n < 0 )
            {
                Console.WriteLine("");
                Console.WriteLine("SIN_POWER_INT - Fatal error!");
                Console.WriteLine("  Power N < 0.");
                return 1.0;
            }

            double sa = Math.Sin ( a );
            double sb = Math.Sin ( b );
            double ca = Math.Cos ( a );
            double cb = Math.Cos ( b );

            if ( ( n % 2 ) == 0 )
            {
                value = b - a;
                mlo = 2;
            }
            else
            {
                value = ca - cb;
                mlo = 3;
            }

            for (int m = mlo; m <= n; m = m + 2 )
            {
                value = ( ( double ) ( m - 1 ) * value
                            + Math.Pow ( sa, (m-1) ) * ca - Math.Pow ( sb, (m-1) ) * cb )
                        / ( double ) ( m );
            }

            return value;
        }
        
    }
}