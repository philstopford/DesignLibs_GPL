using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8_log_10(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LOG_10 returns the logarithm base 10 of the absolute value of an R8.
        //
        //  Discussion:
        //
        //    value = Log10 ( |X| )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose base 2 logarithm is desired.
        //    X should not be 0.
        //
        //    Output, double R8_LOG_10, the logarithm base 10 of the absolute
        //    value of X.  It should be true that |X| = 10^R_LOG_10.
        //
    {
        double value = x switch
        {
            0.0 => -r8_big(),
            _ => Math.Log10(Math.Abs(x))
        };

        return value;
    }

    public static double r8_log_2(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LOG_2 returns the logarithm base 2 of the absolute value of an R8.
        //
        //  Discussion:
        //
        //    value = Log ( |X| ) / Log ( 2.0 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose base 2 logarithm is desired.
        //    X should not be 0.
        //
        //    Output, double R8_LOG_2, the logarithm base 2 of the absolute
        //    value of X.  It should be true that |X| = 2^R_LOG_2.
        //
    {
        double value = x switch
        {
            0.0 => -r8_big(),
            _ => Math.Log(Math.Abs(x)) / Math.Log(2.0)
        };

        return value;
    }

    public static double r8_log_b(double x, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LOG_B returns the logarithm base B of an R8.
        //
        //  Discussion:
        //
        //    value = log ( |X| ) / log ( |B| )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose base B logarithm is desired.
        //    X should not be 0.
        //
        //    Input, double B, the base, which should not be 0, 1 or -1.
        //
        //    Output, double R8_LOG_B, the logarithm base B of the absolute
        //    value of X.  It should be true that |X| = |B|^R_LOG_B.
        //
    {
        double value = 0;

        switch (b)
        {
            case 0.0:
            case 1.0:
            case -1.0:
                value = -r8_big();
                break;
            default:
            {
                value = Math.Abs(x) switch
                {
                    0.0 => -r8_big(),
                    _ => Math.Log(Math.Abs(x)) / Math.Log(Math.Abs(b))
                };

                break;
            }
        }

        return value;
    }
}