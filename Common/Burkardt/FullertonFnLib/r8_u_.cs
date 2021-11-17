using System;

namespace Burkardt.FullertonFnLib;

public static partial class FullertonLib
{
    public static void r8_upak ( double x, ref double y, ref int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_UPAK unpacks an R8 into a mantissa and exponent.
        //
        //  Discussion:
        //
        //    This function unpacks a floating point number x so that
        //
        //      x = y * 2.0^n
        //
        //    where
        //
        //      0.5 <= abs ( y ) < 1.0 .
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 September 2011
        //
        //  Author:
        //
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, double X, the number to be unpacked.
        //
        //    Output, double &Y, the mantissa.
        //
        //    Output, int &N, the exponent.
        //
    {
        double absx;

        absx = Math.Abs ( x );
        n = 0;
        y = 0.0;

        switch (x)
        {
            case 0.0:
                return;
        }

        while ( absx < 0.5 )
        {
            n -= 1;
            absx *= 2.0;
        }

        while ( 1.0 <= absx )
        {
            n += 1;
            absx *= 0.5;
        }

        y = x switch
        {
            < 0.0 => -absx,
            _ => +absx
        };
    }

}