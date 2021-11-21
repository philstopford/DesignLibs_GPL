using System;

namespace Burkardt.Function;

public static class Franke
{
    public static double franke ( double x, double y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FRANKE returns the value of the Franke function #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Richard Franke,
        //    Scattered Data Interpolation: Tests of Some Methods,
        //    Mathematics of Computation,
        //    Volume 38, Number 157, January 1982, pages 181-200.
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evalution points.
        //
        //    Output, double FRANKE, the function values.
        //
    {
        double f = 0.75 * Math.Exp ( 
                       - ( Math.Pow ( 9.0 * x - 2.0, 2 )
                           + Math.Pow ( 9.0 * y - 2.0, 2 ) ) / 4.0 ) 
                   + 0.75 * Math.Exp ( 
                       - Math.Pow ( 9.0 * x + 1.0, 2 ) / 49.0 
                       - ( 9.0 * y + 1.0 )      / 10.0 )       
                   + 0.5  * Math.Exp ( 
                       - ( Math.Pow ( 9.0 * x - 7.0, 2 ) 
                           + Math.Pow ( 9.0 * y - 3.0, 2 ) ) / 4.0 ) 
                   - 0.2  * Math.Exp ( 
                       - Math.Pow ( 9.0 * x - 4.0, 2 ) 
                       - Math.Pow ( 9.0 * y - 7.0, 2 ) );

        return f;
    }
}