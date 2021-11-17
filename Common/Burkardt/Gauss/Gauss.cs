using System;

namespace Burkardt;

public static class Gauss
{
    public static double gauss ( double t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAUSS returns the area of the lower tail of the normal curve.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T, the evaluation point.
        //
        //    Output, double GAUSS, the lower normal tail area.
        //
    {
        double value = 0;

        value = ( 1.0 + Helpers.Erf ( t / Math.Sqrt ( 2.0 ) ) ) / 2.0;

        return value;
    }
}