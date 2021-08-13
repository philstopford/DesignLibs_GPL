using System;

namespace Burkardt.Function
{
    public static class Gudermannian
    {
        public static double gud ( double x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GUD evaluates the Gudermannian function.
            //
            //  Definition:
            //
            //    The Gudermannian function relates the hyperbolic and trigonometric
            //    functions.  For any argument X, there is a corresponding value
            //    GAMMA so that
            //
            //      sinh(x) = tan(gamma).
            //
            //    The value GAMMA is called the Gudermannian of X.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument of the Gudermannian.
            //
            //    Output, double GUD, the value of the Gudermannian.
            //
        {
            double value;

            value = 2.0 * Math.Atan ( Math.Tanh ( 0.5 * x ) );

            return value;
        }
    }
}