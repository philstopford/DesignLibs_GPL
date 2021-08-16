using System;

namespace Burkardt.Function
{
    public static class Gudermannian
    {
        public static double agud ( double g )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AGUD evaluates the inverse Gudermannian function.
            //
            //  Discussion:
            //
            //    The Gudermannian function relates the hyperbolic and trigonomentric
            //    functions.  For any argument X, there is a corresponding value
            //    G so that
            //
            //      SINH(X) = TAN(G).
            //
            //    This value G(X) is called the Gudermannian of X.  The inverse
            //    Gudermannian function is given as input a value G and computes
            //    the corresponding value X.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double G, the value of the Gudermannian.
            //
            //    Output, double AGUD, the argument of the Gudermannian.
            //
        {
            const double r8_pi = 3.141592653589793;
            double value;

            value = Math.Log ( Math.Tan ( 0.25 * r8_pi + 0.5 * g ) );

            return value;
        }
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