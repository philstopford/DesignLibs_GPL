using System;
using System.Numerics;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static Complex cartesian_to_c8(double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CARTESIAN_TO_C8 converts a Cartesian form to a C8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 March 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, Y, the Cartesian form.
            //
            //    Output, Complex CARTESIAN_TO_C8, the complex number.
            //
        {
            Complex c;

            c = new Complex(x, y);

            return c;
        }

        public static Complex polar_to_c8(double r, double theta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLAR_TO_C8 converts a polar form to a C8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 March 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, THETA, the polar form.
            //
            //    Output, Complex POLAR_TO_C8, the complex number.
            //
        {
            Complex c;

            c = new Complex(r * Math.Cos(theta), r * Math.Sin(theta));

            return c;
        }

        public static Complex r8_csqrt(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CSQRT returns the complex square root of an R8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the number whose square root is desired.
            //
            //    Output, Complex R8_CSQRT, the square root of X:
            //
        {
            double argument = 0;
            double magnitude = 1;
            const double r8_pi = 3.141592653589793;
            Complex value;

            if (0.0 < x)
            {
                magnitude = x;
                argument = 0.0;
            }
            else if (0.0 == x)
            {
                magnitude = 0.0;
                argument = 0.0;
            }
            else if (x < 0.0)
            {
                magnitude = -x;
                argument = r8_pi;
            }

            magnitude = Math.Sqrt(magnitude);
            argument = argument / 2.0;

            value = magnitude * new Complex(Math.Cos(argument), Math.Sin(argument));

            return value;
        }

    }
}