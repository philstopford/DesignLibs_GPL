using System;
using System.Numerics;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static Complex cartesian_to_c8 ( double x, double y )

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

            c = new Complex ( x, y );

            return c;
        }

        public static Complex polar_to_c8 ( double r, double theta )

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

            c = new Complex ( r * Math.Cos ( theta ), r * Math.Sin ( theta ) );

            return c;
        }

    }
}