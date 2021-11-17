using System;

namespace Burkardt.CircleNS;

public static class QuadratureRule
{
    public static void circle_rule ( int nt, ref double[] w, ref double[] t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_RULE computes a quadrature rule for the unit circle.
        //
        //  Discussion:
        //
        //    The unit circle is the region:
        //
        //      x * x + y * y = 1.
        //
        //    The integral I(f) is then approximated by
        //
        //      Q(f) = 2 * Math.PI * sum ( 1 <= i <= NT ) W(i) * F ( cos(T(i)), sin(T(i)) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NT, the number of angles to use.
        //
        //    Output, double W[NT], the weights for the rule.
        //
        //    Output, double T[NT], the angles for the rule.
        //
    {
        for ( int it = 0; it < nt; it++ )
        {
            w[it] = 1.0 / nt;
            t[it] = 2.0 * Math.PI * it / nt;
        }
    }
}