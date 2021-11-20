using System;
using Burkardt.IntegralNS;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class CN_Leg
{
    public static void cn_leg_01_1 ( int n, int o, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_LEG_01_1 implements the midpoint rule for region CN_LEG.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 1.
        //
        //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
        //
        //      w(x) = 1. 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int O, the order.
        //
        //    Output, double X[N*O], the abscissas.
        //
        //    Output, double W[O], the weights.
        //
    {
        const int expon = 0;
        double volume = C1.c1_leg_monomial_integral ( expon );
        volume = Math.Pow ( volume, n );

        typeMethods.r8vec_zero ( n * o, ref x );

        int k = - 1;
        //
        //  1 point.
        //
        k += 1;
        w[k] = volume;
    }

    public static int cn_leg_01_1_size ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CN_LEG_01_1_SIZE sizes the midpoint rule for region CN_LEG.
        //
        //  Discussion:
        //
        //    The rule has order O = 1.
        //
        //    The rule has precision P = 1.
        //
        //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
        //
        //      w(x) = 1. 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Output, int CN_LEG_01_1_SIZE, the order.
        //
    {
        const int o = 1;

        return o;
    }

}