using System;
using Burkardt.Types;

namespace Burkardt.Interpolation;

public static class Smolyak
{
    public static void smolyak_coefficients(int l_max, int m, ref int[] c, ref int[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SMOLYAK_COEFFICIENTS returns the Smolyak coefficients and counts.
        //
        //  Discussion:
        //
        //    The Smolyak sparse interpolant can be written as:
        //
        //      A(L,M)(X) = sum ( L-M+1 <= |L| <= L_max ) 
        //        C(|L|) * g(l1)(x1) * g(l2)(x2) * ... * g(lm)(xm).
        //
        //    where:
        //
        //    * L=(l1,l2,...,lm) is a vector of M nonnegative integers;
        //    * |L| is the sum of the entries of L;
        //    * X=(x1,x2,...,xm) is an M-dimensional point in a product space;
        //    * g(i)(xj) is the i-th 1-d interpolation function in dimension j;
        //
        //    Note that:
        //
        //    * W(|L|) will represent the number of distinct interpolants for which
        //      the sublevel, or sum of the L vector entries, is |L|;
        //
        //    * the coefficients C and counts W will be zero for sublevels 
        //      0 through L_MAX - M (and MATLAB indices 1 through L_MAX-M+1).
        //
        //    * it will be the case that W' * C = 1, essentially because the interpolant
        //      to the identity function must be the identity function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int L_MAX, the (maximum) level.
        //    0 <= L_MA.
        //
        //    Input, int M, the spatial dimension.
        //    1 <= M.
        //
        //    Output, int C[L_MAX+1], the coefficients for objects 
        //    at sublevels 0 through L_MAX.
        //
        //    Output, int W[L_MAX+1], the number of objects at 
        //    sublevels 0 through L_MAX.
        //
    {
        int l;
        int l_min;

        l_min = Math.Max(l_max - m + 1, 0);

        for (l = 0; l < l_min; l++)
        {
            c[l] = 0;
        }

        for (l = l_min; l <= l_max; l++)
        {
            c[l] = typeMethods.i4_mop(l_max - l) * typeMethods.i4_choose(m - 1, l_max - l);
        }

        for (l = 0; l < l_min; l++)
        {
            w[l] = 0;
        }

        for (l = l_min; l <= l_max; l++)
        {
            w[l] = typeMethods.i4_choose(l + m - 1, m - 1);
        }
    }

}