using Burkardt.Types;

namespace Burkardt.Square;

public static partial class MinimalRule
{
    public static double[] smr01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SMR01 returns the SMR rule of degree 1.
        //
        //  Discussion:
        //
        //    DEGREE: 1
        //    SYMM.: (X,Y), (-X,-Y)
        //    POINTS CARDINALITY: 1
        //    NORM INF MOMS. RESIDUAL: 0.00000e+00,
        //    SUM NEGATIVE WEIGHTS: 0.00000e+00,
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 February 2018
        //
        //  Author:
        //
        //    Original MATLAB version by Mattia Festa, Alvise Sommariva,
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Mattia Festa, Alvise Sommariva,
        //    Computing almost minimal formulas on the square,
        //    Journal of Computational and Applied Mathematics,
        //    Volume 17, Number 236, November 2012, pages 4296-4302.
        //
        //  Parameters:
        //
        //    Output, double *SMR01[3*1], the requested rule.
        //
    {
        int degree = 1;
        int order;
        double[] xw =
        {
            0.000000000000000e+00, 0.000000000000000e+00, 4.000000000000000e+00
        };
        double[] xw_copy;

        order = square_minimal_rule_order(degree);
        xw_copy = typeMethods.r8mat_copy_new(3, order, xw);

        return xw_copy;
    }
}