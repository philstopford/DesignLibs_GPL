using Burkardt.Types;

namespace Burkardt.Square;

public static partial class MinimalRule
{
    public static double[] smr03 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SMR03 returns the SMR rule of degree 3.
        //
        //  Discussion:
        //
        //    DEGREE: 3
        //    SYMMETRY: (X,Y), (-Y,X), (-X,-Y), (Y,-X)
        //    POINTS CARDINALITY: 4
        //    NORM INF MOMS. RESIDUAL: 2.22045e-16
        //    SUM NEGATIVE WEIGHTS: 0.00000e+00,
        // 
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 February 2018
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
        //    Output, double *SMR03[3*4], the requested rule.
        //
    {
        int degree = 3;
        int order;
        double[] xw = {
            -5.773502691896257e-01, -5.773502691896257e-01, 1.000000000000000e+00,
            -5.773502691896257e-01,  5.773502691896257e-01, 1.000000000000000e+00,
            5.773502691896257e-01,  5.773502691896257e-01, 1.000000000000000e+00,
            5.773502691896257e-01, -5.773502691896257e-01, 1.000000000000000e+00 };
        double[] xw_copy;

        order = square_minimal_rule_order ( degree );
        xw_copy = typeMethods.r8mat_copy_new ( 3, order, xw );

        return xw_copy;
    }
        
}