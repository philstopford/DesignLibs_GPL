using Burkardt.Types;

namespace Burkardt.Square;

public static partial class MinimalRule
{
    public static double[] smr04 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SMR04 returns the SMR rule of degree 4.
        //
        //  Discussion:
        //
        //    DEGREE: 4
        //    POINTS CARDINALITY: 6
        //    NORM INF MOMS. RESIDUAL: 8.88178e-16
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
        //    Output, double *SMR04[3*6], the requested rule.
        //
    {
        int degree = 4;
        int order;
        double[] xw = {
            9.298664473826397e-01,  6.361197473108544e-02,  4.979283660841867e-01,
            -7.329012618874027e-01,  5.903145258425608e-01,  6.883081069413867e-01,
            3.895446419719248e-01,  8.325323327063485e-01,  6.340849824642651e-01,
            5.169294362169509e-01, -8.804002381721473e-01,  4.856847322376568e-01,
            -5.223159975544114e-02, -1.540167862605936e-01,  1.180646405191258e+00,
            -7.693563599017555e-01, -6.943687766134327e-01,  5.133474070812475e-01 };
        double[] xw_copy;

        order = square_minimal_rule_order ( degree );
        xw_copy = typeMethods.r8mat_copy_new ( 3, order, xw );

        return xw_copy;
    }
}