using Burkardt.Types;

namespace Burkardt.Square
{
    public static partial class MinimalRule
    {
        public static double[] smr05 ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SMR05 returns the SMR rule of degree 5.
            //
            //  Discussion:
            //
            //    DEGREE: 5
            //    SYMMETRY: (X,Y), (-X,-Y)
            //    POINTS CARDINALITY: 7
            //    NORM INF MOMS. RESIDUAL: 1.66533e-16
            //    SUM NEGATIVE WEIGHTS: 0
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
            //    Output, double *SMR05[3*7], the requested rule.
            //
        {
            int degree = 5;
            int order;
            double[] xw = {
                -9.660917830792960e-01,  0.000000000000000e+00, 3.174603174603175e-01,
                -5.773502691896258e-01, -7.745966692414834e-01, 5.555555555555556e-01,
                -5.773502691896258e-01,  7.745966692414834e-01, 5.555555555555556e-01,
                0.000000000000000e+00,  0.000000000000000e+00, 1.142857142857143e+00,
                5.773502691896258e-01,  7.745966692414834e-01, 5.555555555555556e-01,
                5.773502691896258e-01, -7.745966692414834e-01, 5.555555555555556e-01,
                9.660917830792960e-01,  0.000000000000000e+00, 3.174603174603175e-01 };
            double[] xw_copy;

            order = square_minimal_rule_order ( degree );
            xw_copy = typeMethods.r8mat_copy_new ( 3, order, xw );

            return xw_copy;
        }
    }
}