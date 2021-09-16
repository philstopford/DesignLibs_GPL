using Burkardt.Types;

namespace Burkardt.Square
{
    public static partial class MinimalRule
    {
        public static double[] smr06 ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SMR06 returns the SMR rule of degree 6.
            //
            //  Discussion:
            //
            //    DEGREE: 6
            //    POINTS CARDINALITY: 10
            //    NORM INF MOMS. RESIDUAL: 9.02056e-16
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
            //    Output, double *SMR06[3*10], the requested rule.
            //
        {
            int degree = 6;
            int order;
            double[] xw = {
                9.785155476851563e-01,  2.032663758845348e-01, 2.115144101443901e-01,
                8.002071174796732e-01, -8.068212789789568e-01, 2.633128596562718e-01,
                6.478447271179293e-01,  8.174437151345265e-01, 3.765392584231271e-01,
                4.527740405918522e-01, -2.317108088207500e-01, 7.128261121734308e-01,
                -4.775668731692919e-01,  9.766597077168154e-01, 1.915080888326935e-01,
                -7.569420164799104e-02,  4.409055401533251e-01, 7.647398979197457e-01,
                -1.241563248035317e-01, -8.409222499092199e-01, 4.422079284502612e-01,
                -9.014158913156406e-01,  5.483121720080090e-01, 2.801510395399764e-01,
                -6.560277628978444e-01, -2.866497529648795e-01, 6.536308666716537e-01,
                -9.537871517453275e-01, -8.861523430050633e-01, 1.035695381884487e-01 };
            double[] xw_copy;

            order = square_minimal_rule_order ( degree );
            xw_copy = typeMethods.r8mat_copy_new ( 3, order, xw );

            return xw_copy;
        }
    }
}