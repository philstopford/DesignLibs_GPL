using Burkardt.Types;

namespace Burkardt.Square;

public static partial class MinimalRule
{
    public static double[] smr02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SMR02 returns the SMR rule of degree 2.
        //
        //  Discussion:
        //
        //    DEGREE: 2
        //    POINTS CARDINALITY: 3
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
        //    Output, double *SMR02[3*3], the requested rule.
        //
    {
        const int degree = 2;
        double[] xw =
        {
            6.519542382482019e-01, 3.623444315022428e-01, 1.498681096511322e+00,
            -8.070876861583226e-01, 5.322309038833022e-01, 1.051530016968792e+00,
            -8.856086946552445e-02, -7.605904084126465e-01, 1.449788886519887e+00
        };

        int order = square_minimal_rule_order(degree);
        double[] xw_copy = typeMethods.r8mat_copy_new(3, order, xw);

        return xw_copy;
    }
}