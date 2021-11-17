using Burkardt.Types;

namespace Burkardt.Square;

public static partial class MinimalRule
{
    public static double[] smr07 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SMR07 returns the SMR rule of degree 7.
        //
        //  Discussion:
        //
        //    DEGREE: 7
        //    ROTATIONALLY INVARIANT: (X,Y),(-Y,X),(-X,-Y),(Y,-X).
        //    POINTS CARDINALITY: 12
        //    NORM INF MOMS. RESIDUAL: 2.77556e-16
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
        //    Output, double *SMR07[3*12], the requested rule.
        //
    {
        int degree = 7;
        int order;
        double[] xw = {
            3.938313610187890e-01, -3.669979343404274e-01, 5.209286223023044e-01,
            8.019684762065925e-01,  8.099259297580919e-01, 2.374024056137991e-01,
            -1.832732766424723e-02,  9.259609259310837e-01, 2.416689720838965e-01,
            3.669979343404274e-01,  3.938313610187890e-01, 5.209286223023044e-01,
            -8.099259297580919e-01,  8.019684762065925e-01, 2.374024056137991e-01,
            -9.259609259310837e-01, -1.832732766424723e-02, 2.416689720838965e-01,
            -3.938313610187890e-01,  3.669979343404274e-01, 5.209286223023044e-01,
            -8.019684762065925e-01, -8.099259297580919e-01, 2.374024056137991e-01,
            1.832732766424723e-02, -9.259609259310837e-01, 2.416689720838965e-01,
            -3.669979343404274e-01, -3.938313610187890e-01, 5.209286223023044e-01,
            8.099259297580919e-01, -8.019684762065925e-01, 2.374024056137991e-01,
            9.259609259310837e-01,  1.832732766424723e-02, 2.416689720838965e-01 };
        double[] xw_copy;

        order = square_minimal_rule_order ( degree );
        xw_copy = typeMethods.r8mat_copy_new ( 3, order, xw );

        return xw_copy;
    }
}