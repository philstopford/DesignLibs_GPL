using Burkardt.Types;

namespace Burkardt.Square;

public static partial class MinimalRule
{
    public static double[] smr09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SMR09 returns the SMR rule of degree 9.
        //
        //  Discussion:
        //
        //    DEGREE: 9
        //    ROTATIONALLY INVARIANT: (X,Y),(-Y,X),(-X,-Y),(Y,-X).-
        //    POINTS CARDINALITY: 17
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
        //    Output, double *SMR09[3*17], the requested rule.
        //
    {
        const int degree = 9;
        double[] xw =
        {
            -6.306801197316687e-01, 9.688499663619777e-01, 8.887937817019871e-02,
            4.533398211356471e-01, -5.237358202144294e-01, 3.982824392620701e-01,
            8.526157293336624e-01, -7.620832819261733e-02, 2.690513376397807e-01,
            -7.502770999789007e-01, -9.279616459595696e-01, 1.120996021295965e-01,
            -9.688499663619777e-01, -6.306801197316687e-01, 8.887937817019871e-02,
            5.237358202144294e-01, 4.533398211356471e-01, 3.982824392620701e-01,
            7.620832819261733e-02, 8.526157293336624e-01, 2.690513376397807e-01,
            9.279616459595696e-01, -7.502770999789007e-01, 1.120996021295965e-01,
            6.306801197316687e-01, -9.688499663619777e-01, 8.887937817019871e-02,
            -4.533398211356471e-01, 5.237358202144294e-01, 3.982824392620701e-01,
            -8.526157293336624e-01, 7.620832819261733e-02, 2.690513376397807e-01,
            7.502770999789007e-01, 9.279616459595696e-01, 1.120996021295965e-01,
            9.688499663619777e-01, 6.306801197316687e-01, 8.887937817019871e-02,
            -5.237358202144294e-01, -4.533398211356471e-01, 3.982824392620701e-01,
            -7.620832819261733e-02, -8.526157293336624e-01, 2.690513376397807e-01,
            -9.279616459595696e-01, 7.502770999789007e-01, 1.120996021295965e-01,
            0.000000000000000e+00, 0.000000000000000e+00, 5.267489711934157e-01
        };

        int order = square_minimal_rule_order(degree);
        double[] xw_copy = typeMethods.r8mat_copy_new(3, order, xw);

        return xw_copy;
    }
}