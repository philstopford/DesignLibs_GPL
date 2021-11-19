﻿using System;
using Burkardt.IntegralNS;
using Burkardt.Types;
using Monomial = Burkardt.MonomialNS.Monomial;

namespace Burkardt.Quadrature;

public static class MonomialQuadrature
{
    public static double monomial_quadrature_chebyshev1 ( int expon, int order, double[] w, 
            double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_QUADRATURE_CHEBYSHEV1 approximates a Chebyshev type 1 monomial integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int EXPON, the exponent.
        //
        //    Input, int ORDER, the number of points in the rule.
        //
        //    Input, double W[ORDER], the quadrature weights.
        //
        //    Input, double X[ORDER], the quadrature points.
        //
        //    Output, double MONOMIAL_QUADRATURE_CHEBYSHEV1, the quadrature error.
        //
    {
        int i;
        //
        //  Get the exact value of the integral.
        //
        double exact = Integral.chebyshev1_integral ( expon );
        //
        //  Evaluate the monomial at the quadrature points
        //  and compute the weighted sum.
        //
        double quad = 0.0;
        for ( i = 0; i < order; i++ )
        {
            quad += w[i] * Math.Pow ( x[i], expon );
        }

        double quad_error = exact switch
        {
            //
            //  Error:
            //
            0.0 => Math.Abs(quad - exact),
            _ => Math.Abs((quad - exact) / exact)
        };

        return quad_error;
    }
        
    public static double monomial_quadrature_chebyshev2 ( int expon, int order, double[] w, 
            double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_QUADRATURE_CHEBYSHEV2 approximates a Chebyshev type 2 monomial integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int EXPON, the exponent.
        //
        //    Input, int ORDER, the number of points in the rule.
        //
        //    Input, double W[ORDER], the quadrature weights.
        //
        //    Input, double X[ORDER], the quadrature points.
        //
        //    Output, double MONOMIAL_QUADRATURE_CHEBYSHEV2, the quadrature error.
        //
    {
        int i;
        //
        //  Get the exact value of the integral.
        //
        double exact = Integral.chebyshev2_integral ( expon );
        //
        //  Evaluate the monomial at the quadrature points
        //  and compute the weighted sum.
        //
        double quad = 0.0;
        for ( i = 0; i < order; i++ )
        {
            quad += w[i] * Math.Pow ( x[i], expon );
        }

        double quad_error = exact switch
        {
            //
            //  Error:
            //
            0.0 => Math.Abs(quad - exact),
            _ => Math.Abs((quad - exact) / exact)
        };

        return quad_error;
    }
        
    public static double monomial_quadrature ( int expon, int order, double[] w, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_QUADRATURE applies a quadrature rule to a monomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int EXPON, the exponent.
        //
        //    Input, int ORDER, the number of points in the rule.
        //
        //    Input, double W[ORDER], the quadrature weights.
        //
        //    Input, double X[ORDER], the quadrature points.
        //
        //    Output, double MONOMIAL_QUADRATURE, the quadrature error.
        //
    {
        int i;
        //
        //  Get the exact value of the integral.
        //
        double exact = 1.0 / (expon + 1);
        //
        //  Evaluate the monomial at the quadrature points.
        //
        double quad = 0.0;
        for ( i = 0; i < order; i++ )
        {
            quad += w[i] * Math.Pow ( x[i], expon );
        }
        //
        //  Relative error:
        //
        double quad_error = Math.Abs ( quad - exact ) / exact;

        return quad_error;
    }

    public static double monomial_quadrature(int dim_num, int[] expon, int point_num,
            double[] weight, double[] x, int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_QUADRATURE applies a quadrature rule to a monomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int EXPON[DIM_NUM], the exponents.
        //
        //    Input, int POINT_NUM, the number of points in the rule.
        //
        //    Input, double WEIGHT[POINT_NUM], the quadrature weights.
        //
        //    Input, double X[DIM_NUM*POINT_NUM], the quadrature points.
        //
        //    Input, int RULE, the index of the rule.
        //    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
        //    2, "F1", Fejer 1 Open Fully Nested rule.
        //    3, "F2", Fejer 2 Open Fully Nested rule.
        //    4, "GP", Gauss Patterson Open Fully Nested rule.
        //    5, "GL", Gauss Legendre Open Weakly Nested rule.
        //    6, "GH", Gauss Hermite Open Weakly Nested rule.
        //    7, "LG", Gauss Laguerre Open Non Nested rule.
        //
        //    Output, double MONOMIAL_QUADRATURE, the quadrature error.
        //
    {
        double exact;
        int point;
        switch (rule)
        {
            //
            //  Get the exact value of the integral of the unscaled monomial.
            //
            case >= 1 and <= 5:
                exact = IntegralNS.Monomial.monomial_integral_legendre(dim_num, expon);
                break;
            case 6:
                exact = IntegralNS.Monomial.monomial_integral_hermite(dim_num, expon);
                break;
            case 7:
                exact = IntegralNS.Monomial.monomial_integral_laguerre(dim_num, expon);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("MONOMIAL_QUADRATURE - Fatal error!");
                Console.WriteLine("  Unrecognized value of RULE.");
                return 1;
        }

        //
        //  Evaluate the monomial at the quadrature points.
        //
        double[] value = Monomial.monomial_value(dim_num, point_num, x, expon);
        //
        //  Compute the quadrature sum.
        //
        double quad = 0.0;
        for (point = 0; point < point_num; point++)
        {
            quad += weight[point] * value[point];
        }

        double quad_error = exact switch
        {
            //
            //  Absolute error if EXACT = 0, relative error otherwise:
            //
            0.0 => Math.Abs(quad - exact),
            _ => Math.Abs(quad - exact) / Math.Abs(exact)
        };

        return quad_error;
    }

    public static double monomial_quadrature ( int dim_num, int[] expon, int point_num, 
            double[] weight, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_QUADRATURE applies a quadrature rule to a monomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int EXPON[DIM_NUM], the exponents.
        //
        //    Input, int POINT_NUM, the number of points in the rule.
        //
        //    Input, double WEIGHT[POINT_NUM], the quadrature weights.
        //
        //    Input, double X[DIM_NUM*POINT_NUM], the quadrature points.
        //
        //    Output, double MONOMIAL_QUADRATURE, the quadrature error.
        //
    {
        int point;
        //
        //  Get the exact value of the integral of the unscaled monomial.
        //
        double scale = IntegralNS.Monomial.monomial_int01 ( dim_num, expon );
        //
        //  Evaluate the monomial at the quadrature points.
        //
        double[] value = Monomial.monomial_value ( dim_num, point_num, x, expon );
        //
        //  Compute the weighted sum and divide by the exact value.
        //
        double quad = 0.0;
        for ( point = 0; point < point_num; point++ )
        {
            quad += weight[point] * value[point];
        }
        quad /= scale;
        //
        //  Error:
        //
        const double exact = 1.0;
        double quad_error = Math.Abs ( quad - exact );
            
        return quad_error;
    }
    public static double monomial_quadrature(int dim_num, int point_num, int[] rule,
            double[] alpha, double[] beta, int[] expon, double[] weight, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_QUADRATURE applies a quadrature rule to a monomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int POINT_NUM, the number of points in the rule.
        //
        //    Input, int RULE[DIM_NUM], the component rules.
        //    1, Gauss-Legendre rule on [-1,+1];
        //    2, Gauss-Jacobi rule on [-1,+1];
        //    3, Gauss-Laguerre rule on [0,+oo);
        //    4, Generalized Gauss-Laguerre rule on [0,+oo);
        //    5, Gauss-Hermite rule on (-oo,+oo);
        //    6, Generalized Gauss-Hermite rule on (-oo,+oo).
        //
        //    Input, double ALPHA[DIM_NUM], BETA[DIM_NUM], parameters that
        //    may be needed for Jacobi, Generalized-Laguerre, or Generalized Hermite rules.
        //
        //    Input, int EXPON[DIM_NUM], the exponents.
        //
        //    Input, double WEIGHT[POINT_NUM], the quadrature weights.
        //
        //    Input, double X[DIM_NUM*POINT_NUM], the quadrature points.
        //
        //    Output, double MONOMIAL_QUADRATURE, the quadrature error.
        //
    {
        //
        //  Get the exact value of the integral of the unscaled monomial.
        //
        double exact = IntegralNS.Monomial.monomial_integral_mixed(dim_num, rule, alpha, beta, expon);
        //
        //  Evaluate the monomial at the quadrature points.
        //
        double[] value = Monomial.monomial_value(dim_num, point_num, expon, x);
        //
        //  Compute the weighted sum and divide by the exact value.
        //
        double quad = typeMethods.r8vec_dot(point_num, weight, value);
        double quad_error = exact switch
        {
            //
            //  Error:
            //
            0.0 => Math.Abs(quad - exact),
            _ => Math.Abs(quad - exact) / Math.Abs(exact)
        };

        return quad_error;
    }
}