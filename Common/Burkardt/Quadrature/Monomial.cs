using System;
using Burkardt.Elliptic;
using Burkardt.Types;

namespace Burkardt.Quadrature
{
    public static class MonomialQuadrature
    {
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
            double exact;
            double quad;
            double quad_error;
            double[] value;
            //
            //  Get the exact value of the integral of the unscaled monomial.
            //
            exact = IntegralNS.Monomial.monomial_integral_mixed(dim_num, rule, alpha, beta, expon);
            //
            //  Evaluate the monomial at the quadrature points.
            //
            value = Monomial.monomial_value(dim_num, point_num, expon, x);
            //
            //  Compute the weighted sum and divide by the exact value.
            //
            quad = typeMethods.r8vec_dot(point_num, weight, value);
            //
            //  Error:
            //
            if (exact == 0.0)
            {
                quad_error = Math.Abs(quad - exact);
            }
            else
            {
                quad_error = Math.Abs(quad - exact) / Math.Abs(exact);
            }

            return quad_error;
        }
    }
}