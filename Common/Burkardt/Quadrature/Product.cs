using System;
using System.Linq;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class Product
{
    public static void product_mixed_growth_weight(int dim_num, int[] order_1d, int order_nd,
            int[] rule, int[] np, double[] p,
            Func<int, int, double[], double[], double[]>[] gw_compute_weights,
            ref double[] weight_nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRODUCT_MIXED_GROWTH_WEIGHT computes the weights of a mixed product rule.
        //
        //  Discussion:
        //
        //    This routine computes the weights for a quadrature rule which is
        //    a product of 1D rules of varying order and kind.
        //
        //    The user must preallocate space for the output array WEIGHT_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
        //
        //    Input, int ORDER_ND, the order of the product rule.
        //
        //    Input, int RULE[DIM_NUM], the rule in each dimension.
        //     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
        //     2, "F2",  Fejer Type 2, Open Fully Nested.
        //     3, "GP",  Gauss Patterson, Open Fully Nested.
        //     4, "GL",  Gauss Legendre, Open Weakly Nested.
        //     5, "GH",  Gauss Hermite, Open Weakly Nested.
        //     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
        //     7, "LG",  Gauss Laguerre, Open Non Nested.
        //     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
        //     9, "GJ",  Gauss Jacobi, Open Non Nested.
        //    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
        //    11, "UO",  User supplied Open, presumably Non Nested.
        //    12, "UC",  User supplied Closed, presumably Non Nested.
        //
        //    Input, int NP[DIM_NUM], the number of parameters used by each rule.
        //
        //    Input, double P[sum(NP[*])], the parameters needed by each rule.
        //
        //    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
        //    an array of pointers to functions which return the 1D quadrature weights 
        //    associated with each spatial dimension for which a Golub Welsch rule 
        //    is used.
        //
        //    Output, double WEIGHT_ND[ORDER_ND], the product rule weights.
        //
    {
        int dim;
        int i;
        typeMethods.r8vecDPData data = new();

        for (i = 0; i < order_nd; i++)
        {
            weight_nd[i] = 1.0;
        }

        int p_index = 0;

        for (dim = 0; dim < dim_num; dim++)
        {
            double[] weight_1d = new double[order_1d[dim]];

            switch (rule[dim])
            {
                case 1:
                    ClenshawCurtis.clenshaw_curtis_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                case 2:
                    Fejer2.fejer2_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                case 3:
                    PattersonQuadrature.patterson_lookup_weights_np(
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                case 4:
                    Legendre.QuadratureRule.legendre_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                case 5:
                    HermiteQuadrature.hermite_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                case 6:
                    HermiteQuadrature.gen_hermite_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                case 7:
                    Laguerre.QuadratureRule.laguerre_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                case 8:
                    Laguerre.QuadratureRule.gen_laguerre_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                case 9:
                    JacobiQuadrature.jacobi_compute_weights_np(
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                case 10:
                    HermiteQuadrature.hermite_genz_keister_lookup_weights_np(
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                case 11:
                case 12:
                    gw_compute_weights[dim](
                        order_1d[dim], np[dim], p.Skip(+p_index).ToArray(), weight_1d);
                    break;
                default:
                    Console.WriteLine("");
                    Console.WriteLine("PRODUCT_MIXED_GROWTH_WEIGHT - Fatal error!");
                    Console.WriteLine("  Unexpected value of RULE[" + dim + "] = "
                                      + rule[dim] + ".");
                    return;
            }

            p_index += np[dim];

            typeMethods.r8vec_direct_product2(ref data, dim, order_1d[dim], weight_1d,
                dim_num, order_nd, ref weight_nd);

        }
    }

    public static double[] product_weights(int dim_num, int[] order_1d, int order_nd, int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRODUCT_WEIGHTS computes the weights of a product rule.
        //
        //  Discussion:
        //
        //    This routine computes the weights for a quadrature rule which is
        //    a product of closed rules of varying order.
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
        //    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
        //
        //    Input, int ORDER_ND, the order of the product rule.
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
        //    Output, double PRODUCT_WEIGHTS_CC[DIM_NUM*ORDER_ND], 
        //    the product rule weights.
        //
    {
        int dim;
        int order;
        typeMethods.r8vecDPData data = new();

        double[] w_nd = new double[order_nd];

        for (order = 0; order < order_nd; order++)
        {
            w_nd[order] = 1.0;
        }

        for (dim = 0; dim < dim_num; dim++)
        {
            double[] w_1d;
            switch (rule)
            {
                case 1:
                    w_1d = ClenshawCurtis.cc_weights(order_1d[dim]);
                    break;
                case 2:
                    w_1d = Fejer1.f1_weights(order_1d[dim]);
                    break;
                case 3:
                    w_1d = Fejer2.f2_weights(order_1d[dim]);
                    break;
                case 4:
                    w_1d = PattersonQuadrature.gp_weights(order_1d[dim]);
                    break;
                case 5:
                    w_1d = GaussQuadrature.gl_weights(order_1d[dim]);
                    break;
                case 6:
                    w_1d = GaussHermite.gh_weights(order_1d[dim]);
                    break;
                case 7:
                    w_1d = Legendre.QuadratureRule.lg_weights(order_1d[dim]);
                    break;
                default:
                    Console.WriteLine("");
                    Console.WriteLine("PRODUCT_WEIGHTS - Fatal error!");
                    Console.WriteLine("  Unrecognized rule number = " + rule + "");
                    return null;
            }

            typeMethods.r8vec_direct_product2(ref data, dim, order_1d[dim], w_1d, dim_num,
                order_nd, ref w_nd);
        }

        return w_nd;
    }

    public static double[] product_weights_open(int dim_num, int[] order_1d, int order_nd,
            int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRODUCT_WEIGHTS_OPEN: weights for an open product rule.
        //
        //  Discussion:
        //
        //    This routine computes the weights for a quadrature rule which is
        //    a product of 1D rules of varying order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
        //
        //    Input, int ORDER_ND, the order of the product rule.
        //
        //    Input, int RULE, the 1D quadrature rule being used.
        //    2, Fejer Type 2 Rule;
        //    3, Gauss-Patterson Rule,
        //    4, Newton-Cotes Open Rule,
        //    5, Tanh-Sinh Rule.
        //
        //    Output, double PRODUCT_WEIGHTS_OPEN[DIM_NUM*ORDER_ND], the product 
        //    rule weights.
        //
    {
        int dim;
        int order;
        double[] w_1d = null;
        typeMethods.r8vecDPData data = new();

        double[] w_nd = new double[order_nd];

        for (order = 0; order < order_nd; order++)
        {
            w_nd[order] = 1.0;
        }

        for (dim = 0; dim < dim_num; dim++)
        {
            w_1d = rule switch
            {
                2 => Fejer2.f2_weights(order_1d[dim]),
                3 => PattersonQuadrature.gp_weights(order_1d[dim]),
                4 => NewtonCotesQuadrature.nco_weights(order_1d[dim]),
                5 => TanhSinh.ts_weights(order_1d[dim]),
                _ => w_1d
            };

            typeMethods.r8vec_direct_product2(ref data, dim, order_1d[dim], w_1d, dim_num,
                order_nd, ref w_nd);
        }

        return w_nd;
    }
}