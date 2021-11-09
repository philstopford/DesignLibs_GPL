using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Types;

namespace Burkardt.Quadrature
{
    public static class Product
    {
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
            double[] w_1d = null;
            double[] w_nd = null;
            typeMethods.r8vecDPData data = new typeMethods.r8vecDPData();

            w_nd = new double[order_nd];

            for (order = 0; order < order_nd; order++)
            {
                w_nd[order] = 1.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                if (rule == 1)
                {
                    w_1d = ClenshawCurtis.cc_weights(order_1d[dim]);
                }
                else if (rule == 2)
                {
                    w_1d = Fejer1.f1_weights(order_1d[dim]);
                }
                else if (rule == 3)
                {
                    w_1d = Fejer2.f2_weights(order_1d[dim]);
                }
                else if (rule == 4)
                {
                    w_1d = PattersonQuadrature.gp_weights(order_1d[dim]);
                }
                else if (rule == 5)
                {
                    w_1d = GaussQuadrature.gl_weights(order_1d[dim]);
                }
                else if (rule == 6)
                {
                    w_1d = GaussHermite.gh_weights(order_1d[dim]);
                }
                else if (rule == 7)
                {
                    w_1d = Legendre.QuadratureRule.lg_weights(order_1d[dim]);
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("PRODUCT_WEIGHTS - Fatal error!");
                    Console.WriteLine("  Unrecognized rule number = " + rule + "");
                    return (null);
                }

                typeMethods.r8vec_direct_product2(ref data, dim, order_1d[dim], w_1d, dim_num,
                    order_nd, ref w_nd);

            }

            return w_nd;
        }
        
        public static double[] product_weights_open ( int dim_num, int[] order_1d, int order_nd, 
        int rule )

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
            double[] w_nd;
            typeMethods.r8vecDPData data = new typeMethods.r8vecDPData();

            w_nd = new double[order_nd];

            for ( order = 0; order < order_nd; order++ )
            {
                w_nd[order] = 1.0;
            }

            for ( dim = 0; dim < dim_num; dim++ )
            {
                if ( rule == 2 )
                {
                    w_1d = Fejer2.f2_weights ( order_1d[dim] );
                }
                else if ( rule == 3 )
                {
                    w_1d = PattersonQuadrature.gp_weights ( order_1d[dim] );
                }
                else if ( rule == 4 )
                {
                    w_1d = NewtonCotesQuadrature.nco_weights ( order_1d[dim] );
                }
                else if ( rule == 5 )
                {
                    w_1d = TanhSinh.ts_weights ( order_1d[dim] );
                }

                typeMethods.r8vec_direct_product2 ( ref data, dim, order_1d[dim], w_1d, dim_num, 
                    order_nd, ref w_nd );

            }

            return w_nd;
        }
    }
}