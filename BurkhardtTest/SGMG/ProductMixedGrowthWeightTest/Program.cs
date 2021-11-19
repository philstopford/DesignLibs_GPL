using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace ProductMixedGrowthWeightTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PRODUCT_MIXED_GROWTH_WEIGHT_TEST.
        //
        //  Discussion:
        //
        //    PRODUCT_MIXED_GROWTH_WEIGHT tests PRODUCT_MIXED_GROWTH_WEIGHT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 December 2009
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
    {
        Console.WriteLine("");
        Console.WriteLine("PRODUCT_MIXED_GROWTH_WEIGHT_TEST");
        //
        //  Make sure the individual product rule weights are computed correctly.
        //
        product_mixed_growth_weight_tests();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("PRODUCT_MIXED_GROWTH_WEIGHT");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");

    }

    private static void product_mixed_growth_weight_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRODUCT_MIXED_GROWTH_WEIGHT_TESTS calls PRODUCT_MIXED_GROWTH_WEIGHT_TEST.
        //
        //  Discussion:
        //
        //    To test Golub Welsch rules for a spatial dimension DIM, we can
        //    set RULE[DIM] = 10, and set the corresponding entry of 
        //    GW_COMPUTE_WEIGHTS to the name of a function that we know is already
        //    available, such as "webbur::clenshaw_curtis_compute_weights".
        //
        //    Note that, for ALL the tests, we set every entry of the GW_COMPUTE_WEIGHTS
        //    array.  However, a particular entry is only inspected if the corresponding
        //    entry of RULE is 10.
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
    {
        int dim_num;
        Func<int, int, double[], double[], double[]>[] gw_compute_weights;
        int[] np;
        int np_sum;
        int[] order_1d;
        int order_nd;
        double[] p;
        int[] rule;

        Console.WriteLine("");
        Console.WriteLine("PRODUCT_MIXED_GROWTH_WEIGHT_TESTS");
        Console.WriteLine("  Call PRODUCT_MIXED_GROWTH_WEIGHT_TEST with various arguments.");

        dim_num = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        order_1d = new int[dim_num];
        order_1d[0] = 3;
        order_1d[1] = 5;
        order_nd = typeMethods.i4vec_product(dim_num, order_1d);
        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 1;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        product_mixed_growth_weight_test(dim_num, order_1d, order_nd, rule,
            np, p, gw_compute_weights);

        dim_num = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        order_1d = new int[dim_num];
        rule = new int[dim_num];
        order_1d[0] = 3;
        order_1d[1] = 7;
        order_nd = typeMethods.i4vec_product(dim_num, order_1d);
        rule[0] = 1;
        rule[1] = 5;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = HermiteQuadrature.hermite_compute_weights_np;
        product_mixed_growth_weight_test(dim_num, order_1d, order_nd, rule,
            np, p, gw_compute_weights);

        dim_num = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        order_1d = new int[dim_num];
        rule = new int[dim_num];
        order_1d[0] = 3;
        order_1d[1] = 3;
        order_nd = typeMethods.i4vec_product(dim_num, order_1d);
        rule[0] = 3;
        rule[1] = 7;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = PattersonQuadrature.patterson_lookup_weights_np;
        gw_compute_weights[1] = Burkardt.Laguerre.QuadratureRule.laguerre_compute_weights_np;
        product_mixed_growth_weight_test(dim_num, order_1d, order_nd, rule,
            np, p, gw_compute_weights);

        dim_num = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 1;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        p[0] = 1.5;
        order_1d = new int[dim_num];
        rule = new int[dim_num];
        order_1d[0] = 5;
        order_1d[1] = 5;
        order_nd = typeMethods.i4vec_product(dim_num, order_1d);
        rule[0] = 1;
        rule[1] = 8;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = Burkardt.Laguerre.QuadratureRule.gen_laguerre_compute_weights_np;
        product_mixed_growth_weight_test(dim_num, order_1d, order_nd, rule,
            np, p, gw_compute_weights);

        dim_num = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        p[0] = 0.5;
        p[1] = 1.5;
        order_1d = new int[dim_num];
        rule = new int[dim_num];
        order_1d[0] = 5;
        order_1d[1] = 5;
        order_nd = typeMethods.i4vec_product(dim_num, order_1d);
        rule[0] = 2;
        rule[1] = 9;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = Fejer2.fejer2_compute_weights_np;
        gw_compute_weights[1] = JacobiQuadrature.jacobi_compute_weights_np;
        product_mixed_growth_weight_test(dim_num, order_1d, order_nd, rule,
            np, p, gw_compute_weights);

        dim_num = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        p[0] = 2.0;
        order_1d = new int[dim_num];
        rule = new int[dim_num];
        order_1d[0] = 7;
        order_1d[1] = 9;
        order_nd = typeMethods.i4vec_product(dim_num, order_1d);
        rule[0] = 6;
        rule[1] = 10;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = HermiteQuadrature.gen_hermite_compute_weights_np;
        gw_compute_weights[1] = HermiteQuadrature.hermite_genz_keister_lookup_weights_np;
        product_mixed_growth_weight_test(dim_num, order_1d, order_nd, rule,
            np, p, gw_compute_weights);

        dim_num = 3;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np[2] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        order_1d = new int[dim_num];
        rule = new int[dim_num];
        order_1d[0] = 2;
        order_1d[1] = 3;
        order_1d[2] = 3;
        order_nd = typeMethods.i4vec_product(dim_num, order_1d);
        rule[0] = 1;
        rule[1] = 4;
        rule[2] = 5;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_weights_np;
        gw_compute_weights[2] = HermiteQuadrature.hermite_compute_weights_np;
        product_mixed_growth_weight_test(dim_num, order_1d, order_nd, rule,
            np, p, gw_compute_weights);
        //
        //  Repeat, treating  rules #2 and #3 as Golub Welsch rules.
        //
        dim_num = 3;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np[2] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        order_1d = new int[dim_num];
        rule = new int[dim_num];
        order_1d[0] = 2;
        order_1d[1] = 3;
        order_1d[2] = 3;
        order_nd = typeMethods.i4vec_product(dim_num, order_1d);
        rule[0] = 1;
        rule[1] = 11;
        rule[2] = 11;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_weights_np;
        gw_compute_weights[2] = HermiteQuadrature.hermite_compute_weights_np;
        product_mixed_growth_weight_test(dim_num, order_1d, order_nd, rule,
            np, p, gw_compute_weights);
        //
        //  Dimension 2, Rule 3
        //
        dim_num = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        order_1d = new int[dim_num];
        rule = new int[dim_num];
        order_1d[0] = 15;
        order_1d[1] = 15;
        order_nd = typeMethods.i4vec_product(dim_num, order_1d);
        rule[0] = 3;
        rule[1] = 3;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = PattersonQuadrature.patterson_lookup_weights_np;
        gw_compute_weights[1] = PattersonQuadrature.patterson_lookup_weights_np;
        product_mixed_growth_weight_test(dim_num, order_1d, order_nd, rule,
            np, p, gw_compute_weights);
    }

    private static void product_mixed_growth_weight_test(int dim_num, int[] order_1d,
            int order_nd, int[] rule, int[] np, double[] p,
            Func<int, int, double[], double[], double[]>[] gw_compute_weights)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    PRODUCT_MIXED_GROWTH_WEIGHT_TEST: weights of a mixed factor product rule.
        //
        //  Discussion:
        //
        //    This routine computes a sparse grid and compares the sum of the weights
        //    to the expected exact value.
        //
        //    The routine cannot produce a result for rules that include one or more
        //    component rules of type 10, that is, Golub-Welsch rules.
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
        //    Input, int NP[RULE_NUM], the number of parameters used by each rule.
        //
        //    Input, double P[sum(NP[*])], the parameters needed by each rule.
        //
        //    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
        //    an array of pointers to functions which return the 1D quadrature weights 
        //    associated with each spatial dimension for which a Golub Welsch rule 
        //    is used.
        //
    {
        double alpha;
        double beta;
        double arg1;
        double arg2;
        double arg3;
        double arg4;
        int dim;
        int i;
        int p_index;
        double value1;
        double value2;
        double[] weight;
        double weight_sum;
        double weight_sum_error;
        double weight_sum_exact;
        //
        //  Determine the integral of 1 over the multidimensional weighted region.
        //
        p_index = 0;

        weight_sum_exact = 1.0;

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (rule[dim])
            {
                case 1:
                case 2:
                case 3:
                case 4:
                    weight_sum_exact *= 2.0;
                    break;
                case 5:
                    weight_sum_exact *= Math.Sqrt(Math.PI);
                    break;
                case 6:
                    alpha = p[p_index];
                    p_index += 1;

                    weight_sum_exact *= typeMethods.r8_gamma(0.5 * (alpha + 1.0));
                    break;
                case 7:
                    weight_sum_exact *= 1.0;
                    break;
                case 8:
                    alpha = p[p_index];
                    p_index += 1;

                    weight_sum_exact *= typeMethods.r8_gamma(alpha + 1.0);
                    break;
                case 9:
                    alpha = p[p_index];
                    p_index += 1;
                    beta = p[p_index];
                    p_index += 1;
                    arg1 = -alpha;
                    arg2 = 1.0;
                    arg3 = beta + 2.0;
                    arg4 = -1.0;
                    value1 = typeMethods.r8_hyper_2f1(arg1, arg2, arg3, arg4);
                    arg1 = -beta;
                    arg2 = 1.0;
                    arg3 = alpha + 2.0;
                    arg4 = -1.0;
                    value2 = typeMethods.r8_hyper_2f1(arg1, arg2, arg3, arg4);
                    weight_sum_exact *= (
                        value1 / (beta + 1.0) + value2 / (alpha + 1.0));
                    break;
                case 10:
                    weight_sum_exact *= Math.Sqrt(Math.PI);
                    break;
                case 11:
                {
                    for (i = 0; i < np[dim]; i++)
                    {
                        alpha = p[p_index];
                        p_index += 1;
                    }

                    weight_sum_exact = 0.0;
                    break;
                }
                case 12:
                {
                    for (i = 0; i < np[dim]; i++)
                    {
                        alpha = p[p_index];
                        p_index += 1;
                    }

                    weight_sum_exact = 0.0;
                    break;
                }
                default:
                    Console.WriteLine("");
                    Console.WriteLine("PRODUCT_MIXED_GROWTH_WEIGHT_TEST - Fatal error!");
                    Console.WriteLine("  Unexpected value of RULE[" + dim + "] = "
                                      + rule[dim] + ".");
                    return;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("PRODUCT_MIXED_GROWTH_WEIGHT_TEST:");
        Console.WriteLine("  Compute the weights of a mixed factor product grid.");
        Console.WriteLine("");
        Console.WriteLine("  As a simple test, sum these weights.");
        Console.WriteLine("  They should sum to exactly " + weight_sum_exact + "");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");

        Console.WriteLine("");
        Console.WriteLine(" Dimension      Rule     Order        Parameters");
        Console.WriteLine("");

        p_index = 0;
        for (dim = 0; dim < dim_num; dim++)
        {
            switch (rule[dim])
            {
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + order_1d[dim].ToString().PadLeft(8) + "");
                    break;
                case 6:
                    alpha = p[p_index];
                    p_index += 1;
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + order_1d[dim].ToString().PadLeft(8)
                                           + "  " + alpha.ToString().PadLeft(14) + "");
                    break;
                case 7:
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + order_1d[dim].ToString().PadLeft(8) + "");
                    break;
                case 8:
                    alpha = p[p_index];
                    p_index += 1;
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + order_1d[dim].ToString().PadLeft(8)
                                           + "  " + alpha.ToString().PadLeft(14) + "");
                    break;
                case 9:
                    alpha = p[p_index];
                    p_index += 1;
                    beta = p[p_index];
                    p_index += 1;
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + order_1d[dim].ToString().PadLeft(8)
                                           + "  " + alpha.ToString().PadLeft(14)
                                           + "  " + beta.ToString().PadLeft(14) + "");
                    break;
                case 10:
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + order_1d[dim].ToString().PadLeft(8) + "");
                    break;
                case 11:
                {
                    string cout = "  " + dim.ToString().PadLeft(8)
                                       + "  " + rule[dim].ToString().PadLeft(8)
                                       + "  " + order_1d[dim].ToString().PadLeft(8);
                    for (i = 0; i < np[dim]; i++)
                    {
                        alpha = p[p_index];
                        p_index += 1;
                        cout += "  " + alpha.ToString().PadLeft(14);
                    }

                    Console.WriteLine(cout);
                    break;
                }
                case 12:
                {
                    string cout = "  " + dim.ToString().PadLeft(8)
                                       + "  " + rule[dim].ToString().PadLeft(8)
                                       + "  " + order_1d[dim].ToString().PadLeft(8);
                    for (i = 0; i < np[dim]; i++)
                    {
                        alpha = p[p_index];
                        p_index += 1;
                        cout += "  " + alpha.ToString().PadLeft(14);
                    }

                    Console.WriteLine(cout);
                    break;
                }
                default:
                    Console.WriteLine("");
                    Console.WriteLine("PRODUCT_MIXED_GROWTH_WEIGHT_TEST - Fatal error!");
                    Console.WriteLine("  Cannot perform test for rule = " + rule[dim] + "");
                    return;
            }
        }

        //
        //  Compute the weights and points of the sparse grid mixed rule.
        //
        weight = new double[order_nd];

        Product.product_mixed_growth_weight(dim_num, order_1d, order_nd, rule,
            np, p, gw_compute_weights, ref weight);
        //
        //  Sum the weights to get the SGM approximation to the integral of 1.
        //
        weight_sum = typeMethods.r8vec_sum(order_nd, weight);
        //
        //  Compare the exact and estimated integrals.
        //
        weight_sum_error = typeMethods.r8_abs(weight_sum - weight_sum_exact);

        Console.WriteLine("");
        Console.WriteLine("    Weight sum  Expected sum    Difference");
        Console.WriteLine("");
        Console.WriteLine("  " + weight_sum.ToString().PadLeft(14)
                               + "  " + weight_sum_exact.ToString().PadLeft(14)
                               + "  " + weight_sum_error.ToString().PadLeft(14) + "");

    }
}