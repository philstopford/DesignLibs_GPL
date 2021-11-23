using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Quadrature;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SGMGWeightTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SGMG_WEIGHT_TEST.
        //
        //  Discussion:
        //
        //    SGMG_WEIGHT_TEST tests SGMG_WEIGHT.
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
        double tol = typeMethods.r8_epsilon();

        Console.WriteLine("");
        Console.WriteLine("SGMG_WEIGHT_TEST");
        //
        //  Generate the weights for sparse grid rules.
        //
        sgmg_weight_tests(tol);
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("SGMG_WEIGHT_TEST");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");

    }

    private static void sgmg_weight_tests(double tol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMG_WEIGHT_TESTS calls SGMG_WEIGHT_TEST.
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
        //    Input, double TOL, a tolerance for point equality.
        //    A value of sqrt ( eps ) is reasonable, and will allow the code to
        //    consolidate points which are equal, or very nearly so.  A value of
        //    -1.0, on the other hand, will force the code to use every point, 
        //    regardless of duplication.
        //
    {
        int dim_num;
        int[] growth;
        Func<int, int, double[], double[], double[]>[] gw_compute_points;
        Func<int, int, double[], double[], double[]>[] gw_compute_weights;
        int level_max_max;
        int level_max_min;
        int[] np;
        int np_sum;
        double[] p;
        int[] rule;

        Console.WriteLine("");
        Console.WriteLine("SGMG_WEIGHT_TESTS");
        Console.WriteLine("  Call SGMG_WEIGHT_TEST with various arguments.");
        Console.WriteLine("");
        Console.WriteLine("  All tests will use a point equality tolerance of " + tol + "");

        dim_num = 2;
        level_max_min = 0;
        level_max_max = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 1;
        growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 6;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_points[1] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);

        dim_num = 2;
        level_max_min = 0;
        level_max_max = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 3;
        growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 6;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_points[1] = Fejer2.fejer2_compute_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = Fejer2.fejer2_compute_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);

        dim_num = 2;
        level_max_min = 0;
        level_max_max = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 4;
        growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 3;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_points[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);

        dim_num = 2;
        level_max_min = 0;
        level_max_max = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 7;
        growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 3;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_points[1] = Burkardt.Laguerre.QuadratureRule.laguerre_compute_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = Burkardt.Laguerre.QuadratureRule.laguerre_compute_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);

        dim_num = 2;
        level_max_min = 0;
        level_max_max = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 1;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        p[0] = 1.5;
        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 8;
        growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 3;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_points[1] = Burkardt.Laguerre.QuadratureRule.gen_laguerre_compute_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = Burkardt.Laguerre.QuadratureRule.gen_laguerre_compute_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);

        dim_num = 2;
        level_max_min = 0;
        level_max_max = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 2;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        p[0] = 0.5;
        p[1] = 1.5;
        rule = new int[dim_num];
        rule[0] = 2;
        rule[1] = 9;
        growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 3;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = Fejer2.fejer2_compute_points_np;
        gw_compute_points[1] = JacobiQuadrature.jacobi_compute_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = Fejer2.fejer2_compute_weights_np;
        gw_compute_weights[1] = JacobiQuadrature.jacobi_compute_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);

        dim_num = 2;
        level_max_min = 0;
        level_max_max = 2;
        np = new int[dim_num];
        np[0] = 1;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        p[0] = 2.0;
        rule = new int[dim_num];
        rule[0] = 6;
        rule[1] = 10;
        growth = new int[dim_num];
        growth[0] = 3;
        growth[1] = 4;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = HermiteQuadrature.gen_hermite_compute_points_np;
        gw_compute_points[1] = HermiteQuadrature.hermite_genz_keister_lookup_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = HermiteQuadrature.gen_hermite_compute_weights_np;
        gw_compute_weights[1] = HermiteQuadrature.hermite_genz_keister_lookup_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);

        dim_num = 3;
        level_max_min = 0;
        level_max_max = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np[2] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 4;
        rule[2] = 5;
        growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 3;
        growth[2] = 3;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_points[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_points_np;
        gw_compute_points[2] = HermiteQuadrature.hermite_compute_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_weights_np;
        gw_compute_weights[2] = HermiteQuadrature.hermite_compute_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);
        //
        //  Repeat, treating  rules #2 and #3 as Golub Welsch rules.
        //
        dim_num = 3;
        level_max_min = 0;
        level_max_max = 2;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np[2] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 11;
        rule[2] = 11;
        growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 3;
        growth[2] = 3;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_points[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_points_np;
        gw_compute_points[2] = HermiteQuadrature.hermite_compute_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_weights_np;
        gw_compute_weights[2] = HermiteQuadrature.hermite_compute_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);
        //
        //  Dimension 2, Level 4, Rule 3, exponential growth.
        //
        dim_num = 2;
        level_max_min = 0;
        level_max_max = 4;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        rule = new int[dim_num];
        rule[0] = 3;
        rule[1] = 3;
        growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 6;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = PattersonQuadrature.patterson_lookup_points_np;
        gw_compute_points[1] = PattersonQuadrature.patterson_lookup_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = PattersonQuadrature.patterson_lookup_weights_np;
        gw_compute_weights[1] = PattersonQuadrature.patterson_lookup_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);
        //
        //  Dimension 2, Level 4, Rule 3, slow exponential growth.
        //
        dim_num = 2;
        level_max_min = 0;
        level_max_max = 4;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        rule = new int[dim_num];
        rule[0] = 3;
        rule[1] = 3;
        growth = new int[dim_num];
        growth[0] = 4;
        growth[1] = 4;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = PattersonQuadrature.patterson_lookup_points_np;
        gw_compute_points[1] = PattersonQuadrature.patterson_lookup_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = PattersonQuadrature.patterson_lookup_weights_np;
        gw_compute_weights[1] = PattersonQuadrature.patterson_lookup_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);
        //
        //  Dimension 2, Level 4, Rule 3, moderate exponential growth.
        //
        dim_num = 2;
        level_max_min = 0;
        level_max_max = 4;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        rule = new int[dim_num];
        rule[0] = 3;
        rule[1] = 3;
        growth = new int[dim_num];
        growth[0] = 5;
        growth[1] = 5;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = PattersonQuadrature.patterson_lookup_points_np;
        gw_compute_points[1] = PattersonQuadrature.patterson_lookup_points_np;
        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = PattersonQuadrature.patterson_lookup_weights_np;
        gw_compute_weights[1] = PattersonQuadrature.patterson_lookup_weights_np;
        sgmg_weight_test(dim_num, level_max_min, level_max_max,
            rule, growth, np, p, gw_compute_points, gw_compute_weights, tol);

    }

    private static void sgmg_weight_test(int dim_num, int level_max_min,
            int level_max_max, int[] rule, int[] growth, int[] np, double[] p,
            Func<int, int, double[], double[], double[]>[] gw_compute_points,
            Func<int, int, double[], double[], double[]>[] gw_compute_weights,
            double tol)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    SGMG_WEIGHT_TEST checks the sum of the quadrature weights.
        //
        //  Discussion:
        //
        //    If any component rule is of Golub-Welsch type, we cannot compute
        //    the exact weight sum, which we set, instead, to zero.
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
        //    Input, int LEVEL_MAX_MIN, LEVEL_MAX_MAX, the minimum and
        //    maximum values of LEVEL_MAX.
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
        //    Input, int GROWTH[DIM_NUM], the growth rule in each dimension. 
        //    0, "DF", default growth associated with this quadrature rule;
        //    1, "SL", slow linear, L+1;
        //    2  "SO", slow linear odd, O=1+2((L+1)/2)
        //    3, "ML", moderate linear, 2L+1;
        //    4, "SE", slow exponential;
        //    5, "ME", moderate exponential;
        //    6, "FE", full exponential.
        //
        //    Input, int NP[RULE_NUM], the number of parameters used by each rule.
        //
        //    Input, double P[sum(NP[*])], the parameters needed by each rule.
        //
        //    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double x[] ),
        //    an array of pointers to functions which return the 1D quadrature points 
        //    associated with each spatial dimension for which a Golub Welsch rule 
        //    is used.
        //
        //    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
        //    an array of pointers to functions which return the 1D quadrature weights 
        //    associated with each spatial dimension for which a Golub Welsch rule 
        //    is used.
        //
        //    Input, double TOL, a tolerance for point equality.
        //
    {
        double alpha;
        double arg1;
        double arg2;
        double arg3;
        double arg4;
        double beta;
        int dim;
        int i;
        int level_max;
        int p_index;
        int point_num;
        int point_total_num;
        int[] sparse_unique_index;
        double[] sparse_weight;
        double value1;
        double value2;
        double weight_sum;
        double weight_sum_error;
        double weight_sum_exact;

        Console.WriteLine("");
        Console.WriteLine("SGMG_WEIGHT_TEST");
        Console.WriteLine("  Compute the weights of a sparse grid.");
        Console.WriteLine("");
        Console.WriteLine("  Each sparse grid is of spatial dimension DIM_NUM,");
        Console.WriteLine("  and is made up of product grids of levels up to LEVEL_MAX.");
        Console.WriteLine("");
        Console.WriteLine(" Dimension      Rule  Growth rate       Parameters");
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
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                    break;
                case 6:
                    alpha = p[p_index];
                    p_index += 1;
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8)
                                           + "  " + alpha.ToString().PadLeft(14) + "");
                    break;
                case 7:
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                    break;
                case 8:
                    alpha = p[p_index];
                    p_index += 1;
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8)
                                           + "  " + alpha.ToString().PadLeft(14) + "");
                    break;
                case 9:
                    alpha = p[p_index];
                    p_index += 1;
                    beta = p[p_index];
                    p_index += 1;
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8)
                                           + "  " + alpha.ToString().PadLeft(14)
                                           + "  " + beta.ToString().PadLeft(14) + "");
                    break;
                case 10:
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                    break;
                case 11:
                {
                    string cout = "  " + dim.ToString().PadLeft(8)
                                       + "  " + rule[dim].ToString().PadLeft(8)
                                       + "  " + growth[dim].ToString().PadLeft(8);
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
                                       + "  " + growth[dim].ToString().PadLeft(8);
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
                    Console.WriteLine("SPARSE_GRID_MIXED_WEIGHT_TEST - Fatal error!");
                    Console.WriteLine("  Unexpected value of RULE = " + rule[dim] + "");
                    return;
            }
        }

        weight_sum_exact = 1.0;

        p_index = 0;
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
                    weight_sum_exact *= value1 / (beta + 1.0) + value2 / (alpha + 1.0);
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
                    Console.WriteLine("SPARSE_GRID_MIXED_WEIGHT_TEST - Fatal error!");
                    Console.WriteLine("  Unexpected value of RULE[" + dim + "] = "
                                      + rule[dim] + ".");
                    return;
            }
        }

        switch (weight_sum_exact)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("  Because this rule includes Golub-Welsch components,");
                Console.WriteLine("  we do not try to compute the exact weight sum.");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  As a simple test, sum these weights.");
                Console.WriteLine("  They should sum to exactly " + weight_sum_exact + "");
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("     Level      Weight sum  Expected sum    Difference");
        Console.WriteLine("");

        for (level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            point_total_num = SGMG.sgmg_size_total(dim_num,
                level_max, rule, growth);

            point_num = SGMG.sgmg_size(dim_num, level_max,
                rule, np, p, gw_compute_points, tol, growth);

            sparse_unique_index = new int[point_total_num];

            SGMG.sgmg_unique_index(dim_num, level_max, rule,
                np, p, gw_compute_points, tol, point_num, point_total_num,
                growth, ref sparse_unique_index);

            sparse_weight = new double[point_num];

            SGMG.sgmg_weight(dim_num, level_max, rule, np,
                p, gw_compute_weights, point_num, point_total_num, sparse_unique_index,
                growth, ref sparse_weight);

            weight_sum = typeMethods.r8vec_sum(point_num, sparse_weight);

            weight_sum_error = typeMethods.r8_abs(weight_sum - weight_sum_exact);

            Console.WriteLine("  " + level_max.ToString().PadLeft(8)
                                   + "  " + weight_sum.ToString().PadLeft(14)
                                   + "  " + weight_sum_exact.ToString().PadLeft(14)
                                   + "  " + weight_sum_error.ToString().PadLeft(14) + "");

        }
    }
}