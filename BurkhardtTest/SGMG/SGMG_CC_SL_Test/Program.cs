using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SGMG_CC_SL_Test;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SGMG_CC_SL.
        //
        //  Discussion:
        //
        //    SGMG_CC_SL tests SGMG with Clenshaw-Curtis rules with slow linear growth.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2013
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
        Console.WriteLine("SGMG_CC_SL");
        //
        //  Generate sparse grid rules and write them to files.
        //
        sgmg_cc_sl_tests();

        Console.WriteLine("");
        Console.WriteLine("SGMG_CC_SL");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void sgmg_cc_sl_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMG_CC_SL_TESTS calls SGMG_CC_SL_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num;
        string file_name;
        int[] growth;
        Func<int, int, double[], double[], double[]>[] gw_compute_points;
        Func<int, int, double[], double[], double[]>[] gw_compute_weights;
        int level_max;
        int[] np;
        int np_sum;
        double[] p;
        int[] rule;
        double tol;

        Console.WriteLine("");
        Console.WriteLine("SGMG_WRITE_TESTS");
        Console.WriteLine("  Call SGMG_WRITE_TEST with various arguments.");

        dim_num = 2;

        level_max = 2;

        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 1;

        growth = new int[dim_num];
        growth[0] = 1;
        growth[1] = 1;

        np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;

        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];

        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_points[1] = ClenshawCurtis.clenshaw_curtis_compute_points_np;

        gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
        gw_compute_weights[1] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;

        tol = 0.0000001;

        file_name = "sgmg_cc_sl";

        sgmg_create_rule(dim_num, level_max, rule, growth, np,
            p, gw_compute_points, gw_compute_weights, tol, file_name);

    }

    private static void sgmg_create_rule(int dim_num, int level_max,
            int[] rule, int[] growth, int[] np, double[] p,
            Func<int, int, double[], double[], double[]>[] gw_compute_points,
            Func<int, int, double[], double[], double[]>[] gw_compute_weights,
            double tol, string file_name)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    SGMG_CREATE_RULE creates the requested rule and writes it to a file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, integer DIM_NUM, the spatial dimension.
        //
        //    Input, integer LEVEL_MAX, the level that defines the grid.
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
        //    Input, string FILE_NAME, the main name of the output files.
        //
    {
        int point_num;
        int point_total_num;
        int[] sparse_index;
        int[] sparse_order;
        double[] sparse_point;
        int[] sparse_unique_index;
        double[] sparse_weight;

        Console.WriteLine("");
        Console.WriteLine("SGMG_CREATE_RULE");
        Console.WriteLine("  Create the requested sparse grid rule, and write it");
        Console.WriteLine("  to X, W and R files.");
        //
        //  Compute necessary data.
        //
        point_total_num = SGMG.sgmg_size_total(dim_num,
            level_max, rule, growth);

        point_num = SGMG.sgmg_size(dim_num, level_max,
            rule, np, p, gw_compute_points, tol, growth);

        sparse_unique_index = new int[point_total_num];

        SGMG.sgmg_unique_index(dim_num, level_max, rule,
            np, p, gw_compute_points, tol, point_num, point_total_num,
            growth, ref sparse_unique_index);

        sparse_order = new int[dim_num * point_num];
        sparse_index = new int[dim_num * point_num];

        SGMG.sgmg_index(dim_num, level_max, rule, point_num,
            point_total_num, sparse_unique_index, growth, ref sparse_order, ref sparse_index);
        //
        //  Compute points and weights.
        //
        sparse_point = new double [dim_num * point_num];

        SGMG.sgmg_point(dim_num, level_max, rule, np,
            p, gw_compute_points, point_num, sparse_order, sparse_index,
            growth, ref sparse_point);

        sparse_weight = new double[point_num];

        SGMG.sgmg_weight(dim_num, level_max, rule, np,
            p, gw_compute_weights, point_num, point_total_num, sparse_unique_index,
            growth, ref sparse_weight);
        //
        //  Write points and weights to files.
        //
        SGMG.sgmg_write(dim_num, rule, np, p,
            point_num, sparse_weight, sparse_point, file_name);

    }
}