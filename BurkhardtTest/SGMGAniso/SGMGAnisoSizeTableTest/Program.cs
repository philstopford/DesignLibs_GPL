using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Quadrature;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SGMGAnisoSizeTableTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SGMGA_SIZE_TABLE.
        //
        //  Discussion:
        //
        //    SGMGA_SIZE_TABLE tests the SGMGA_SIZE function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        DateTime ctime;
        int dim_max;
        int dim_min;
        int growth_1d;
        int level_max_max;
        int level_max_min;
        int np_1d;
        double[] p_1d;
        int rule_1d;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_SIZE_TABLE");
        Console.WriteLine("  Make tables of point counts.");
        Console.WriteLine("  Measure the CPU time for each table.");
        //
        //  Clenshaw-Curtis Grid (1), slow exponential growth, rule (4).
        //
        rule_1d = 1;
        growth_1d = 4;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, ClenshawCurtis.clenshaw_curtis_compute_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Clenshaw-Curtis Grid (1), exponential growth (6).
        //
        rule_1d = 1;
        growth_1d = 6;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, ClenshawCurtis.clenshaw_curtis_compute_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Clenshaw-Curtis Grid (1), exponential growth (6).
        //
        rule_1d = 1;
        growth_1d = 6;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, ClenshawCurtis.clenshaw_curtis_compute_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Clenshaw Curtis Grid (1), exponential growth (6)
        //
        rule_1d = 1;
        growth_1d = 6;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 100;
        dim_max = 100;
        level_max_min = 0;
        level_max_max = 2;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, ClenshawCurtis.clenshaw_curtis_compute_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss-Patterson Grid (3), slow exponential growth (4).
        //
        rule_1d = 3;
        growth_1d = 4;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, PattersonQuadrature.patterson_lookup_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss-Patterson Grid (3), moderate exponential growth (5).
        //
        rule_1d = 3;
        growth_1d = 5;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, PattersonQuadrature.patterson_lookup_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss-Patterson Grid (3), exponential growth (6).
        //
        rule_1d = 3;
        growth_1d = 6;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 6;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, PattersonQuadrature.patterson_lookup_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss Legendre Grid (4), slow linear odd growth (2)
        //
        rule_1d = 4;
        growth_1d = 2;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, Burkardt.Legendre.QuadratureRule.legendre_compute_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss Legendre Grid (4), moderate linear growth (3)
        //
        rule_1d = 4;
        growth_1d = 3;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, Burkardt.Legendre.QuadratureRule.legendre_compute_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss Laguerre Grid (7),  moderate linear growth (3).
        //
        rule_1d = 7;
        growth_1d = 3;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, Burkardt.Laguerre.QuadratureRule.laguerre_compute_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Hermite Genz Keister (10), slow exponential growth (4).
        //
        rule_1d = 10;
        growth_1d = 4;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmga_size_tabulate(rule_1d, growth_1d, np_1d, p_1d, dim_min, dim_max,
            level_max_min, level_max_max, HermiteQuadrature.hermite_genz_keister_lookup_points_np);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + (DateTime.Now - ctime).TotalSeconds + "");

        Console.WriteLine("");
        Console.WriteLine("SGMGA_SIZE_TABLE");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void sgmga_size_tabulate(int rule_1d, int growth_1d, int np_1d, double[] p_1d,
            int dim_min, int dim_max, int level_max_min, int level_max_max,
            Func<int, int, double[], double[], double[]> gw_compute_points_1d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_SIZE_TABULATE tests SGMGA_SIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int RULE_1D, the 1D rule.
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
        //    Input, int GROWTH_1D, the desired growth in each dimension.
        //    0, "DF", default growth associated with this quadrature rule;
        //    1, "SL", slow linear, L+1;
        //    2  "SO", slow linear odd, O=1+2((L+1)/2)
        //    3, "ML", moderate linear, 2L+1;
        //    4, "SE", slow exponential;
        //    5, "ME", moderate exponential;
        //    6, "FE", full exponential.
        //
        //    Input, int NP_1D, the number of parameters in the 1D rule.
        //
        //    Input, double P_1D[NP_1D], the parameters.
        //
        //    Input, int DIM_MIN, the minimum spatial dimension to consider.
        //
        //    Input, int DIM_MAX, the maximum spatial dimension to consider.
        //
        //    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
        //
        //    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
        //
        //    Input, GW_COMPUTE_POINTS_1D ( int order, int np, double p[], double x[] ),
        //    a function which return the 1D quadrature points.
        //
    {
        int dim;
        int dim_num;
        int[] growth;
        Func<int, int, double[], double[], double[]>[] gw_compute_points;
        int i;
        int level_max;
        double[] level_weight;
        int[] np;
        int np_sum;
        double[] p;
        int point_num;
        int[] rule;
        double tol;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_SIZE_TABULATE");
        Console.WriteLine("  SGMGA_SIZE returns the number of distinct");
        Console.WriteLine("  points in a sparse grid.");
        Console.WriteLine("");
        Console.WriteLine("  We use the same rule in all dimensions, and count the points,");
        Console.WriteLine("  for a range of dimensions and levels.");
        Console.WriteLine("");
        Console.WriteLine("  1D rule index = " + rule_1d + "");
        Console.WriteLine("  1D growth rule = " + growth_1d + "");
        Console.WriteLine("");

        tol = Math.Sqrt(typeMethods.r8_epsilon());

        string cout = "   DIM: ";
        for (dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString().PadLeft(8);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString().PadLeft(4);
            for (dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                level_weight = new double[dim_num];
                rule = new int[dim_num];
                growth = new int[dim_num];
                np = new int[dim_num];
                np_sum = dim_num * np_1d;
                p = new double[np_sum];
                gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];

                for (dim = 0; dim < dim_num; dim++)
                {
                    level_weight[dim] = 1.0;
                    rule[dim] = rule_1d;
                    growth[dim] = growth_1d;
                    np[dim] = np_1d;
                    for (i = 0; i < np_1d; i++)
                    {
                        p[i + dim * np_1d] = p_1d[i];
                    }

                    gw_compute_points[dim] = gw_compute_points_1d;
                }

                point_num = SGMGAniso.sgmga_size(dim_num, level_weight, level_max, rule,
                    np, p, gw_compute_points, tol, growth);

                cout += "  " + point_num.ToString().PadLeft(8);

            }

            Console.WriteLine(cout);
        }
    }
}