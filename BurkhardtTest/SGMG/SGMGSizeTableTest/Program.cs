using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Quadrature;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SGMGSizeTableTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SGMG_SIZE_TABLE.
        //
        //  Discussion:
        //
        //    SGMG_SIZE_TABLE_PRB makes a point count table.
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
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
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
        Console.WriteLine("SGMG_SIZE_TABLE");
        //
        //  Clenshaw-Curtis (1), slow exponential growth (4).
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
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, ClenshawCurtis.clenshaw_curtis_compute_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Clenshaw-Curtis (1), slow exponential growth (4).
        //
        rule_1d = 1;
        growth_1d = 4;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 5;
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, ClenshawCurtis.clenshaw_curtis_compute_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
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
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, ClenshawCurtis.clenshaw_curtis_compute_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
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
        level_max_max = 5;
        ctime = DateTime.Now;
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, ClenshawCurtis.clenshaw_curtis_compute_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Clenshaw-Curtis Grid (1), exponential growth (6).
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
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, ClenshawCurtis.clenshaw_curtis_compute_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss-Patterson (3), slow exponential growth (4).
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
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, PattersonQuadrature.patterson_lookup_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss-Patterson (3), slow exponential growth (4).
        //
        rule_1d = 3;
        growth_1d = 4;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 5;
        ctime = DateTime.Now;
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, PattersonQuadrature.patterson_lookup_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss-Patterson (3), Moderate Exponential Growth (5).
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
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, PattersonQuadrature.patterson_lookup_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss-Patterson (3), Moderate Exponential Growth (5).
        //
        rule_1d = 3;
        growth_1d = 5;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 5;
        ctime = DateTime.Now;
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, PattersonQuadrature.patterson_lookup_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss-Patterson (3), exponential growth (6).
        //
        rule_1d = 3;
        growth_1d = 6;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, PattersonQuadrature.patterson_lookup_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss-Patterson (3), exponential growth (6).
        //
        rule_1d = 3;
        growth_1d = 6;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 5;
        ctime = DateTime.Now;
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, PattersonQuadrature.patterson_lookup_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU Time = " + (DateTime.Now - ctime).TotalSeconds + "");
        //
        //  Gauss Legendre Grid (4), slow linear growth (1)
        //
        rule_1d = 4;
        growth_1d = 1;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 1;
        dim_max = 5;
        level_max_min = 0;
        level_max_max = 10;
        ctime = DateTime.Now;
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, Burkardt.Legendre.QuadratureRule.legendre_compute_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + ctime + "");
        //
        //  Gauss Legendre Grid (4), slow linear growth (1)
        //
        rule_1d = 4;
        growth_1d = 1;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 7;
        ctime = DateTime.Now;
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, Burkardt.Legendre.QuadratureRule.legendre_compute_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + ctime + "");
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
        level_max_max = 10;
        ctime = DateTime.Now;
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, Burkardt.Legendre.QuadratureRule.legendre_compute_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + ctime + "");
        //
        //  Gauss Legendre Grid (4), slow linear odd growth (2)
        //
        rule_1d = 4;
        growth_1d = 2;
        np_1d = 0;
        p_1d = new double[np_1d];
        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 8;
        ctime = DateTime.Now;
        sgmg_size_tabulate(rule_1d, growth_1d, np_1d,
            p_1d, Burkardt.Legendre.QuadratureRule.legendre_compute_points_np,
            dim_min, dim_max, level_max_min, level_max_max);
        Console.WriteLine("");
        Console.WriteLine("  CPU_TIME = " + ctime + "");
        Console.WriteLine("");
        Console.WriteLine("SGMG_SIZE_TABLE");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void sgmg_size_tabulate(int rule_1d, int growth_1d,
            int np_1d, double[] p_1d,
            Func<int, int, double[], double[], double[]> gw_compute_points_1d,
            int dim_min, int dim_max, int level_max_min, int level_max_max)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    SGMG_SIZE_TABULATE tests SGMG_SIZE.
        //
        //  Discussion:
        //
        //    We do NOT consider mixed rules.  Instead, we are looking at sparse grid
        //    rules for which all dimensions use the same 1D rule family.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 April 2014
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
        //    Input, int GROWTH_1D, the 1D growth rule. 
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
        //    Input, double P_1D[NP_1D], the parameters for the 1D rule.
        //
        //    Input, void gw_compute_points_1d ( int order, int np, double p[], double w[] ),
        //    the function to be used to compute points for the 1D rule.
        //
        //    Input, int DIM_MIN, the minimum spatial dimension.
        //
        //    Input, int DIM_MAX, the maximum spatial dimension.
        //
        //    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX.
        //
        //    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX.
        //
    {
        int dim;
        int dim_num;
        int[] growth;
        Func<int, int, double[], double[], double[]>[] gw_compute_points;
        int i;
        int j;
        int level_max;
        int[] np;
        int np_sum;
        double[] p;
        int point_num;
        int[] rule;
        double tol;

        Console.WriteLine("");
        Console.WriteLine("SGMG_SIZE_TABULATE");
        Console.WriteLine("  SGMG_SIZE returns the number of distinct");
        Console.WriteLine("  points in a sparse grid.");
        Console.WriteLine("");
        Console.WriteLine("  We use the same rule in all dimensions, and count the points");
        Console.WriteLine("  for a range of dimensions and levels.");
        Console.WriteLine("");
        Console.WriteLine("  1D rule index is " + rule_1d + "");
        Console.WriteLine("  1D growth index is " + growth_1d + "");
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
                rule = new int[dim_num];
                growth = new int[dim_num];
                np = new int[dim_num];
                np_sum = dim_num * np_1d;
                p = new double[np_sum];
                gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];

                j = 0;
                for (dim = 0; dim < dim_num; dim++)
                {
                    rule[dim] = rule_1d;
                    growth[dim] = growth_1d;
                    np[dim] = np_1d;
                    for (i = 0; i < np_sum; i++)
                    {
                        p[j] = p_1d[i];
                    }

                    gw_compute_points[dim] = gw_compute_points_1d;
                }

                point_num = SGMG.sgmg_size(dim_num, level_max, rule,
                    np, p, gw_compute_points, tol, growth);

                cout += "  " + point_num.ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
        }
    }
}