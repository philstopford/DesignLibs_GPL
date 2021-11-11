using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Quadrature;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SGMGSizeTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SGMG_SIZE_TEST.
            //
            //  Discussion:
            //
            //    SGMG_TEST tests the SGMG routines.
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
            double tol;

            Console.WriteLine("");
            Console.WriteLine("SGMG_SIZE_TEST");
            //
            //  1) Using a tolerance that is less than 0 means that there will be no
            //  consolidation of duplicate points.
            //
            //  2) Using a small positive tolerance means there will be consolidation of
            //  points whose maximum difference is no more than TOL.
            //
            tol = -1.0;
            sgmg_size_tests(tol);

            tol = Math.Sqrt(typeMethods.r8_epsilon());
            sgmg_size_tests(tol);
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("SGMG_SIZE_TEST");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");

        }

        static void sgmg_size_tests(double tol)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SGMG_SIZE_TESTS calls SGMG_SIZE_TEST.
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
            int level_max_max;
            int level_max_min;
            int[] np;
            int np_sum;
            int[] order_1d;
            int order_nd;
            double[] p;
            int[] rule;

            Console.WriteLine("");
            Console.WriteLine("SGMG_SIZE_TESTS");
            Console.WriteLine("  Call SGMG_SIZE_TEST with various arguments.");
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
            sgmg_size_test(dim_num, level_max_min, level_max_max,
                rule, growth, np, p, gw_compute_points, tol);

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
            gw_compute_points[1] = PattersonQuadrature.patterson_lookup_points_np;
            sgmg_size_test(dim_num, level_max_min, level_max_max,
                rule, growth, np, p, gw_compute_points, tol);

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
            sgmg_size_test(dim_num, level_max_min, level_max_max,
                rule, growth, np, p, gw_compute_points, tol);

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
            sgmg_size_test(dim_num, level_max_min, level_max_max,
                rule, growth, np, p, gw_compute_points, tol);

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
            sgmg_size_test(dim_num, level_max_min, level_max_max,
                rule, growth, np, p, gw_compute_points, tol);

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
            sgmg_size_test(dim_num, level_max_min, level_max_max,
                rule, growth, np, p, gw_compute_points, tol);

            dim_num = 2;
            level_max_min = 0;
            level_max_max = 2;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 1;
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
            sgmg_size_test(dim_num, level_max_min, level_max_max,
                rule, growth, np, p, gw_compute_points, tol);

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
            sgmg_size_test(dim_num, level_max_min, level_max_max,
                rule, growth, np, p, gw_compute_points, tol);
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
            sgmg_size_test(dim_num, level_max_min, level_max_max,
                rule, growth, np, p, gw_compute_points, tol);

            return;
        }

        static void sgmg_size_test(int dim_num, int level_max_min,
                int level_max_max, int[] rule, int[] growth, int[] np, double[] p,
                Func<int, int, double[], double[], double[]>[] gw_compute_points,
                double tol)

            //***************************************************************************80
            //
            //  Purpose:
            //
            //    SGMG_SIZE_TEST tests SGMG_SIZE, SGMG_SIZE_TOTAL.
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
            //    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double w[] ),
            //    an array of pointers to functions which return the 1D quadrature points 
            //    associated with each spatial dimension for which a Golub Welsch rule 
            //    is used.
            //
            //    Input, double TOL, a tolerance for point equality.
            //
        {
            double alpha;
            double beta;
            int dim;
            int i;
            int level_max;
            int level_min;
            int p_index;
            int point_num;
            int point_total_num;

            Console.WriteLine("");
            Console.WriteLine("SGMG_SIZE_TEST");
            Console.WriteLine("  SGMG_SIZE returns the number of distinct");
            Console.WriteLine("  points in a multidimensional sparse grid with mixed factors.");
            Console.WriteLine("");
            Console.WriteLine("  SGMG_SIZE_TOTAL returns the TOTAL number of");
            Console.WriteLine("  points in a multidimensional sparse grid with mixed factors,");
            Console.WriteLine("  without checking for duplication.");
            Console.WriteLine("");
            Console.WriteLine("  Each sparse grid is of spatial dimension DIM_NUM,");
            Console.WriteLine("  and is made up of product grids of levels up to LEVEL_MAX.");
            Console.WriteLine("");
            Console.WriteLine(" Dimension      Rule  Growth rate       Parameters");
            Console.WriteLine("");

            p_index = 0;

            for (dim = 0; dim < dim_num; dim++)
            {
                if (rule[dim] == 1)
                {
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                }
                else if (rule[dim] == 2)
                {
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                }
                else if (rule[dim] == 3)
                {
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                }
                else if (rule[dim] == 4)
                {
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                }
                else if (rule[dim] == 5)
                {
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                }
                else if (rule[dim] == 6)
                {
                    alpha = p[p_index];
                    p_index = p_index + 1;
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8)
                                           + "  " + alpha.ToString().PadLeft(14) + "");
                }
                else if (rule[dim] == 7)
                {
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                }
                else if (rule[dim] == 8)
                {
                    alpha = p[p_index];
                    p_index = p_index + 1;
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8)
                                           + "  " + alpha.ToString().PadLeft(14) + "");
                }
                else if (rule[dim] == 9)
                {
                    alpha = p[p_index];
                    p_index = p_index + 1;
                    beta = p[p_index];
                    p_index = p_index + 1;
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8)
                                           + "  " + alpha.ToString().PadLeft(14)
                                           + "  " + beta.ToString().PadLeft(14) + "");
                }
                else if (rule[dim] == 10)
                {
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                }
                else if (rule[dim] == 11)
                {
                    string cout = "  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8);
                    for (i = 0; i < np[dim]; i++)
                    {
                        alpha = p[p_index];
                        p_index = p_index + 1;
                        cout += "  " + alpha.ToString().PadLeft(14);
                    }

                    Console.WriteLine(cout);
                }
                else if (rule[dim] == 12)
                {
                    string cout = "  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8);
                    for (i = 0; i < np[dim]; i++)
                    {
                        alpha = p[p_index];
                        p_index = p_index + 1;
                        cout += "  " + alpha.ToString().PadLeft(14);
                    }

                    Console.WriteLine(cout);
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("SGMG_SIZE_TEST - Fatal error!");
                    Console.WriteLine("  Unexpected value of RULE = " + rule[dim] + "");
                    return;
                }
            }

            Console.WriteLine("");
            Console.WriteLine(" LEVEL_MIN LEVEL_MAX POINT_NUM POINT_NUM");
            Console.WriteLine("                        Unique     Total");
            Console.WriteLine("");

            for (level_max = level_max_min; level_max <= level_max_max; level_max++)
            {
                point_total_num = SGMG.sgmg_size_total(dim_num,
                    level_max, rule, growth);

                point_num = SGMG.sgmg_size(dim_num, level_max,
                    rule, np, p, gw_compute_points, tol, growth);

                level_min = Math.Max(0, level_max + 1 - dim_num);

                Console.WriteLine("  " + level_min.ToString().PadLeft(8)
                                       + "  " + level_max.ToString().PadLeft(8)
                                       + "  " + point_num.ToString().PadLeft(8)
                                       + "  " + point_total_num.ToString().PadLeft(8) + "");
            }

        }
    }
}