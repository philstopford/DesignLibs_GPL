using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Quadrature;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SGMGAnisoWriteTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SGMGA_WRITE_TEST.
            //
            //  Discussion:
            //
            //    SGMGA_WRITE_TEST tests the SGMGA_WRITE function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 November 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("SGMGA_WRITE_TEST");
            Console.WriteLine("  Test the SGMGA_WRITE function.");
            //
            //  Generate sparse grid rules and write them to files.
            //
            sgmga_write_tests();
            //
            //  That's all.
            //
            Console.WriteLine("");
            Console.WriteLine("SGMGA_WRITE_TEST");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }

        static void sgmga_write_tests()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SGMGA_WRITE_TESTS calls SGMGA_WRITE_TEST.
            //
            //  Discussion:
            //  
            //    We can't test Golub-Welsch rules in this routine, because the program
            //    that writes out the files needs to know the integration region for each
            //    component, and we have not specified how that would be done with 
            //    Golub Welsch rules.
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
            //  Local Parameters:
            //
            //    Local, double TOL, a tolerance for point equality.
            //    A value of sqrt ( eps ) is reasonable, and will allow the code to
            //    consolidate points which are equal, or very nearly so.  A value of
            //    -1.0, on the other hand, will force the code to use every point, 
            //    regardless of duplication.
            //
        {
            int dim;
            int dim_num;
            string file_name;
            int[] growth;
            Func<int, int, double[], double[], double[]>[] gw_compute_points;
            Func<int, int, double[], double[], double[]>[] gw_compute_weights;
            double[] importance;
            int level_max;
            int level_max_max;
            int level_max_min;
            double[] level_weight;
            int[] np;
            int np_sum;
            int[] order_1d;
            int order_nd;
            double[] p;
            int[] rule;
            double tol;

            Console.WriteLine("");
            Console.WriteLine("SGMGA_WRITE_TESTS");
            Console.WriteLine("  Call SGMGA_WRITE_TEST with various arguments.");
            //
            //  Set the point equality tolerance.
            //
            tol = Math.Sqrt(typeMethods.r8_epsilon());
            Console.WriteLine("");
            Console.WriteLine("  All tests will use a point equality tolerance of " + tol + "");

            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = 1.0;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 1;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 6;
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
            file_name = "sgmga_d2_l2_ccxcc_iso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);

            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 1;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 6;
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
            file_name = "sgmga_d2_l2_ccxcc_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);

            dim_num = 3;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = 1.0;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 1;
            rule[2] = 1;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 6;
            growth[2] = 6;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 0;
            np[2] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_points[1] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_points[2] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            gw_compute_weights[1] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            gw_compute_weights[2] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            file_name = "sgmga_d3_l2_ccxccxcc_iso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);

            dim_num = 3;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 1;
            rule[2] = 1;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 6;
            growth[2] = 6;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 0;
            np[2] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_points[1] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_points[2] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            gw_compute_weights[1] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            gw_compute_weights[2] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            file_name = "sgmga_d3_l2_ccxccxcc_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);

            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 3;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 3;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 6;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_points[1] = PattersonQuadrature.patterson_lookup_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            gw_compute_weights[1] = PattersonQuadrature.patterson_lookup_weights_np;
            file_name = "sgmga_d2_l3_ccxgp_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);

            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 4;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 3;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_points[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            gw_compute_weights[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_weights_np;
            file_name = "sgmga_d2_l2_ccxgl_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);

            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 7;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 3;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_points[1] = Burkardt.Laguerre.QuadratureRule.laguerre_compute_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            gw_compute_weights[1] = Burkardt.Laguerre.QuadratureRule.laguerre_compute_weights_np;
            file_name = "sgmga_d2_l2_ccxlg_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);

            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 8;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 3;
            np = new int[dim_num];
            np[0] = 1;
            np[1] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            p[0] = 1.5;
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_points[1] = Burkardt.Laguerre.QuadratureRule.gen_laguerre_compute_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            gw_compute_weights[1] = Burkardt.Laguerre.QuadratureRule.gen_laguerre_compute_weights_np;
            file_name = "sgmga_d2_l2_ccxglg_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);

            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            rule = new int[dim_num];
            rule[0] = 2;
            rule[1] = 9;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 3;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 2;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            p[0] = 0.5;
            p[1] = 1.5;
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = Fejer2.fejer2_compute_points_np;
            gw_compute_points[1] = JacobiQuadrature.jacobi_compute_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = Fejer2.fejer2_compute_weights_np;
            gw_compute_weights[1] = JacobiQuadrature.jacobi_compute_weights_np;
            file_name = "sgmga_d2_l2_f2xgj_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);

            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            rule = new int[dim_num];
            rule[0] = 6;
            rule[1] = 10;
            growth = new int[dim_num];
            growth[0] = 3;
            growth[1] = 4;
            np = new int[dim_num];
            np[0] = 1;
            np[1] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            p[0] = 2.0;
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = HermiteQuadrature.gen_hermite_compute_points_np;
            gw_compute_points[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = HermiteQuadrature.gen_hermite_compute_weights_np;
            gw_compute_weights[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_weights_np;
            file_name = "sgmga_d2_l2_gghxhgk_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);
            //
            //  LEVEL_MAX = 1
            //
            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 1;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 1;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 6;
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
            file_name = "sgmga_d2_l1_ccxcc_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);
            //
            //  LEVEL_MAX = 2 (already done)
            //
            //  LEVEL_MAX = 3
            //
            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 3;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 1;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 6;
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
            file_name = "sgmga_d2_l3_ccxcc_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);
            //
            //  LEVEL_MAX = 4
            //
            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 4;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 1;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 6;
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
            file_name = "sgmga_d2_l4_ccxcc_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);
            //
            //  LEVEL_MAX = 5
            //
            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 5;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 1;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 6;
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
            file_name = "sgmga_d2_l5_ccxcc_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);
            //
            //  Dimension 3
            //
            dim_num = 3;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = (double) (dim + 1);
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 2;
            rule = new int[dim_num];
            rule[0] = 1;
            rule[1] = 4;
            rule[2] = 5;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 3;
            growth[2] = 3;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 0;
            np[2] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
            gw_compute_points[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_points_np;
            gw_compute_points[2] = HermiteQuadrature.hermite_compute_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = ClenshawCurtis.clenshaw_curtis_compute_weights_np;
            gw_compute_weights[1] = Burkardt.Legendre.QuadratureRule.legendre_compute_weights_np;
            gw_compute_weights[2] = HermiteQuadrature.hermite_compute_weights_np;
            file_name = "sgmga_d3_l2_ccxglxgh_aniso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);
            //
            //  Rule 3, LEVEL_MAX = 4
            //
            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = 1.0;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 4;
            rule = new int[dim_num];
            rule[0] = 3;
            rule[1] = 3;
            growth = new int[dim_num];
            growth[0] = 6;
            growth[1] = 6;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = PattersonQuadrature.patterson_lookup_points_np;
            gw_compute_points[1] = PattersonQuadrature.patterson_lookup_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = PattersonQuadrature.patterson_lookup_weights_np;
            gw_compute_weights[1] = PattersonQuadrature.patterson_lookup_weights_np;
            file_name = "sgmga_d2_l4_gpxgp_iso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);
            //
            //  Rule 3, Slow Exponential Growth, LEVEL_MAX = 4
            //
            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = 1.0;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 4;
            rule = new int[dim_num];
            rule[0] = 3;
            rule[1] = 3;
            growth = new int[dim_num];
            growth[0] = 4;
            growth[1] = 4;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = PattersonQuadrature.patterson_lookup_points_np;
            gw_compute_points[1] = PattersonQuadrature.patterson_lookup_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = PattersonQuadrature.patterson_lookup_weights_np;
            gw_compute_weights[1] = PattersonQuadrature.patterson_lookup_weights_np;
            file_name = "sgmga_d2_l4_gpsexgpse_iso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);
            //
            //  Rule 3, Moderate Exponential Growth, LEVEL_MAX = 4
            //
            dim_num = 2;
            importance = new double[dim_num];
            for (dim = 0; dim < dim_num; dim++)
            {
                importance[dim] = 1.0;
            }

            level_weight = new double[dim_num];
            SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
            level_max = 4;
            rule = new int[dim_num];
            rule[0] = 3;
            rule[1] = 3;
            growth = new int[dim_num];
            growth[0] = 5;
            growth[1] = 5;
            np = new int[dim_num];
            np[0] = 0;
            np[1] = 0;
            np_sum = typeMethods.i4vec_sum(dim_num, np);
            p = new double[np_sum];
            gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_points[0] = PattersonQuadrature.patterson_lookup_points_np;
            gw_compute_points[1] = PattersonQuadrature.patterson_lookup_points_np;
            gw_compute_weights = new Func<int, int, double[], double[], double[]>[dim_num];
            gw_compute_weights[0] = PattersonQuadrature.patterson_lookup_weights_np;
            gw_compute_weights[1] = PattersonQuadrature.patterson_lookup_weights_np;
            file_name = "sgmga_d2_l4_gpmexgpme_iso";
            sgmga_write_test(dim_num, level_weight, level_max, rule, growth, np,
                p, gw_compute_points, gw_compute_weights, tol, file_name);

            return;
        }

        static void sgmga_write_test(int dim_num, double[] level_weight, int level_max,
                int[] rule, int[] growth, int[] np, double[] p,
                Func<int, int, double[], double[], double[]>[] gw_compute_points,
                Func<int, int, double[], double[], double[]>[] gw_compute_weights,
                double tol, string file_name)

            //***************************************************************************80
            //
            //  Purpose:
            //
            //    SGMGA_WRITE_TEST tests SGMGA_WRITE.
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
            //    Input, integer DIM_NUM, the spatial dimension.
            //
            //    Input, double LEVEL_WEIGHT[DIM_NUM], the weights for each dimension.
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
            //    Input, int GROWTH[DIM_NUM], the desired growth in each dimension.
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
            Console.WriteLine("SGMGA_WRITE_TEST");
            Console.WriteLine("  SGMGA_WRITE writes a sparse grid rule to files.");
            //
            //  Compute necessary data.
            //
            point_total_num = SGMGAniso.sgmga_size_total(dim_num, level_weight,
                level_max, rule, growth);

            point_num = SGMGAniso.sgmga_size(dim_num, level_weight, level_max,
                rule, np, p, gw_compute_points, tol, growth);

            sparse_unique_index = new int[point_total_num];

            SGMGAniso.sgmga_unique_index(dim_num, level_weight, level_max, rule,
                np, p, gw_compute_points, tol, point_num, point_total_num,
                growth, ref sparse_unique_index);

            sparse_order = new int[dim_num * point_num];
            sparse_index = new int[dim_num * point_num];

            SGMGAniso.sgmga_index(dim_num, level_weight, level_max, rule, point_num,
                point_total_num, sparse_unique_index, growth, ref sparse_order, ref sparse_index);
            //
            //  Compute points and weights.
            //
            sparse_point = new double [dim_num * point_num];

            SGMGAniso.sgmga_point(dim_num, level_weight, level_max, rule, np,
                p, gw_compute_points, point_num, sparse_order, sparse_index,
                growth, ref sparse_point);

            sparse_weight = new double[point_num];

            SGMGAniso.sgmga_weight(dim_num, level_weight, level_max, rule, np,
                p, gw_compute_weights, point_num, point_total_num, sparse_unique_index,
                growth, sparse_weight);
            //
            //  Write points and weights to files.
            //
            SGMGAniso.sgmga_write(dim_num, level_weight, rule, np, p,
                point_num, sparse_weight, sparse_point, file_name);

        }
    }
}