﻿using System;
using System.Globalization;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Quadrature;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SGMGAnisoSizeTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SGMGA_SIZE_TEST.
        //
        //  Discussion:
        //
        //    SGMGA_SIZE_TEST tests the SGMGA_SIZE and SGMGA_SIZE_TOTAL functions.
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
        Console.WriteLine("SGMGA_SIZE_TEST");
        Console.WriteLine("  Test the SGMGA_SIZE and SGMGA_SIZE_TOTAL functions.");

        sgmga_size_tests();

        Console.WriteLine("");
        Console.WriteLine("SGMGA__SIZE_TEST");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void sgmga_size_tests()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_SIZE_TESTS calls SGMGA_SIZE_TEST.
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

        Console.WriteLine("");
        Console.WriteLine("SGMGA_SIZE_TESTS");
        Console.WriteLine("  Call SGMGA_SIZE_TEST with various arguments.");
        //
        //  Set the point equality tolerance.
        //
        double tol = Math.Sqrt(typeMethods.r8_epsilon());
        Console.WriteLine("");
        Console.WriteLine("  Point equality tolerance = " + tol + "");

        int dim_num = 2;
        double[] importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = 1.0;
        }

        double[] level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        int level_max_min = 0;
        int level_max_max = 5;
        int[] rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 1;
        int[] growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 6;
        int[] np = new int[dim_num];
        np[0] = 0;
        np[1] = 0;
        int np_sum = typeMethods.i4vec_sum(dim_num, np);
        double[] p = new double[np_sum];
        Func<int, int, double[], double[], double[]>[] gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_points[1] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

        dim_num = 2;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 5;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

        dim_num = 3;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = 1.0;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 5;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

        dim_num = 3;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 5;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

        dim_num = 2;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 2;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

        dim_num = 2;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 2;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

        dim_num = 2;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 2;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

        dim_num = 2;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 2;
        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 8;
        growth = new int[dim_num];
        growth[0] = 6;
        growth[1] = 3;
        np = new int[dim_num];
        np[0] = 0;
        np[1] = 1;
        np_sum = typeMethods.i4vec_sum(dim_num, np);
        p = new double[np_sum];
        p[0] = 1.5;
        gw_compute_points = new Func<int, int, double[], double[], double[]>[dim_num];
        gw_compute_points[0] = ClenshawCurtis.clenshaw_curtis_compute_points_np;
        gw_compute_points[1] = Burkardt.Laguerre.QuadratureRule.gen_laguerre_compute_points_np;
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

        dim_num = 2;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 2;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

        dim_num = 2;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 2;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

        dim_num = 3;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 2;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);
        //
        //  Repeat, treating rules #2 and #3 as Golub Welsch rules.
        //
        dim_num = 3;
        importance = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            importance[dim] = dim + 1;
        }

        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 2;
        rule = new int[dim_num];
        rule[0] = 1;
        rule[1] = 11;
        rule[2] = 11;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);
        //
        //  Try a case involving a dimension of "0" importance.
        //
        dim_num = 3;
        importance = new double[dim_num];
        importance[0] = 1.0;
        importance[1] = 0.0;
        importance[2] = 1.0;
        level_weight = new double[dim_num];
        SGMGAniso.sgmga_importance_to_aniso(dim_num, importance, ref level_weight);
        level_max_min = 0;
        level_max_max = 3;
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
        sgmga_size_test(dim_num, importance, level_weight, level_max_min,
            level_max_max, rule, growth, np, p, gw_compute_points, tol);

    }

    private static void sgmga_size_test(int dim_num, double[] importance, double[] level_weight,
            int level_max_min, int level_max_max, int[] rule, int[] growth, int[] np,
            double[] p,
            Func<int, int, double[], double[], double[]>[] gw_compute_points,
            double tol)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    SGMGA_SIZE_TEST tests SGMG_SIZE, SGMG_SIZE_TOTAL.
        //
        //  Discussion:
        //
        //    The "level_to_order" argument for the SGMG_SIZE and SGMG_SIZE_TOTAL
        //    functions is set to "LEVEL_TO_ORDER_DEFAULT", which uses exponential
        //    growth for fully nested rules, and linear otherwise.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double IMPORTANCE[DIM_NUM], the importance for each dimension.
        //
        //    Input, double LEVEL_WEIGHT[DIM_NUM], the weights for each dimension.
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
        //    Input, double TOL, a tolerance for point equality.
        //
    {
        int dim;
        int level_max;

        Console.WriteLine("");
        Console.WriteLine("SGMGA_SIZE_TEST");
        Console.WriteLine("  SGMGA_SIZE_TOTAL counts the total number of points,");
        Console.WriteLine("  including duplications, in an SGMGA sparse grid.");
        Console.WriteLine("  SGMGA_SIZE counts the total number of points,");
        Console.WriteLine("  excluding duplications, in an SGMGA sparse grid.");
        Console.WriteLine("");
        string cout = "  IMPORTANCE:  ";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + importance[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
        cout = "  LEVEL_WEIGHT:";
        for (dim = 0; dim < dim_num; dim++)
        {
            cout += "  " + level_weight[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine(" Dimension      Rule  Growth rate       Parameters");
        Console.WriteLine("");

        int p_index = 0;

        for (dim = 0; dim < dim_num; dim++)
        {
            int i;
            double alpha;
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
                                           + "  " + alpha.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
                                           + "  " + alpha.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
                case 9:
                    alpha = p[p_index];
                    p_index += 1;
                    double beta = p[p_index];
                    p_index += 1;
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8)
                                           + "  " + alpha.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + beta.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
                case 10:
                    Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                           + "  " + rule[dim].ToString().PadLeft(8)
                                           + "  " + growth[dim].ToString().PadLeft(8) + "");
                    break;
                case 11:
                {
                    cout = "  " + dim.ToString().PadLeft(8)
                                + "  " + rule[dim].ToString().PadLeft(8)
                                + "  " + growth[dim].ToString().PadLeft(8);
                    for (i = 0; i < np[dim]; i++)
                    {
                        alpha = p[p_index];
                        p_index += 1;
                        cout += "  " + alpha.ToString(CultureInfo.InvariantCulture).PadLeft(14);
                    }

                    Console.WriteLine(cout);
                    break;
                }
                case 12:
                {
                    cout = "  " + dim.ToString().PadLeft(8)
                                + "  " + rule[dim].ToString().PadLeft(8)
                                + "  " + growth[dim].ToString().PadLeft(8);
                    for (i = 0; i < np[dim]; i++)
                    {
                        alpha = p[p_index];
                        p_index += 1;
                        cout += "  " + alpha.ToString(CultureInfo.InvariantCulture).PadLeft(14);
                    }

                    Console.WriteLine(cout);
                    break;
                }
                default:
                    Console.WriteLine("");
                    Console.WriteLine("SGMGA_SIZE_TEST - Fatal error!");
                    Console.WriteLine("  Unexpected value of RULE = " + rule[dim] + "");
                    return;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("   DIM_NUM LEVEL_MAX POINT_NUM POINT_NUM");
        Console.WriteLine("                        Unique     Total");
        Console.WriteLine("");

        for (level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            int point_total_num = SGMGAniso.sgmga_size_total(dim_num, level_weight,
                level_max, rule, growth);

            int point_num = SGMGAniso.sgmga_size(dim_num, level_weight, level_max,
                rule, np, p, gw_compute_points, tol, growth);

            Console.WriteLine("  " + dim_num.ToString().PadLeft(8)
                                   + "  " + level_max.ToString().PadLeft(8)
                                   + "  " + point_num.ToString().PadLeft(8)
                                   + "  " + point_total_num.ToString().PadLeft(8) + "");
        }
    }
}