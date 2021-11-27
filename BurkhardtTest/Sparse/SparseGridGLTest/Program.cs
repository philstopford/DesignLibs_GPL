﻿using System;
using System.Globalization;
using Burkardt.Composition;
using Burkardt.Quadrature;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SparseGridGLTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPARSE_GRID_GL_TEST.
        //
        //  Discussion:
        //
        //    SPARSE_GRID_GL_TEST tests the SPARSE_GRID_GL library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 October 2007
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
        Console.WriteLine("SPARSE_GRID_GL_TEST");
        Console.WriteLine("  Test the SPARSE_GRID_GL library.");
        //
        //  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
        //
        int dim_min = 1;
        int dim_max = 6;
        int level_max_min = 0;
        int level_max_max = 6;

        test01(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 6;
        dim_max = 10;
        level_max_min = 0;
        level_max_max = 5;

        test01(dim_min, dim_max, level_max_min, level_max_max);

        dim_min = 100;
        dim_max = 100;
        level_max_min = 0;
        level_max_max = 2;

        test01(dim_min, dim_max, level_max_min, level_max_max);
        //
        //  Compute abstract grid indices of sparse grid points as selected from product grid
        //  for DIMENSION, LEVEL_MAX.
        //
        int dim_num = 2;
        int level_max = 3;
        test02(dim_num, level_max);

        dim_num = 2;
        level_max = 4;
        test02(dim_num, level_max);

        dim_num = 3;
        level_max = 0;
        test02(dim_num, level_max);

        dim_num = 3;
        level_max = 2;
        test02(dim_num, level_max);

        dim_num = 6;
        level_max = 2;
        test02(dim_num, level_max);
        //
        //  Compute sparse Gauss-Legendre rule for DIMENSION, LEVEL_MAX.
        //
        dim_num = 2;
        level_max = 3;
        test03(dim_num, level_max);

        dim_num = 2;
        level_max = 4;
        test03(dim_num, level_max);

        dim_num = 3;
        level_max = 0;
        test03(dim_num, level_max);

        dim_num = 3;
        level_max = 2;
        test03(dim_num, level_max);
        //
        //  Test sum of weights for DIMENSION, LEVEL_MAX.
        //
        dim_num = 2;
        level_max = 4;
        test04(dim_num, level_max);

        dim_num = 3;
        level_max = 0;
        test04(dim_num, level_max);

        dim_num = 3;
        level_max = 1;
        test04(dim_num, level_max);

        dim_num = 3;
        level_max = 6;
        test04(dim_num, level_max);

        dim_num = 10;
        level_max = 3;
        test04(dim_num, level_max);
        //
        //  Test monomial exactness for DIMENSION, LEVEL_MAX, DEGREE_MAX.
        //
        dim_num = 2;
        level_max = 0;
        int degree_max = 3;
        test05(dim_num, level_max, degree_max);

        dim_num = 2;
        level_max = 1;
        degree_max = 5;
        test05(dim_num, level_max, degree_max);

        dim_num = 2;
        level_max = 2;
        degree_max = 10;
        test05(dim_num, level_max, degree_max);

        dim_num = 2;
        level_max = 3;
        degree_max = 14;
        test05(dim_num, level_max, degree_max);

        dim_num = 2;
        level_max = 4;
        degree_max = 22;
        test05(dim_num, level_max, degree_max);

        dim_num = 2;
        level_max = 5;
        degree_max = 32;
        test05(dim_num, level_max, degree_max);

        dim_num = 3;
        level_max = 0;
        degree_max = 2;
        test05(dim_num, level_max, degree_max);

        dim_num = 3;
        level_max = 1;
        degree_max = 4;
        test05(dim_num, level_max, degree_max);

        dim_num = 3;
        level_max = 2;
        degree_max = 6;
        test05(dim_num, level_max, degree_max);

        dim_num = 3;
        level_max = 3;
        degree_max = 12;
        test05(dim_num, level_max, degree_max);
        //
        //  Show how to write a rule to a file.
        //
        dim_num = 2;
        level_max = 3;

        test06(dim_num, level_max);

        Console.WriteLine("");
        Console.WriteLine("SPARSE_GRID_GL_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int dim_min, int dim_max, int level_max_min, int level_max_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests SPARSE_GRID_GL_SIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 August 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_MIN, the minimum spatial dimension to consider.
        //
        //    Input, int DIM_MAX, the maximum spatial dimension to consider.
        //
        //    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
        //
        //    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
        //
    {
        int dim_num;
        int level_max;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  SPARSE_GRID_GL_SIZE returns the number of distinct");
        Console.WriteLine("  points in a Gauss-Legendre sparse grid.");
        Console.WriteLine("");
        Console.WriteLine("  Note that, unlike most sparse grids, a sparse grid based on");
        Console.WriteLine("  Gauss-Legendre points is almost entirely NOT nested.");
        Console.WriteLine("");
        Console.WriteLine("  Hence the point counts should be much higher than for a grid of");
        Console.WriteLine("  the same level, but using rules such as Fejer1 or Fejer2 or");
        Console.WriteLine("  Gauss-Patterson or Newton-Cotes-Open or Newton-Cotes-Open-Half.");
        Console.WriteLine("");
        Console.WriteLine("  Each sparse grid is of spatial dimension DIM,");
        Console.WriteLine("  and is made up of all product grids of levels up to LEVEL_MAX.");
        Console.WriteLine("");
        string cout = "   DIM: ";
        for (dim_num = dim_min; dim_num <= dim_max; dim_num++)
        {
            cout += "  " + dim_num.ToString(CultureInfo.InvariantCulture).PadLeft(8);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("   LEVEL_MAX");
        Console.WriteLine("");

        for (level_max = level_max_min; level_max <= level_max_max; level_max++)
        {
            cout = "    " + level_max.ToString(CultureInfo.InvariantCulture).PadLeft(4);
            for (dim_num = dim_min; dim_num <= dim_max; dim_num++)
            {
                int point_num = Grid_GaussLegendre.sparse_grid_gl_size(dim_num, level_max);
                cout += "  " + point_num.ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test02(int dim_num, int level_max)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests SPARSE_GRID_GL_INDEX.
        //
        //  Discussion:
        //
        //    The routine computes abstract indices that describe the sparse grid
        //    of GL points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the level.
        //
    {
        int point;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  SPARSE_GRID_GL_INDEX returns abstract indices for the");
        Console.WriteLine("  points that make up a Gauss-Legendre sparse grid.");

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        Console.WriteLine("");
        Console.WriteLine("  LEVEL_MIN = " + level_min + "");
        Console.WriteLine("  LEVEL_MAX = " + level_max + "");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");

        int point_num = Grid_GaussLegendre.sparse_grid_gl_size(dim_num, level_max);

        Console.WriteLine("");
        Console.WriteLine("  Number of unique points in the grid = " + point_num + "");
        //
        //  Compute the orders and points.
        //
        int[] grid_index = new int[dim_num * point_num];
        int[] grid_base = new int[dim_num * point_num];

        Grid_GaussLegendre.sparse_grid_gl_index(dim_num, level_max, point_num, ref grid_index,
            ref grid_base);
        //
        //  Now we're done.  Print the merged grid data.
        //
        Console.WriteLine("");
        Console.WriteLine("  Grid index/base:");
        Console.WriteLine("");
        for (point = 0; point < point_num; point++)
        {
            string cout = "  " + point.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
            int dim;
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += grid_index[dim + point * dim_num].ToString(CultureInfo.InvariantCulture).PadLeft(6);
            }

            Console.WriteLine(cout);
            cout = "        ";
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += grid_base[dim + point * dim_num].ToString(CultureInfo.InvariantCulture).PadLeft(6);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test03(int dim_num, int level_max)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 call SPARSE_GRID_GL to create a sparse Gauss-Legendre grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 September 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the level.
        //
    {
        int point;

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  SPARSE_GRID_GL makes a sparse Gauss-Legendre grid.");

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        Console.WriteLine("");
        Console.WriteLine("  LEVEL_MIN = " + level_min + "");
        Console.WriteLine("  LEVEL_MAX = " + level_max + "");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        //
        //  Determine the number of points.
        //
        int point_num = Grid_GaussLegendre.sparse_grid_gl_size(dim_num, level_max);

        Console.WriteLine("");
        Console.WriteLine("  Number of unique points in the grid = " + point_num + "");
        //
        //  Allocate space for the weights and points.
        //
        double[] grid_weight = new double[point_num];
        double[] grid_point = new double[dim_num * point_num];
        //
        //  Compute the weights and points.
        //
        Grid_GaussLegendre.sparse_grid_gl(dim_num, level_max, point_num, ref grid_weight, ref grid_point);
        //
        //  Print them out.
        //
        Console.WriteLine("");
        Console.WriteLine("  Grid weights:");
        Console.WriteLine("");
        for (point = 0; point < point_num; point++)
        {
            Console.WriteLine("  " + point.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  "
                                   + grid_weight[point].ToString("0.######").PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Grid points:");
        Console.WriteLine("");
        for (point = 0; point < point_num; point++)
        {
            string cout = "  " + point.ToString(CultureInfo.InvariantCulture).PadLeft(4);
            int dim;
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += "  "
                        + grid_point[dim + point * dim_num].ToString("0.######").PadLeft(9);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test04(int dim_num, int level_max)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 sums the weights and compares them to 2^DIM_NUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 September 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the level.
        //
    {
        int point;

        Console.WriteLine("");
        Console.WriteLine("TEST04:");
        Console.WriteLine("  Compute the weights of a Gauss-Legendre sparse grid .");
        Console.WriteLine("");
        Console.WriteLine("  As a simple test, sum these weights.");
        Console.WriteLine("  They should sum to exactly 2^DIM_NUM.");

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        Console.WriteLine("");
        Console.WriteLine("  LEVEL_MIN = " + level_min + "");
        Console.WriteLine("  LEVEL_MAX = " + level_max + "");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        //
        //  Determine the number of points.
        //
        int point_num = Grid_GaussLegendre.sparse_grid_gl_size(dim_num, level_max);

        Console.WriteLine("");
        Console.WriteLine("  Number of unique points in the grid = " + point_num + "");
        //
        //  Allocate space for the weights and points.
        //
        double[] grid_weight = new double[point_num];
        double[] grid_point = new double[dim_num * point_num];
        //
        //  Compute the weights and points.
        //
        Grid_GaussLegendre.sparse_grid_gl(dim_num, level_max, point_num, ref grid_weight, ref grid_point);
        //
        //  Sum the weights.
        //
        double weight_sum = 0.0;
        for (point = 0; point < point_num; point++)
        {
            weight_sum += grid_weight[point];
        }

        double weight_sum_exact = Math.Pow(2, dim_num);

        double weight_sum_error = Math.Abs(weight_sum - weight_sum_exact);

        Console.WriteLine("");
        Console.WriteLine("    Weight sum     Exact sum    Difference");
        Console.WriteLine("");
        Console.WriteLine(
            "  " + weight_sum.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                 + "  " + weight_sum_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                 + "  " + weight_sum_error.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void test05(int dim_num, int level_max, int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests a Gauss-Legendre sparse grid rule for monomial exactness.
        //
        //  Discussion:
        //
        //    This test is going to check EVERY monomial of total degree DEGREE_MAX
        //    or less.  Even for a moderately high dimension of DIM_NUM = 10, you
        //    do NOT want to use a large value of DEGREE_MAX, since there are
        //
        //      1         monomials of total degree 0,
        //      DIM_NUM   monomials of total degree 1,
        //      DIM_NUM^2 monomials of total degree 2,
        //      DIM_NUM^3 monomials of total degree 3, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the level.
        //
        //    Input, int DEGREE_MAX, the maximum monomial total degree to check.
        //
    {
        int degree;
        int dim;
        int point;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Check the exactness of a Gauss-Legendre sparse");
        Console.WriteLine("  grid quadrature rule, applied to all monomials ");
        Console.WriteLine("  of orders 0 to DEGREE_MAX.");

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        Console.WriteLine("");
        Console.WriteLine("  LEVEL_MIN = " + level_min + "");
        Console.WriteLine("  LEVEL_MAX = " + level_max + "");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("");
        Console.WriteLine("  The maximum total degree to be checked is DEGREE_MAX = " + degree_max + "");
        //
        //  Determine the number of points in the rule.
        //
        int point_num = Grid_GaussLegendre.sparse_grid_gl_size(dim_num, level_max);

        Console.WriteLine("");
        Console.WriteLine("  Number of unique points in the grid = " + point_num + "");
        //
        //  Allocate space for the weights and points.
        //
        double[] grid_weight = new double[point_num];
        double[] grid_point = new double[dim_num * point_num];
        //
        //  Compute the weights and points.
        //
        Grid_GaussLegendre.sparse_grid_gl(dim_num, level_max, point_num, ref grid_weight, ref grid_point);
        //
        //  Rescale the weights, and translate the abscissas.
        //
        double volume = Math.Pow(2, dim_num);

        for (point = 0; point < point_num; point++)
        {
            grid_weight[point] /= volume;
        }

        for (dim = 0; dim < dim_num; dim++)
        {
            for (point = 0; point < point_num; point++)
            {
                grid_point[dim + point * dim_num] = (grid_point[dim + point * dim_num] + 1.0) / 2.0;
            }
        }

        //
        //  Explore the monomials.
        //
        int[] expon = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine("      Error      Total   Monomial");
        Console.WriteLine("                 Degree  Exponents");

        for (degree = 0; degree <= degree_max; degree++)
        {
            bool more = false;
            int h = 0;
            int t = 0;

            Console.WriteLine("");
            for (;;)
            {
                Comp.comp_next(degree, dim_num, ref expon, ref more, ref h, ref t);

                double quad_error = MonomialQuadrature.monomial_quadrature(dim_num, expon, point_num,
                    grid_weight, grid_point);

                string cout = "  " + quad_error.ToString("0.#").PadLeft(12)
                                   + "     " + degree.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "      ";

                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += expon[dim].ToString(CultureInfo.InvariantCulture).PadLeft(3);
                }

                Console.WriteLine(cout);

                if (!more)
                {
                    break;
                }
            }
        }
    }

    private static void test06(int dim_num, int level_max)

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 creates a sparse Gauss-Legendre grid and writes it to a file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL_MAX, the level.
        //
    {
        int dim;

        Console.WriteLine("");
        Console.WriteLine("TEST06:");
        Console.WriteLine("  Call SPARSE_GRID_GL to make a sparse Gauss-Legendre grid.");
        Console.WriteLine("  Write the data to a set of quadrature files.");

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        Console.WriteLine("");
        Console.WriteLine("  LEVEL_MIN = " + level_min + "");
        Console.WriteLine("  LEVEL_MAX = " + level_max + "");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        //
        //  Determine the number of points.
        //
        int point_num = Grid_GaussLegendre.sparse_grid_gl_size(dim_num, level_max);
        //
        //  Allocate space for the weights and points.
        //
        double[] r = new double[dim_num * 2];
        double[] w = new double[point_num];
        double[] x = new double[dim_num * point_num];
        //
        //  Compute the weights and points.
        //
        for (dim = 0; dim < dim_num; dim++)
        {
            r[dim + 0 * dim_num] = -1.0;
            r[dim + 1 * dim_num] = +1.0;
        }

        Grid_GaussLegendre.sparse_grid_gl(dim_num, level_max, point_num, ref w, ref x);
        //
        //  Write the data out.
        //
        string r_filename = "gl_d" + dim_num + "_level"
                            + level_max + "_r.txt";
        string w_filename = "gl_d" + dim_num + "_level"
                            + level_max + "_w.txt";
        string x_filename = "gl_d" + (dim_num, "%d") + "_level"
                            + level_max + "_x.txt";

        typeMethods.r8mat_write(r_filename, dim_num, 2, r);
        typeMethods.r8mat_write(w_filename, 1, point_num, w);
        typeMethods.r8mat_write(x_filename, dim_num, point_num, x);

        Console.WriteLine("");
        Console.WriteLine("  R data written to \"" + r_filename + "\".");
        Console.WriteLine("  W data written to \"" + w_filename + "\".");
        Console.WriteLine("  X data written to \"" + x_filename + "\".");
    }
}