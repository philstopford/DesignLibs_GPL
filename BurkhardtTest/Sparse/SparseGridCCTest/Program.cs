using System;
using Burkardt.Composition;
using Burkardt.Quadrature;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SparseGridCCTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SPARSE_GRID_CC_TEST.
            //
            //  Discussion:
            //
            //    SPARSE_GRID_CC_TEST tests the SPARSE_GRID_CC library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 March 2013
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
            int dim_max;
            int dim_min;
            int dim_num;
            int level_max;
            int level_max_max;
            int level_max_min;

            Console.WriteLine("");
            Console.WriteLine("SPARSE_GRID_CC_TEST");
            Console.WriteLine("  Test the SPARSE_GRID_CC library.");
            //
            //  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
            //
            dim_min = 1;
            dim_max = 5;
            level_max_min = 0;
            level_max_max = 10;
            test01(dim_min, dim_max, level_max_min, level_max_max);

            Console.WriteLine("");

            dim_min = 6;
            dim_max = 10;
            level_max_min = 0;
            level_max_max = 10;
            test01(dim_min, dim_max, level_max_min, level_max_max);

            Console.WriteLine("");

            dim_min = 100;
            dim_max = 100;
            level_max_min = 0;
            level_max_max = 4;
            test01(dim_min, dim_max, level_max_min, level_max_max);

            Console.WriteLine("");
            //
            //  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
            //
            dim_min = 1;
            dim_max = 5;
            level_max_min = 0;
            level_max_max = 10;
            test015(dim_min, dim_max, level_max_min, level_max_max);

            Console.WriteLine("");

            dim_min = 6;
            dim_max = 10;
            level_max_min = 0;
            level_max_max = 10;
            test015(dim_min, dim_max, level_max_min, level_max_max);

            Console.WriteLine("");

            dim_min = 100;
            dim_max = 100;
            level_max_min = 0;
            level_max_max = 4;
            test015(dim_min, dim_max, level_max_min, level_max_max);

            Console.WriteLine("");
            //
            //  Compute abstract grid indices of sparse grid points as selected from product grid
            //  for DIMENSION, LEVEL_MAX.
            //
            test02(2, 3);
            test02(2, 4);
            test02(3, 0);
            test02(3, 2);
            test02(6, 2);
            //
            //  Compute sparse Clenshaw-Curtis rule for DIMENSION, LEVEL_MAX.
            //
            test03(2, 3);
            test03(3, 0);
            test03(3, 1);
            //
            //  Test sum of weights for DIMENSION, LEVEL_MAX.
            //
            test04(2, 4);
            test04(3, 0);
            test04(3, 1);
            test04(3, 6);
            test04(10, 3);
            //
            //  Test monomial exactness for DIMENSION, LEVEL_MAX, DEGREE_MAX.
            // 
            test05(2, 0, 3);
            test05(2, 1, 5);
            test05(2, 2, 7);
            test05(2, 3, 9);
            test05(2, 4, 11);
            test05(2, 5, 13);

            test05(3, 0, 2);
            test05(3, 1, 4);
            test05(3, 2, 6);
            test05(3, 3, 8);
            //
            //  Show how to write a rule to a file.
            //
            dim_num = 2;
            level_max = 3;

            test06(dim_num, level_max);

            Console.WriteLine("");
            Console.WriteLine("SPARSE_GRID_CC_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01(int dim_min, int dim_max, int level_max_min, int level_max_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests SPARSE_GRID_CFN_SIZE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 December 2009
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
            int point_num;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  SPARSE_GRID_CFN_SIZE returns the number of distinct");
            Console.WriteLine("  points in a sparse grid of Closed Fully Nested rules.");
            Console.WriteLine("");
            Console.WriteLine("  Each sparse grid is of spatial dimension DIM,");
            Console.WriteLine("  and is made up of all product grids of levels up to LEVEL_MAX.");
            Console.WriteLine("");
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
                    point_num = Grid_ClenshawCurtis.sparse_grid_cfn_size(dim_num, level_max);
                    cout += "  " + point_num.ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
            }
        }

        static void test015(int dim_min, int dim_max, int level_max_min, int level_max_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST015 tests SPARSE_GRID_CCS_SIZE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 December 2009
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
            int point_num;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  SPARSE_GRID_CCS_SIZE returns the number of distinct");
            Console.WriteLine("  points in a Clenshaw Curtis Slow-Growth sparse grid.");
            Console.WriteLine("");
            Console.WriteLine("  Each sparse grid is of spatial dimension DIM,");
            Console.WriteLine("  and is made up of all product grids of levels up to LEVEL_MAX.");
            Console.WriteLine("");
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
                    point_num = Grid_ClenshawCurtis.sparse_grid_ccs_size(dim_num, level_max);
                    cout += "  " + point_num.ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
            }
        }

        static void test02(int dim_num, int level_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests SPARSE_GRID_CC_INDEX.
            //
            //  Discussion:
            //
            //    The routine computes the indices of the unique points used in a sparse 
            //    multidimensional grid whose size is controlled by a parameter LEVEL_MAX.
            //
            //    Once these indices are returned, they can be converted into
            //    Clenshaw Curtis points.
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
            //    Input, int LEVEL_MAX, the level.
            //
        {
            int dim;
            int[] grid_index;
            int point;
            int point_num;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  SPARSE_GRID_CC_INDEX returns all grid indexes");
            Console.WriteLine("  whose level value satisfies");
            Console.WriteLine("    0 <= LEVEL <= LEVEL_MAX.");
            Console.WriteLine("  Here, LEVEL is the sum of the levels of the 1D rules,");
            Console.WriteLine("  and the order of the rule is 2**LEVEL + 1.");

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");

            point_num = Grid_ClenshawCurtis.sparse_grid_cfn_size(dim_num, level_max);

            Console.WriteLine("");
            Console.WriteLine("  Number of unique points in the grid = " + point_num + "");
            //
            //  Compute the orders and points.
            //
            grid_index = Grid_ClenshawCurtis.sparse_grid_cc_index(dim_num, level_max, point_num);
            //
            //  Now we're done.  Print the merged grid data.
            //
            Console.WriteLine("");
            Console.WriteLine("  Grid index:");
            Console.WriteLine("");
            for (point = 0; point < point_num; point++)
            {
                string cout = "  " + point.ToString().PadLeft(4) + "  ";
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += grid_index[dim + point * dim_num].ToString().PadLeft(6);
                }

                Console.WriteLine(cout);
            }
        }

        static void test03(int dim_num, int level_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 call SPARSE_GRID_CC to create a Clenshaw Curtis grid.
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
            //    Input, int LEVEL_MAX, the level.
            //
        {
            int dim;
            double[] grid_point;
            double[] grid_weight;
            int point;
            int point_num;

            Console.WriteLine("");
            Console.WriteLine("TEST03:");
            Console.WriteLine("  SPARSE_GRID_CC makes a sparse Clenshaw Curtis grid.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
            //
            //  Determine the number of points.
            //
            point_num = Grid_ClenshawCurtis.sparse_grid_cfn_size(dim_num, level_max);

            Console.WriteLine("");
            Console.WriteLine("  Number of unique points in the grid = " + point_num + "");
            //
            //  Allocate space for the weights and points.
            //
            grid_weight = new double[point_num];
            grid_point = new double[dim_num * point_num];
            //
            //  Compute the weights and points.
            //
            Grid_ClenshawCurtis.sparse_grid_cc(dim_num, level_max, point_num, ref grid_weight, ref grid_point);
            //
            //  Print them out.
            //
            Console.WriteLine("");
            Console.WriteLine("  Grid weights:");
            Console.WriteLine("");
            for (point = 0; point < point_num; point++)
            {
                Console.WriteLine("  " + point.ToString().PadLeft(4)
                                       + " " + grid_weight[point].ToString().PadLeft(10) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Grid points:");
            Console.WriteLine("");
            for (point = 0; point < point_num; point++)
            {
                string cout = "  " + point.ToString().PadLeft(4);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += " " + grid_point[dim + point * dim_num].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

        }

        static void test04(int dim_num, int level_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 sums the weights and compares them to 2^DIM_NUM.
            //
            //  Discussion:
            //
            //    This routine gets the sparse grid indices and determines the 
            //    corresponding sparse grid abscissas.
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
            //    Input, int LEVEL_MAX, the level.
            //
        {
            double[] grid_point;
            double[] grid_weight;
            int point;
            int point_num;
            double weight_sum;
            double weight_sum_error;
            double weight_sum_exact;

            Console.WriteLine("");
            Console.WriteLine("TEST04:");
            Console.WriteLine("  Compute the weights of a Clenshaw Curtis sparse grid .");
            Console.WriteLine("");
            Console.WriteLine("  As a simple test, sum these weights.");
            Console.WriteLine("  They should sum to exactly 2^DIM_NUM.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
            //
            //  Determine the number of points.
            //
            point_num = Grid_ClenshawCurtis.sparse_grid_cfn_size(dim_num, level_max);

            Console.WriteLine("");
            Console.WriteLine("  Number of unique points in the grid = " + point_num + "");
            //
            //  Allocate space for the weights and points.
            //
            grid_weight = new double[point_num];
            grid_point = new double[dim_num * point_num];
            //
            //  Compute the weights and points.
            //
            Grid_ClenshawCurtis.sparse_grid_cc(dim_num, level_max, point_num, ref grid_weight, ref grid_point);
            //
            //  Sum the weights.
            //
            weight_sum = 0.0;
            for (point = 0; point < point_num; point++)
            {
                weight_sum = weight_sum + grid_weight[point];
            }

            weight_sum_exact = Math.Pow(2.0, dim_num);

            weight_sum_error = Math.Abs(weight_sum - weight_sum_exact);

            Console.WriteLine("");
            Console.WriteLine("    Weight sum     Exact sum    Difference");
            Console.WriteLine("");
            Console.WriteLine("  " + weight_sum.ToString().PadLeft(12)
                                   + "  " + weight_sum_exact.ToString().PadLeft(12)
                                   + "  " + weight_sum_error.ToString().PadLeft(12) + "");

        }

        static void test05(int dim_num, int level_max, int degree_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 tests a Clenshaw Curtis sparse grid rule for monomial exactness.
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
            //    02 July 2008
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
            int[] expon;
            double[] grid_point;
            double[] grid_weight;
            int h = 0;
            bool more;
            int point;
            int point_num;
            double quad_error;
            int t = 0;
            double volume;

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  Check the exactness of a Clenshaw Curtis sparse grid quadrature rule,");
            Console.WriteLine("  applied to all monomials of orders 0 to DEGREE_MAX.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
            Console.WriteLine("");
            Console.WriteLine("  The maximum total degree to be checked is DEGREE_MAX = " + degree_max + "");
            Console.WriteLine("");
            Console.WriteLine("  We expect this rule to be accurate up to and including total degree "
                              + 2 * level_max + 1 + "");
            //
            //  Determine the number of points in the rule.
            //
            point_num = Grid_ClenshawCurtis.sparse_grid_cfn_size(dim_num, level_max);

            Console.WriteLine("");
            Console.WriteLine("  Number of unique points in the grid = " + point_num + "");
            //
            //  Allocate space for the weights and points.
            //
            grid_weight = new double[point_num];
            grid_point = new double[dim_num * point_num];
            //
            //  Compute the weights and points.
            //
            Grid_ClenshawCurtis.sparse_grid_cc(dim_num, level_max, point_num, ref grid_weight,
                ref grid_point);
            //
            //  Rescale the weights, and translate the abscissas.
            //
            volume = Math.Pow(2.0, dim_num);

            for (point = 0; point < point_num; point++)
            {
                grid_weight[point] = grid_weight[point] / volume;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                for (point = 0; point < point_num; point++)
                {
                    grid_point[dim + point * dim_num] = (grid_point[dim + point * dim_num] + 1.0)
                                                        / 2.0;
                }
            }

            //
            //  Explore the monomials.
            //
            expon = new int[dim_num];

            Console.WriteLine("");
            Console.WriteLine("      Error      Total   Monomial");
            Console.WriteLine("                 Degree  Exponents");
            Console.WriteLine("");

            for (degree = 0; degree <= degree_max; degree++)
            {
                more = false;

                for (;;)
                {
                    Comp.comp_next(degree, dim_num, ref expon, ref more, ref h, ref t);

                    quad_error = MonomialQuadrature.monomial_quadrature(dim_num, expon, point_num,
                        grid_weight, grid_point);

                    string cout = "  " + quad_error.ToString().PadLeft(12)
                                       + "     " + degree.ToString().PadLeft(2)
                                       + "      ";

                    for (dim = 0; dim < dim_num; dim++)
                    {
                        cout += expon[dim].ToString().PadLeft(2);
                    }

                    Console.WriteLine(cout);

                    if (!more)
                    {
                        break;
                    }
                }

                Console.WriteLine("");
            }

        }

        static void test06(int dim_num, int level_max)

            //***************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 creates a sparse Clenshaw-Curtis grid and writes it to a file.
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
            int point_num;
            double[] r;
            string r_filename;
            double[] w;
            string w_filename;
            double[] x;
            string x_filename;

            Console.WriteLine("");
            Console.WriteLine("TEST06:");
            Console.WriteLine("  Call SPARSE_GRID_CC to make a sparse Clenshaw-Curtis grid.");
            Console.WriteLine("  Write the data to a set of quadrature files.");

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
            //
            //  Determine the number of points.
            //
            point_num = Grid_ClenshawCurtis.sparse_grid_cfn_size(dim_num, level_max);
            //
            //  Allocate space for the weights and points.
            //
            r = new double[dim_num * 2];
            w = new double[point_num];
            x = new double[dim_num * point_num];
            //
            //  Compute the weights and points.
            //
            for (dim = 0; dim < dim_num; dim++)
            {
                r[dim + 0 * dim_num] = -1.0;
                r[dim + 1 * dim_num] = +1.0;
            }

            Grid_ClenshawCurtis.sparse_grid_cc(dim_num, level_max, point_num, ref w, ref x);
            //
            //  Write the data out.
            //
            r_filename = "cc_d" + dim_num + "_level"
                         + level_max + "_r.txt";
            w_filename = "cc_d" + dim_num + "_level"
                         + level_max + "_w.txt";
            x_filename = "cc_d" + dim_num + "_level"
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
}