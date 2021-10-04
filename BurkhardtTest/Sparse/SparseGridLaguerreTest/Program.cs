using System;
using Burkardt.Composition;
using Burkardt.Sparse;
using Burkardt.TriangleNS;
using Burkardt.Types;

namespace SparseGridLaguerreTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SPARSE_GRID_LAGUERRE_TEST.
            //
            //  Discussion:
            //
            //    SPARSE_GRID_LAGUERRE_TEST tests the SPARSE_GRID_LAGUERRE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 October 2007
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
            int degree_max;
            int dim_max;
            int dim_min;
            int dim_num;
            int level_max;
            int level_max_max;
            int level_max_min;

            Console.WriteLine("");
            Console.WriteLine("SPARSE_GRID_LAGUERRE_TEST");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Test the SPARSE_GRID_LAGUERRE library.");
            //
            //  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
            //
            dim_min = 1;
            dim_max = 5;
            level_max_min = 0;
            level_max_max = 10;

            test01(dim_min, dim_max, level_max_min, level_max_max);

            dim_min = 6;
            dim_max = 10;
            level_max_min = 0;
            level_max_max = 10;

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
            dim_num = 2;
            level_max = 3;
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
            //  Compute sparse Gauss-Laguerre rule for DIMENSION, LEVEL_MAX.
            //
            dim_num = 2;
            level_max = 0;
            test03(dim_num, level_max);

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
            degree_max = 3;
            test05(dim_num, level_max, degree_max);

            dim_num = 2;
            level_max = 1;
            degree_max = 5;
            test05(dim_num, level_max, degree_max);

            dim_num = 2;
            level_max = 2;
            degree_max = 7;
            test05(dim_num, level_max, degree_max);

            dim_num = 2;
            level_max = 3;
            degree_max = 9;
            test05(dim_num, level_max, degree_max);

            dim_num = 2;
            level_max = 4;
            degree_max = 11;
            test05(dim_num, level_max, degree_max);

            dim_num = 2;
            level_max = 5;
            degree_max = 13;
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
            degree_max = 8;
            test05(dim_num, level_max, degree_max);
            //
            //  Show how to write a rule to a file.
            //
            dim_num = 2;
            level_max = 3;

            test06(dim_num, level_max);
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("SPARSE_GRID_LAGUERRE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01(int dim_min, int dim_max, int level_max_min, int level_max_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests SPARSE_GRID_LAGUERRE_SIZE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 December 2009
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
            Console.WriteLine("  SPARSE_GRID_LAGUERRE_SIZE returns the number of distinct");
            Console.WriteLine("  points in a Gauss-Laguerre sparse grid.");
            Console.WriteLine("");
            Console.WriteLine("  Note that, unlike most sparse grids, a sparse grid based on");
            Console.WriteLine("  Gauss-Laguerre points is NOT nested.");
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
                cout += "  " + dim_num.ToString().PadLeft(10);
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
                    point_num = Grid_Laguerre.sparse_grid_laguerre_size(dim_num, level_max);
                    cout += "  " + point_num.ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

        }

        static void test02(int dim_num, int level_max)

            //***************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests SPARSE_GRID_LAGUERRE_INDEX.
            //
            //  Discussion:
            //
            //    The routine computes abstract indices that describe the sparse grid
            //    of Gauss-Laguerre points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 October 2007
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
            int[] grid_base;
            int[] grid_index;
            int level_min;
            int point;
            int point_num;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  SPARSE_GRID_LAGUERRE_INDEX returns abstract indices for the");
            Console.WriteLine("  points that make up a Gauss-Laguerre sparse grid.");

            level_min = Math.Max(0, level_max + 1 - dim_num);

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MIN = " + level_min + "");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");

            point_num = Grid_Laguerre.sparse_grid_laguerre_size(dim_num, level_max);

            Console.WriteLine("");
            Console.WriteLine("  Number of unique points in the grid = " + point_num + "");
            //
            //  Compute the orders and points.
            //
            grid_index = new int[dim_num * point_num];
            grid_base = new int[dim_num * point_num];

            Grid_Laguerre.sparse_grid_laguerre_index(dim_num, level_max, point_num, ref grid_index,
                ref grid_base);
            //
            //  Now we're done.  Print the merged grid data.
            //
            Console.WriteLine("");
            Console.WriteLine("  Grid index/base:");
            Console.WriteLine("");
            for (point = 0; point < point_num; point++)
            {
                string cout = "  " + point.ToString().PadLeft(4) + "  ";
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += grid_index[dim + point * dim_num].ToString().PadLeft(6);
                }

                Console.WriteLine(cout);
                cout = "        ";
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += grid_base[dim + point * dim_num].ToString().PadLeft(6);
                }

                Console.WriteLine(cout);
            }
        }

        static void test03(int dim_num, int level_max)

            //***************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 call SPARSE_GRID_LAGUERRE to create a sparse Gauss-Laguerre grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 October 2007
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
            int level_min;
            int point;
            int point_num;

            Console.WriteLine("");
            Console.WriteLine("TEST03:");
            Console.WriteLine("  SPARSE_GRID_LAGUERRE makes a sparse Gauss-Laguerre grid.");

            level_min = Math.Max(0, level_max + 1 - dim_num);

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MIN = " + level_min + "");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
            //
            //  Determine the number of points.
            //
            point_num = Grid_Laguerre.sparse_grid_laguerre_size(dim_num, level_max);

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
            Grid_Laguerre.sparse_grid_laguerre(dim_num, level_max, point_num, ref grid_weight,
                ref grid_point);
            //
            //  Print them out.
            //
            Console.WriteLine("");
            Console.WriteLine("  Grid weights:");
            Console.WriteLine("");
            string cout = "";
            for (point = 0; point < point_num; point++)
            {
                Console.WriteLine(point.ToString().PadLeft(4)
                                  + "  " + grid_weight[point].ToString("0.######").PadLeft(10) + "");
            }

            Console.WriteLine();
            Console.WriteLine("  Grid points:");
            Console.WriteLine("");
            for (point = 0; point < point_num; point++)
            {
                cout = "  " + point.ToString().PadLeft(4);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  "
                            + grid_point[dim + point * dim_num].ToString("0.######").PadLeft(9);
                }

                Console.WriteLine(cout);
            }
        }

        static void test04(int dim_num, int level_max)

            //***************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 sums the weights and compares them to 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 October 2007
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
            int level_min;
            int point;
            int point_num;
            double weight_sum;
            double weight_sum_error;
            double weight_sum_exact;

            Console.WriteLine("");
            Console.WriteLine("TEST04:");
            Console.WriteLine("  Compute the weights of a Gauss-Laguerre sparse grid .");
            Console.WriteLine("");
            Console.WriteLine("  As a simple test, sum these weights.");
            Console.WriteLine("  They should sum to exactly 1.");

            level_min = Math.Max(0, level_max + 1 - dim_num);

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MIN = " + level_min + "");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
            //
            //  Determine the number of points.
            //
            point_num = Grid_Laguerre.sparse_grid_laguerre_size(dim_num, level_max);

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
            Grid_Laguerre.sparse_grid_laguerre(dim_num, level_max, point_num, ref grid_weight,
                ref grid_point);
            //
            //  Sum the weights.
            //
            weight_sum = 0.0;
            for (point = 0; point < point_num; point++)
            {
                weight_sum = weight_sum + grid_weight[point];
            }

            weight_sum_exact = 1.0;

            weight_sum_error = Math.Abs(weight_sum - weight_sum_exact);

            Console.WriteLine("");
            Console.WriteLine("    Weight sum     Exact sum    Difference");
            Console.WriteLine("");
            Console.WriteLine("  " + weight_sum.ToString().PadLeft(14)
                                   + "  " + weight_sum_exact.ToString().PadLeft(14)
                                   + "  " + weight_sum_error.ToString().PadLeft(14) + "");

        }

        static void test05(int dim_num, int level_max, int degree_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 tests a Gauss-Laguerre sparse grid rule for monomial exactness.
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
            //    05 July 2008
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
            int h;
            int level_min;
            bool more;
            int point_num;
            double quad_error;
            int t;

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  Check the exactness of a Gauss-Laguerre sparse");
            Console.WriteLine("  grid quadrature rule, applied to all monomials ");
            Console.WriteLine("  of orders 0 to DEGREE_MAX.");

            level_min = Math.Max(0, level_max + 1 - dim_num);

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MIN = " + level_min + "");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
            Console.WriteLine("");
            Console.WriteLine("  The maximum total degree to be checked is DEGREE_MAX = " + degree_max + "");
            //
            //  Determine the number of points in the rule.
            //
            point_num = Grid_Laguerre.sparse_grid_laguerre_size(dim_num, level_max);

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
            Grid_Laguerre.sparse_grid_laguerre(dim_num, level_max, point_num, ref grid_weight,
                ref grid_point);
            //
            //  Explore the monomials.
            //
            expon = new int[dim_num];

            Console.WriteLine("");
            Console.WriteLine("      Error      Total   Monomial");
            Console.WriteLine("                 Degree  Exponents");

            for (degree = 0; degree <= degree_max; degree++)
            {
                more = false;
                h = 0;
                t = 0;

                Console.WriteLine("");
                for (;;)
                {
                    Comp.comp_next(degree, dim_num, ref expon, ref more, ref h, ref t);

                    quad_error = Burkardt.Quadrature.MonomialQuadrature.monomial_quadrature(dim_num, expon, point_num,
                        grid_weight, grid_point);

                    string cout = "  " + quad_error.ToString("0.#").PadLeft(12)
                                       + "     " + degree.ToString().PadLeft(2)
                                       + "      ";

                    for (dim = 0; dim < dim_num; dim++)
                    {
                        cout += expon[dim].ToString().PadLeft(3);
                    }

                    Console.WriteLine(cout);

                    if (!more)
                    {
                        break;
                    }
                }
            }
        }

        static void test06(int dim_num, int level_max)

            //***************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 creates a sparse Gauss-Laguerre grid and writes it to a file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 July 2009
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
            int level_min;
            int point_num;
            double[] r;
            string r_filename;
            double[] w;
            string w_filename;
            double[] x;
            string x_filename;

            Console.WriteLine("");
            Console.WriteLine("TEST06:");
            Console.WriteLine("  Call SPARSE_GRID_LAGUERRE to make a sparse Gauss-Laguerre grid.");
            Console.WriteLine("  Write the data to a set of quadrature files.");

            level_min = Math.Max(0, level_max + 1 - dim_num);

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_MIN = " + level_min + "");
            Console.WriteLine("  LEVEL_MAX = " + level_max + "");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
            //
            //  Determine the number of points.
            //
            point_num = Grid_Laguerre.sparse_grid_laguerre_size(dim_num, level_max);
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
                r[dim + 0 * dim_num] = 0.0;
                r[dim + 1 * dim_num] = +typeMethods.r8_huge();
            }

            Grid_Laguerre.sparse_grid_laguerre(dim_num, level_max, point_num, ref w, ref x);
            //
            //  Write the data out.
            //
            r_filename = "lg_d" + dim_num
                                + "_level" + level_max + "_r.txt";
            w_filename = "lg_d" + dim_num
                                + "_level" + level_max + "_w.txt";
            x_filename = "lg_d" + dim_num
                                + "_level" + level_max + "_x.txt";

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