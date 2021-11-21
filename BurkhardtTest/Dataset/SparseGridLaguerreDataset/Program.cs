using System;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SparseGridLaguerreDataset;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPARSE_GRID_LAGUERRE_DATASET.
        //
        //  Discussion:
        //
        //    This program computes a sparse grid quadrature rule based on a 1D
        //    Gauss-Laguerre rule and writes it to a file.. 
        //
        //    The user specifies:
        //    * the spatial dimension of the quadrature region,
        //    * the level that defines the Smolyak grid.
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
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
    {
        int dim;
        int dim_num;
        int level_max;
        int level_min;
        int point;
        int point_num;
        double[] r;
        string r_filename;
        double[] w;
        string w_filename;
        double weight_sum;
        double[] x;
        string x_filename;

        Console.WriteLine("");
        Console.WriteLine("SPARSE_GRID_LAGUERRE_DATASET");
        Console.WriteLine("");
        Console.WriteLine("  Compute the abscissas and weights of a quadrature rule");
        Console.WriteLine("  associated with a sparse grid derived from a Smolyak");
        Console.WriteLine("  construction based on a 1D Gauss-Laguerre rule.");
        Console.WriteLine("");
        Console.WriteLine("  Inputs to the program include:");
        Console.WriteLine("");
        Console.WriteLine("    DIM_NUM, the spatial dimension.");
        Console.WriteLine("    (typically in the range of 2 to 10)");
        Console.WriteLine("");
        Console.WriteLine("    LEVEL_MAX, the level of the sparse grid.");
        Console.WriteLine("    (typically in the range of 0, 1, 2, 3, ...");
        Console.WriteLine("");
        Console.WriteLine("  Output from the program includes:");
        Console.WriteLine("");
        Console.WriteLine("    * A printed table of the abscissas and weights.");
        Console.WriteLine("");
        Console.WriteLine("    * A set of 3 files that define the quadrature rule.");
        Console.WriteLine("");
        Console.WriteLine("    (1) lag_d?_level?_r.txt, the ranges;");
        Console.WriteLine("    (2) lag_d?_level?_w.txt, the weights;");
        Console.WriteLine("    (3) lag_d?_level?_x.txt, the abscissas.");
        //
        //  Get the spatial dimension.
        //
        try
        {
            dim_num = Convert.ToInt32(args[0]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of DIM_NUM (1 or greater)");
            dim_num = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension requested is = " + dim_num + "");
        //
        //  Get the level.
        //
        try
        {
            level_max = Convert.ToInt32(args[1]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of LEVEL_MAX (0 or greater).");
            level_max = Convert.ToInt32(Console.ReadLine());
            ;
        }

        level_min = Math.Max(0, level_max + 1 - dim_num);

        Console.WriteLine("");
        Console.WriteLine("  LEVEL_MIN is = " + level_min + "");
        Console.WriteLine("  LEVEL_MAX is = " + level_max + "");
        // 
        //  How many distinct points will there be?
        //
        point_num = Grid_Laguerre.sparse_grid_laguerre_size(dim_num, level_max);

        Console.WriteLine("");
        Console.WriteLine("  The number of distinct abscissas in the");
        Console.WriteLine("  quadrature rule is determined from the spatial");
        Console.WriteLine("  dimension DIM_NUM and the level LEVEL_MAX.");
        Console.WriteLine("  For the given input, this value will be = " + point_num + "");
        //
        //  Allocate memory.
        //
        r = new double[dim_num * 2];
        w = new double[point_num];
        x = new double[dim_num * point_num];
        //
        //  Compute the weights and points.
        //
        for (dim = 0; dim < dim_num; dim++)
        {
            r[dim + 0 * dim_num] = 0.0E+00;
            r[dim + 1 * dim_num] = +typeMethods.r8_huge();
        }

        Grid_Laguerre.sparse_grid_laguerre(dim_num, level_max, point_num, ref w, ref x);

        typeMethods.r8mat_transpose_print_some(dim_num, point_num, x, 1, 1, dim_num,
            10, "  First 10 grid points:");

        typeMethods.r8vec_print_some(point_num, w, 1, 10, "  First 10 grid weights:");

        weight_sum = 0.0;
        for (point = 0; point < point_num; point++)
        {
            weight_sum += w[point];
        }

        Console.WriteLine("");
        Console.WriteLine("  Weights sum to   " + weight_sum + "");
        Console.WriteLine("  Correct value is " + 1.0 + "");
        //
        //  Construct appropriate file names.
        //
        r_filename = "lag_d" + dim_num
                             + "_level" + level_max + "_r.txt";
        w_filename = "lag_d" + dim_num
                             + "_level" + level_max + "_w.txt";
        x_filename = "lag_d" + dim_num
                             + "_level" + level_max + "_x.txt";
        //
        //  Write the rule to files.
        //
        Console.WriteLine("");
        Console.WriteLine("  Creating R file = \"" + r_filename + "\".");

        typeMethods.r8mat_write(r_filename, dim_num, 2, r);

        Console.WriteLine("  Creating W file = \"" + w_filename + "\".");

        typeMethods.r8mat_write(w_filename, 1, point_num, w);

        Console.WriteLine("  Creating X file = \"" + x_filename + "\".");

        typeMethods.r8mat_write(x_filename, dim_num, point_num, x);

        Console.WriteLine("");
        Console.WriteLine("SPARSE_GRID_LAGUERRE_DATASET");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}