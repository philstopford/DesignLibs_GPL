﻿using System;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SparseGridGLDataset;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPARSE_GRID_GL_DATASET.
        //
        //  Discussion:
        //
        //    This program computes a sparse grid quadrature rule based on a 1D
        //    Gauss-Legendre rule and writes it to a file.. 
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
        //    09 July 2009
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
        int point;

        Console.WriteLine("");
        Console.WriteLine("SPARSE_GRID_GL_DATASET");
        Console.WriteLine("");
        Console.WriteLine("  Compute the abscissas and weights of a quadrature rule");
        Console.WriteLine("  associated with a sparse grid derived from a Smolyak");
        Console.WriteLine("  construction based on a 1D Gauss-Legendre rule.");
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
        Console.WriteLine("    (1) gl_d?_level?_r.txt, the ranges;");
        Console.WriteLine("    (2) gl_d?_level?_w.txt, the weights;");
        Console.WriteLine("    (3) gl_d?_level?_x.txt, the abscissas.");
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
        }

        int level_min = Math.Max(0, level_max + 1 - dim_num);

        Console.WriteLine("");
        Console.WriteLine("  LEVEL_MIN is = " + level_min + "");
        Console.WriteLine("  LEVEL_MAX is = " + level_max + "");
        // 
        //  How many distinct points will there be?
        //
        int point_num = Grid_GaussLegendre.sparse_grid_gl_size(dim_num, level_max);

        Console.WriteLine("");
        Console.WriteLine("  The number of distinct abscissas in the");
        Console.WriteLine("  quadrature rule is determined from the spatial");
        Console.WriteLine("  dimension DIM_NUM and the level LEVEL_MAX.");
        Console.WriteLine("  For the given input, this value will be = " + point_num + "");
        //
        //  Allocate memory.
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

        typeMethods.r8mat_transpose_print_some(dim_num, point_num, x, 1, 1, dim_num,
            10, "  First 10 grid points:");

        typeMethods.r8vec_print_some(point_num, w, 1, 10, "  First 10 grid weights:");

        double weight_sum = 0.0;
        for (point = 0; point < point_num; point++)
        {
            weight_sum += w[point];
        }

        Console.WriteLine("");
        Console.WriteLine("  Weights sum to   " + weight_sum + "");
        Console.WriteLine("  Correct value is " + Math.Pow(2.0, dim_num) + "");
        //
        //  Construct appropriate file names.
        //
        string r_filename = "gl_d" + dim_num
                                   + "_level" + level_max + "_r.txt";
        string w_filename = "gl_d" + dim_num
                                   + "_level" + level_max + "_w.txt";
        string x_filename = "gl_d" + dim_num
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
        Console.WriteLine("SPARSE_GRID_GL_DATASET");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}