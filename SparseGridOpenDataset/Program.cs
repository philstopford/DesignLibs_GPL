using System;
using Burkardt.Quadrature;
using Burkardt.Sparse;
using Burkardt.Types;

namespace SparseGridOpenDataset
{
 class Program
 {
  static void Main(string[] args)
   //****************************************************************************80
   //
   //  Purpose:
   //
   //    MAIN is the main program for SPARSE_GRID_OPEN_DATASET.
   //
   //  Discussion:
   //
   //    This program computes a quadrature rule and writes it to a file.
   //
   //    The quadrature rule is associated with a sparse grid derived from
   //    a Smolyak construction using an open 1D quadrature rule. 
   //
   //    The user specifies:
   //    * the spatial dimension of the quadrature region,
   //    * the level that defines the Smolyak grid.
   //    * the open 1D quadrature rule.
   //
   //  Licensing:
   //
   //    This code is distributed under the GNU LGPL license. 
   //
   //  Modified:
   //
   //    01 February 2009
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
   int[] grid_index;
   double[] grid_point;
   double[] grid_region;
   double[] grid_weight;
   double h;
   int level_max;
   int m;
   int n;
   int order_max;
   int point;
   int point_num;
   string r_filename = "";
   int rule;
   string w_filename = "";
   double weight_sum;
   string x_filename = "";

   Console.WriteLine("");
   Console.WriteLine("SPARSE_GRID_OPEN_DATASET");
   Console.WriteLine("");
   Console.WriteLine("  Compute the abscissas and weights of a quadrature rule");
   Console.WriteLine("  associated with a sparse grid derived from a Smolyak");
   Console.WriteLine("  construction based on an open quadrature rule.");
   Console.WriteLine("");
   Console.WriteLine("  Inputs to the program include:");
   Console.WriteLine("");
   Console.WriteLine("    DIM_NUM, the spatial dimension.");
   Console.WriteLine("    (typically in the range of 2 to 10)");
   Console.WriteLine("");
   Console.WriteLine("    LEVEL_MAX, the \"level\" of the sparse grid.");
   Console.WriteLine("    (typically in the range of 0, 1, 2, 3, ...");
   Console.WriteLine("");
   Console.WriteLine("    RULE, the 1D quadrature rule");
   Console.WriteLine("    2: Fejer Type 2 (\"F2\").");
   Console.WriteLine("    3: Gauss-Patterson (\"GP\");");
   Console.WriteLine("    4: Newton-Cotes Open (\"NCO\").");
   Console.WriteLine("    5: Tanh-Sinh (\"TS\").");
   Console.WriteLine("");
   Console.WriteLine("  Output from the program includes:");
   Console.WriteLine("");
   Console.WriteLine("    A printed table of the abscissas and weights.");
   Console.WriteLine("");
   Console.WriteLine("    A set of files defining the quadrature rules.");
   Console.WriteLine("");
   Console.WriteLine("    \"***_d?_level?_x.txt\", a file of the abscissas;");
   Console.WriteLine("    \"***_d?_level?_w.txt\", a file of the weights;");
   Console.WriteLine("    \"***_d?_level?_r.txt\", a file of the ranges.");
   //
   //  Get the spatial dimension:
   //
   try
   {
    dim_num = Convert.ToInt32(args[0]);
   }
   catch
   {
    Console.WriteLine("");
    Console.WriteLine("SPARSE_GRID_OPEN_DATASET:");
    Console.WriteLine("  Enter the value of DIM_NUM.");

    dim_num = Convert.ToInt32(Console.ReadLine());
   }

   Console.WriteLine("");
   Console.WriteLine("  Spatial dimension requested is = " + dim_num + "");
   //
   //  Get the product file root name:
   //
   try
   {
    level_max = Convert.ToInt32(args[1]);
   }
   catch
   {
    Console.WriteLine("");
    Console.WriteLine("SPARSE_GRID_OPEN_DATASET:");
    Console.WriteLine("  Enter the value of LEVEL_MAX.");
    level_max = Convert.ToInt32(Console.ReadLine());
   }

   Console.WriteLine("");
   Console.WriteLine("  The sparse grid level is = " + level_max + "");
   //
   //  Get the rule index:
   //
   try
   {
    rule = Convert.ToInt32(args[2]);
   }
   catch
   {
    Console.WriteLine("");
    Console.WriteLine("SPARSE_GRID_OPEN_DATASET:");
    Console.WriteLine("  Enter the value of RULE.");
    Console.WriteLine("  2 = F2   = Fejer Type 2 Rule,");
    Console.WriteLine("  3 = GP   = Gauss-Patterson,");
    Console.WriteLine("  4 = NCO  = Newton-Cotes Open,");
    Console.WriteLine("  5 = TS   = Tanh-Sinh.");
    rule = Convert.ToInt32(Console.ReadLine());
   }

   Console.WriteLine("");
   Console.WriteLine("  The 1D quadrature rule index = " + rule + "");

   if (rule == 2)
   {
    Console.WriteLine("  F2:   Fejer Type 2 Rule.");
   }
   else if (rule == 3)
   {
    Console.WriteLine("  GP:   Gauss-Patterson Rule.");
   }
   else if (rule == 4)
   {
    Console.WriteLine("  NCO:  Newton-Cotes Open Rule.");
   }
   else if (rule == 5)
   {
    Console.WriteLine("  TS:   Tanh-Sinh Rule.");
   }
   else
   {
    Console.WriteLine("");
    Console.WriteLine("SPARSE_GRID_OPEN_DATASET - Fatal error!");
    Console.WriteLine("  Illegal value of RULE.");
    return;
   }

   //
   //  How many distinct points will there be?
   //
   point_num = Grid.sparse_grid_ofn_size(dim_num, level_max);

   Console.WriteLine("");
   Console.WriteLine("  The number of distinct abscissas in the");
   Console.WriteLine("  quadrature rule is determined from the spatial");
   Console.WriteLine("  dimension DIM_NUM and the level LEVEL_MAX.");
   Console.WriteLine("  For the given input, this value will be = " + point_num + "");

   grid_point = new double[dim_num * point_num];
   //
   //  Determine the index vector, relative to the full product grid,
   //  that identifies the points in the sparse grid.
   //
   grid_index = Grid.spgrid_open_index(dim_num, level_max, point_num);

   typeMethods.i4mat_transpose_print_some(dim_num, point_num, grid_index, 1, 1,
    dim_num, 10, "  First 10 entries of grid index:");
   //
   //  Compute the physical coordinates of the abscissas.
   //
   order_max = (int) Math.Pow(2, level_max + 1) - 1;

   if (rule == 5)
   {
    m = level_max - 3;
    n = ((order_max + 1) / 2) - 1;
    h = 4.0 / (double) (order_max + 1);

    Console.WriteLine("  M = " + m
                               + "  ORDER_MAX = " + order_max
                               + "  N = " + n
                               + "  H = " + h + "");
   }

   if (rule == 2)
   {
    for (point = 0; point < point_num; point++)
    {
     for (dim = 0; dim < dim_num; dim++)
     {
      grid_point[dim + point * dim_num] =
       Fejer2.f2_abscissa(order_max, grid_index[dim + point * dim_num]);
     }
    }
   }
   else if (rule == 3)
   {
    for (point = 0; point < point_num; point++)
    {
     for (dim = 0; dim < dim_num; dim++)
     {
      grid_point[dim + point * dim_num] =
       PattersonQuadrature.gp_abscissa(order_max, grid_index[dim + point * dim_num]);
     }
    }
   }
   else if (rule == 4)
   {
    for (point = 0; point < point_num; point++)
    {
     for (dim = 0; dim < dim_num; dim++)
     {
      grid_point[dim + point * dim_num] =
       NewtonCotesQuadrature.nco_abscissa(order_max, grid_index[dim + point * dim_num]);
     }
    }
   }
   else if (rule == 5)
   {
    for (point = 0; point < point_num; point++)
    {
     for (dim = 0; dim < dim_num; dim++)
     {
      grid_point[dim + point * dim_num] =
       TanhSinh.ts_abscissa(order_max, grid_index[dim + point * dim_num]);
     }
    }
   }

   typeMethods.r8mat_transpose_print_some(dim_num, point_num, grid_point, 1, 1,
    dim_num, 10, "  First 10 entries of grid point:");
   //
   //  Gather the weights.
   //
   grid_weight = Grid.spgrid_open_weights(dim_num, level_max, point_num,
    grid_index, rule);

   typeMethods.r8vec_print_some(point_num, grid_weight, 1, 10,
    "  First 10 grid weights:");

   weight_sum = typeMethods.r8vec_sum(point_num, grid_weight);

   Console.WriteLine("");
   Console.WriteLine("  Weights sum to   "
                     + weight_sum.ToString("0.################").PadLeft(24) + "");
   Console.WriteLine("  Correct value is "
                     + Math.Pow(2.0, dim_num).ToString("0.################").PadLeft(24) + "");
   //
   //  Write the rule to files.
   //
   if (rule == 2)
   {
    r_filename = "f2_d" + dim_num
                        + "_level" + level_max + "_r.txt";
    w_filename = "f2_d" + dim_num
                        + "_level" + level_max + "_w.txt";
    x_filename = "f2_d" + dim_num
                        + "_level" + level_max + "_x.txt";
   }
   else if (rule == 3)
   {
    r_filename = "gp_d" + dim_num
                        + "_level" + level_max + "_r.txt";
    w_filename = "gp_d" + dim_num
                        + "_level" + level_max + "_w.txt";
    x_filename = "gp_d" + dim_num
                        + "_level" + level_max + "_x.txt";
   }
   else if (rule == 4)
   {
    r_filename = "nco_d" + dim_num
                         + "_level" + level_max + "_r.txt";
    w_filename = "nco_d" + dim_num
                         + "_level" + level_max + "_w.txt";
    x_filename = "nco_d" + dim_num
                         + "_level" + level_max + "_x.txt";
   }
   else if (rule == 5)
   {
    r_filename = "ts_d" + dim_num
                        + "_level" + level_max + "_r.txt";
    w_filename = "ts_d" + dim_num
                        + "_level" + level_max + "_w.txt";
    x_filename = "ts_d" + dim_num
                        + "_level" + level_max + "_x.txt";
   }

   Console.WriteLine("");
   Console.WriteLine("  Creating X file = \"" + x_filename + "\".");

   typeMethods.r8mat_write(x_filename, dim_num, point_num, grid_point);

   Console.WriteLine("  Creating W file = \"" + w_filename + "\".");

   typeMethods.r8mat_write(w_filename, 1, point_num, grid_weight);

   grid_region = new double[dim_num * 2];

   for (dim = 0; dim < dim_num; dim++)
   {
    grid_region[dim + 0 * dim_num] = -1.0;
    grid_region[dim + 1 * dim_num] = +1.0;
   }

   Console.WriteLine("  Creating R file = \"" + r_filename + "\".");

   typeMethods.r8mat_write(r_filename, dim_num, 2, grid_region);

   Console.WriteLine("");
   Console.WriteLine("SPARSE_GRID_OPEN_DATASET:");
   Console.WriteLine("  Normal end of execution.");
   Console.WriteLine("");
  }
 }
}