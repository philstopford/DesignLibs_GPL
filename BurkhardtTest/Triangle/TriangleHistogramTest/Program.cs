using System;
using Burkardt.Table;
using Burkardt.Types;

namespace TriangleHistogramTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE_HISTOGRAM.
        //
        //  Discussion:
        //
        //    TRIANGLE_HISTOGRAM reads a file of data, which is presumed to be a 
        //    list of points in the unit triangle.
        //
        //    The unit triangle has vertices (1,0), (0,1) and (0,0).
        //
        //    The program divides the unit triangle into N * ( N + 1 ) / 2
        //    subtriangles of equal area, and counts the number of data points
        //    that are contained in each.
        //
        //    A report is produced.
        //
        //    One use for this program is to examine the uniformity of distribution
        //    of a sampling of points in the triangle, produced by an algorithm
        //    that is attempting a uniform sampling.
        //
        //  Usage:
        //
        //    triangle_histogram data_filename n
        //
        //    where
        //
        //    data_filename is the name of the file containing the sample points,
        //    n is the number of subintervals each triangle side is divided into.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        string data_filename;
        int data_num;
        double[] data_xy;
        int dim_num;
        int[] histo;
        double histo_ave;
        int histo_max;
        int histo_min;
        double histo_std;
        int i;
        int j;
        int k;
        int n;
        int sub_num;
        int t;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_HISTOGRAM");
        Console.WriteLine("  Compute a histogram for data in the reference triangle.");
        Console.WriteLine("");
        //
        //  Argument 1 is the data filename.
        //
        try
        {
            data_filename = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_HISTOGRAM:");
            Console.WriteLine("  Please enter the data filename.");

            data_filename = Console.ReadLine();
        }

        //
        //  Argument 2 is the value of N.
        //
        try
        {
            n = Convert.ToInt32(args[1]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_HISTOGRAM:");
            Console.WriteLine("  Please enter the value of N.");

            n = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Read the DATA_XY data.
        //
        TableHeader h = typeMethods.r8mat_header_read(data_filename);
        dim_num = h.m;
        data_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + data_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + dim_num + "");
        Console.WriteLine("  Number of points DATA_NUM  = " + data_num + "");

        if (dim_num != 2)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_HISTOGRAM - Fatal error!");
            Console.WriteLine("  Dataset must have spatial dimension 2.");
            return;
        }

        data_xy = typeMethods.r8mat_data_read(data_filename, dim_num, data_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + data_filename + "\".");

        typeMethods.r8mat_transpose_print_some(dim_num, data_num, data_xy, 1, 1, 5, 5,
            "  5 by 5 portion of data read from file:");
        //
        //  Check the points for legality.
        //
        for (d = 0; d < data_num; d++)
        {
            if (data_xy[0 + d * 2] < 0.0 ||
                data_xy[1 + d * 2] < 0.0 ||
                1.0 < data_xy[0 + d * 2] + data_xy[1 + d * 2])
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_HISTOGRAM - Fatal error!");
                Console.WriteLine("  Point " + d + " does not lie inside the triangle.");
                Console.WriteLine("  X = " + data_xy[0 + d * 2] + "");
                Console.WriteLine("  Y = " + data_xy[1 + d * 2] + "");
                return;
            }
        }

        //
        //  Prepare the histogram.
        //
        sub_num = n * n;
        histo = new int[sub_num + 1];
        for (i = 0; i < sub_num + 1; i++)
        {
            histo[i] = 0;
        }

        //
        //  For each point, compute its barycentric coordinates, and hence
        //  its triangle.
        //
        for (d = 0; d < data_num; d++)
        {
            i = (int) (data_xy[0 + d * 2] * n) + 1;
            j = (int) (data_xy[1 + d * 2] * n) + 1;
            k = (int) ((1.0 - data_xy[0 + d * 2] - data_xy[1 + d * 2]) * n) + 1;

            if (i < 1 || n < i ||
                j < 1 || n < j ||
                k < 1 || n < k)
            {
                t = n * n + 1;
            }
            else
            {
                t = 2 * i + (n - j) * (n - j) + (k - 1) % 2 - 1;
            }

            //
            //  Increment the histogram vector.
            //
            histo[t - 1] += 1;
        }

        //
        //  Histogram statistics.
        //
        histo_ave = typeMethods.i4vec_mean(sub_num, histo);
        histo_max = typeMethods.i4vec_max(sub_num, histo);
        histo_min = typeMethods.i4vec_min(sub_num, histo);
        histo_std = typeMethods.i4vec_std(sub_num, histo);
        //
        //  Print the histogram.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data from file \"" + data_filename + "\".");
        Console.WriteLine("  Number of points = " + data_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Triangle refinement level N = " + n + "");
        Console.WriteLine("  Number of triangles = " + sub_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Counting number of points in each subtriangle:");
        Console.WriteLine("");
        Console.WriteLine("  Minimum = " + histo_min + "");
        Console.WriteLine("  Average = " + histo_ave + "");
        Console.WriteLine("  Maximum = " + histo_max + "");
        Console.WriteLine("  STD     = " + histo_std + "");
        Console.WriteLine("");
        Console.WriteLine("     COUNT");
        Console.WriteLine("");
        for (k = 0; k < n * n; k++)
        {
            Console.WriteLine("  " + histo[k].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of out-or-range points = " + histo[n * n] + "");

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_HISTOGRAM:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}