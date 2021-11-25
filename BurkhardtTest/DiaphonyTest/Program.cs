﻿using System;
using System.Globalization;
using Burkardt;
using Burkardt.Table;
using Burkardt.Types;

namespace DiaphonyTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for DIAPHONY.
        //
        //  Discussion:
        //
        //    DIAPHONY reads a table dataset and applies the diaphony test.
        //
        //  Usage:
        //
        //    diaphony table_file
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string input_filename;

        //
        //  Get the input filename
        //
        try
        {
            input_filename = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the input filename:");
            input_filename = Console.ReadLine();
        }

        //
        //  Get the data size.
        //
        TableHeader h = typeMethods.r8mat_header_read(input_filename);
        int dim_num = h.m;
        int point_num = h.n;

        //
        //  Read the data.
        //
        double[] points = typeMethods.r8mat_data_read(input_filename, dim_num, point_num);

        if (typeMethods.r8mat_min(dim_num, point_num, points) < 0.0)
        {
            Console.WriteLine("");
            Console.WriteLine("DIAPHONY - Fatal error!");
            Console.WriteLine("  At least one coordinate of a point is less than 0.");
            return;
        }

        if (1.0 < typeMethods.r8mat_max(dim_num, point_num, points))
        {
            Console.WriteLine("");
            Console.WriteLine("DIAPHONY - Fatal error!");
            Console.WriteLine("  At least one coordinate of a point is greater than 1.");
            return;
        }

        //
        //  Analyze the data.
        //
        double d = Diaphony.diaphony_compute(dim_num, point_num, points);

        double e = 1.0 / Math.Sqrt(point_num);
        double de = d / e;

        Console.WriteLine("");
        Console.WriteLine("  File  M  N  Diaphony  1/sqrt(N)  D/sqrt(N)");
        Console.WriteLine("  " + input_filename
                               + "  " + dim_num.ToString().PadLeft(2)
                               + "  " + point_num.ToString().PadLeft(5)
                               + "  " + d.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + de.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }
}