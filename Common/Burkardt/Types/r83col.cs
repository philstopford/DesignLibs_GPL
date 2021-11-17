using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r83col_print_part(int n, double[] a, int max_print, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83COL_PRINT_PART prints "part" of an R83COL.
        //
        //  Discussion:
        //
        //    An R83COL is an (N,3) array of R8's.
        //
        //    The user specifies MAX_PRINT, the maximum number of lines to print.
        //
        //    If N, the size of the vector, is no more than MAX_PRINT, then
        //    the entire vector is printed, one entry per line.
        //
        //    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
        //    followed by a line of periods suggesting an omission,
        //    and the last entry.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries of the vector.
        //
        //    Input, double A[N*3], the vector to be printed.
        //
        //    Input, int MAX_PRINT, the maximum number of lines
        //    to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        switch (max_print)
        {
            case <= 0:
                return;
        }

        switch (n)
        {
            case <= 0:
                return;
        }

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");

        if (n <= max_print)
        {
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + a[i + 0 * n].ToString().PadLeft(14)
                                       + "  " + a[i + 1 * n].ToString().PadLeft(14)
                                       + "  " + a[i + 2 * n].ToString().PadLeft(14) + "");
            }
        }
        else
        {
            switch (max_print)
            {
                case >= 3:
                {
                    for (i = 0; i < max_print - 2; i++)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + ": " + a[i + 0 * n].ToString().PadLeft(14)
                                               + "  " + a[i + 1 * n].ToString().PadLeft(14)
                                               + "  " + a[i + 2 * n].ToString().PadLeft(14) + "");
                    }

                    Console.WriteLine("  ........  ..............  ..............  ..............");
                    i = n - 1;
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + ": " + a[i + 0 * n].ToString().PadLeft(14)
                                           + "  " + a[i + 1 * n].ToString().PadLeft(14)
                                           + "  " + a[i + 2 * n].ToString().PadLeft(14) + "");
                    break;
                }
                default:
                {
                    for (i = 0; i < max_print - 1; i++)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + ": " + a[i + 0 * n].ToString().PadLeft(14)
                                               + "  " + a[i + 1 * n].ToString().PadLeft(14)
                                               + "  " + a[i + 2 * n].ToString().PadLeft(14) + "");
                    }

                    i = max_print - 1;
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + ": " + a[i + 0 * n].ToString().PadLeft(14)
                                           + "  " + a[i + 1 * n].ToString().PadLeft(14)
                                           + "  " + a[i + 2 * n].ToString().PadLeft(14)
                                           + "  " + "...more entries...");
                    break;
                }
            }
        }
    }

}