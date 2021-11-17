using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Burkardt.Types;

namespace Burkardt.Table;

public class TableHeader
{
    public int m { get; set; }
    public int n { get; set; }
    public int code { get; set; }

    public TableHeader()
    {
        code = 1;
    }
}

public static class TableSummary
{
    public static void r8block_print ( int l, int m, int n, double[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLOCK_PRINT prints a double precision block (a 3D matrix).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int L, M, N, the dimensions of the block.
        //
        //    Input, double A[L*M*N], the matrix to be printed.
        //
        //    Input, char *TITLE, a title to be printed first.
        //    TITLE may be blank.
        //
    {
        if ( 0 < typeMethods.s_len_trim ( title ) )
        {
            Console.WriteLine();
            Console.WriteLine(title);
        }

        for (int k = 1; k <= n; k++ )
        {
            Console.WriteLine();
            Console.WriteLine("  K = " + k);
            Console.WriteLine();
            for (int jlo = 1; jlo <= m; jlo += 5 )
            {
                int jhi = Math.Min( jlo + 4, m );
                Console.WriteLine();
                string cout = "      ";
                for (int j = jlo; j <= jhi; j++ )
                {
                    cout += j.ToString().PadLeft(7) + "       ";
                }
                Console.WriteLine(cout);
                Console.WriteLine();
                cout = "";
                for (int i = 1; i <= l; i++ )
                {
                    string t = i.ToString().PadLeft(4);
                    cout += "  " + t;
                    for (int j = jlo; j <= jhi; j++ )
                    {
                        t = a[i-1+(j-1)*l+(k-1)*l*m].ToString().PadLeft(12);
                        cout += "  " +t;
                    }
                    Console.WriteLine(cout);
                }
            }
        }
    }
}