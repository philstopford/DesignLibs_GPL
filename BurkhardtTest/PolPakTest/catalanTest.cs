using System;
using Burkardt.Sequence;

namespace PolPakTest;

public static class catalanTest
{

    public static void catalan_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CATALAN_TEST tests CATALAN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c = 0;
        int[] c2;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("CATALAN_TEST");
        Console.WriteLine("  CATALAN computes Catalan numbers.");
        Console.WriteLine("");
        Console.WriteLine("  N  exact C(I)  computed C(I)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Catalan.catalan_values(ref n_data, ref n, ref c);

            if (n_data == 0)
            {
                break;
            }

            c2 = Catalan.catalan(n);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + c.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + c2[n].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

        }

    }

    public static void catalan_row_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CATALAN_ROW_NEXT_TEST tests CATALAN_ROW_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 10;

        int[] c = new int[N_MAX + 1];
        int i;
        int n;
        bool next;

        Console.WriteLine("");
        Console.WriteLine("CATALAN_ROW_NEXT_TEST");
        Console.WriteLine("  CATALAN_ROW_NEXT computes a row of Catalan''s triangle.");
        Console.WriteLine("");
        Console.WriteLine("  First, compute row 7:");
        Console.WriteLine("");

        next = false;
        n = 7;
        Catalan.catalan_row_next(next, n, ref c);

        string cout = n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
        for (i = 0; i <= n; i++)
        {
            cout += c[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  ";
        }

        Console.WriteLine(cout);

        Console.WriteLine("");
        Console.WriteLine("  Now compute rows consecutively, one at a time:");
        Console.WriteLine("");

        next = false;

        for (n = 0; n <= N_MAX; n++)
        {
            Catalan.catalan_row_next(next, n, ref c);
            next = true;

            cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
            for (i = 0; i <= n; i++)
            {
                cout += c[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  ";
            }

            Console.WriteLine(cout);

        }

    }

}