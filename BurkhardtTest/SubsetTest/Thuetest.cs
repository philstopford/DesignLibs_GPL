using System;
using Burkardt;
using Burkardt.Sequence;

namespace SubsetTestNS;

public static class Thuetest
{
    public static void thue_binary_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    THUE_BINARY_NEXT_TEST tests THUE_BINARY_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 100;

        int i;
        int j;
        int n;
        int[] thue = new int[N_MAX];

        Console.WriteLine("");
        Console.WriteLine("THUE_BINARY_NEXT_TEST");
        Console.WriteLine("  THUE_BINARY_NEXT returns the next");
        Console.WriteLine("  Thue binary sequence.");
        Console.WriteLine("");

        n = 1;
        thue[0] = 0;
        string cout = n.ToString().PadLeft(4) + "    ";
        for (i = 0; i < n; i++)
        {
            cout += thue[i];
        }

        Console.WriteLine(cout);

        for (i = 1; i <= 6; i++)
        {
            Thue.thue_binary_next(ref n, ref thue);

            cout = n.ToString().PadLeft(4) + "    ";
            for (j = 0; j < n; j++)
            {
                cout += thue[j];
            }

            Console.WriteLine(cout);
        }
    }

    public static void thue_ternary_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    THUE_TERNARY_NEXT_TEST tests THUE_TERNARY_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 100;

        int i;
        int j;
        int n;
        int[] thue = new int[N_MAX];

        Console.WriteLine("");
        Console.WriteLine("THUE_TERNARY_NEXT_TEST");
        Console.WriteLine("  THUE_TERNARY_NEXT returns the next");
        Console.WriteLine("  Thue ternary sequence.");
        Console.WriteLine("");

        n = 1;
        thue[0] = 1;
        string cout = n.ToString().PadLeft(4) + "    ";
        for (i = 0; i < n; i++)
        {
            cout += thue[i];
        }

        Console.WriteLine(cout);

        for (i = 1; i <= 5; i++)
        {
            Thue.thue_ternary_next(ref n, ref thue);

            cout = n.ToString().PadLeft(4) + "    ";
            for (j = 0; j < n; j++)
            {
                cout += thue[j];
            }

            Console.WriteLine(cout);
        }

    }





}