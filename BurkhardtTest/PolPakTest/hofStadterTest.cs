using System;
using Burkardt.Sequence;

namespace PolPakTest;

public static class hofStadterTest
{
    public static void f_hofstadter_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_HOFSTADTER_TEST tests F_HOFSTADTER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("F_HOFSTADTER_TEST");
        Console.WriteLine("  F_HOFSTADTER evaluates Hofstadter's recursive");
        Console.WriteLine("  F function.");
        Console.WriteLine("");
        Console.WriteLine("     N   F(N)");
        Console.WriteLine("");

        for (i = 0; i <= 30; i++)
        {
            int f = Hofstadter.f_hofstadter(i);

            Console.WriteLine("  "
                              + i.ToString().PadLeft(6) + "  "
                              + f.ToString().PadLeft(6) + "");
        }

    }

    public static void g_hofstadter_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    G_HOFSTADTER_TEST tests G_HOFSTADTER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("G_HOFSTADTER_TEST");
        Console.WriteLine("  G_HOFSTADTER evaluates Hofstadter's recursive");
        Console.WriteLine("  G function.");
        Console.WriteLine("");
        Console.WriteLine("     N   G(N)");
        Console.WriteLine("");

        for (i = 0; i <= 30; i++)
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(6) + "  "
                              + Hofstadter.g_hofstadter(i).ToString().PadLeft(6) + "");
        }

    }

    public static void h_hofstadter_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    H_HOFSTADTER_TEST tests H_HOFSTADTER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("H_HOFSTADTER_TEST");
        Console.WriteLine("  H_HOFSTADTER evaluates Hofstadter's recursive");
        Console.WriteLine("  H function.");

        Console.WriteLine("");
        Console.WriteLine("     N   H(N)");
        Console.WriteLine("");

        for (i = 0; i <= 30; i++)
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(6) + "  "
                              + Hofstadter.h_hofstadter(i).ToString().PadLeft(6) + "");
        }

    }

    public static void v_hofstadter_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_HOFSTADTER_TEST tests V_HOFSTADTER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("V_HOFSTADTER_TEST");
        Console.WriteLine("  V_HOFSTADTER evaluates Hofstadter's recursive");
        Console.WriteLine("  V function.");
        Console.WriteLine("");
        Console.WriteLine("     N   V(N)");
        Console.WriteLine("");

        for (i = 0; i <= 30; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(6)
                                   + "  " + Hofstadter.v_hofstadter(i).ToString().PadLeft(6) + "");
        }

    }

}