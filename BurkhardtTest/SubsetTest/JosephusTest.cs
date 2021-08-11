using System;
using Burkardt;

namespace SubsetTest
{
    public static class JosephusTest
    {
        public static void josephus_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JOSEPHUS_TEST tests JOSEPHUS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int k;
            int m;
            int n;
            int x;

            Console.WriteLine("");
            Console.WriteLine("JOSEPHUS_TEST");
            Console.WriteLine("  JOSEPHUS solves Josephus problems.");
            Console.WriteLine("");
            Console.WriteLine("    N    M    K	 X");
            Console.WriteLine("");

            m = 3;
            n = 41;
            k = 41;
            x = Josephus.josephus(n, m, k);

            Console.WriteLine(n.ToString().PadLeft(5) + "  "
                                                      + m.ToString().PadLeft(5) + "  "
                                                      + k.ToString().PadLeft(5) + "  "
                                                      + x.ToString().PadLeft(5) + "");

            m = -38;
            n = 41;
            k = 41;
            x = Josephus.josephus(n, m, k);

            Console.WriteLine(n.ToString().PadLeft(5) + "  "
                                                      + m.ToString().PadLeft(5) + "  "
                                                      + k.ToString().PadLeft(5) + "  "
                                                      + x.ToString().PadLeft(5) + "");

            m = 3;
            n = 41;
            k = 40;
            x = Josephus.josephus(n, m, k);

            Console.WriteLine(n.ToString().PadLeft(5) + "  "
                                                      + m.ToString().PadLeft(5) + "  "
                                                      + k.ToString().PadLeft(5) + "  "
                                                      + x.ToString().PadLeft(5) + "");

            m = 2;
            n = 64;
            k = 64;
            x = Josephus.josephus(n, m, k);

            Console.WriteLine(n.ToString().PadLeft(5) + "  "
                                                      + m.ToString().PadLeft(5) + "  "
                                                      + k.ToString().PadLeft(5) + "  "
                                                      + x.ToString().PadLeft(5) + "");

            m = 2;
            n = 1000;
            k = 1000;
            x = Josephus.josephus(n, m, k);

            Console.WriteLine(n.ToString().PadLeft(5) + "  "
                                                      + m.ToString().PadLeft(5) + "  "
                                                      + k.ToString().PadLeft(5) + "  "
                                                      + x.ToString().PadLeft(5) + "");
        }

    }
}