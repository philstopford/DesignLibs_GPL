using System;
using Burkardt.Graph;
using Burkardt.Types;

namespace SubsetTestNS;

public static class DigraphTest
{
    public static void digraph_arc_euler_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIGRAPH_ARC_EULER_TEST calls DIGRAPH_ARC_EULER.
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
        const int NEDGE = 7;
        const int NNODE = 5;

        int[] inode = { 2, 1, 2, 1, 3, 5, 4 };
        int[] jnode = { 5, 4, 3, 2, 1, 1, 2 };
        bool success = false;
        int[] trail = new int[NEDGE];

        Console.WriteLine("");
        Console.WriteLine("DIGRAPH_ARC_EULER_TEST");
        Console.WriteLine("  DIGRAPH_ARC_EULER finds an Euler circuit of a digraph.");

        Digraph.digraph_arc_print(NEDGE, inode, jnode, "  The arc list of the digraph:");

        Digraph.digraph_arc_euler(NNODE, NEDGE, inode, jnode, ref success, ref trail);

        switch (success)
        {
            case true:
            {
                typeMethods.i4vec1_print(NEDGE, trail, "  The edge list of the Euler circuit:");

                Console.WriteLine("");
                Console.WriteLine("  The node list of the Euler circuit:");
                Console.WriteLine("");
                Console.WriteLine("	 I  Edge  Node");
                Console.WriteLine("");

                int i;
                for (i = 0; i < NEDGE; i++)
                {
                    int j = trail[i];

                    int jp1 = i + 1 == NEDGE ? trail[0] : trail[i + 1];

                    int in_;
                    if (jnode[j - 1] == inode[jp1 - 1])
                    {
                        in_ = jnode[j - 1];
                    }
                    else
                    {
                        Console.WriteLine("");
                        Console.WriteLine("The circuit has failed!");
                        Console.WriteLine("  JNODE[" + (j - 1) + "] = " + jnode[j - 1] + "");
                        Console.WriteLine("  INODE[" + (jp1 - 1) + "] = " + inode[jp1 - 1] + "");
                        break;
                    }

                    Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                                              + j.ToString().PadLeft(6) + "  "
                                                              + in_.ToString().PadLeft(6) + "");
                }

                break;
            }
            default:
                Console.WriteLine("");
                Console.WriteLine("  The digraph is not eulerian.");
                Console.WriteLine("");
                break;
        }
    }

    public static void digraph_arc_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIGRAPH_ARC_PRINT_TEST calls DIGRAPH_ARC_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] inode = { 2, 1, 2, 1, 3, 5, 4 };
        int[] jnode = { 5, 4, 3, 2, 1, 1, 2 };

        Console.WriteLine("");
        Console.WriteLine("DIGRAPH_ARC_PRINT_TEST");
        Console.WriteLine("  DIGRAPH_ARC_PRINT prints a digraph.");

        int nedge = 7;

        Digraph.digraph_arc_print(nedge, inode, jnode, "  The arc list of the digraph:");
    }

}