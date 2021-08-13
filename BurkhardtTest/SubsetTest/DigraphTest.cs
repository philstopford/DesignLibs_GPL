using System;
using Burkardt;
using Burkardt.Graph;
using Burkardt.Types;

namespace SubsetTestNS
{
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
            int NEDGE = 7;
            int NNODE = 5;

            int i;
            int in_;
            int[] inode = { 2, 1, 2, 1, 3, 5, 4 };
            int j;
            int[] jnode = { 5, 4, 3, 2, 1, 1, 2 };
            int jp1;
            bool success = false;
            int[] trail = new int[NEDGE];

            Console.WriteLine("");
            Console.WriteLine("DIGRAPH_ARC_EULER_TEST");
            Console.WriteLine("  DIGRAPH_ARC_EULER finds an Euler circuit of a digraph.");

            Digraph.digraph_arc_print(NEDGE, inode, jnode, "  The arc list of the digraph:");

            Digraph.digraph_arc_euler(NNODE, NEDGE, inode, jnode, ref success, ref trail);

            if (success)
            {
                typeMethods.i4vec1_print(NEDGE, trail, "  The edge list of the Euler circuit:");

                Console.WriteLine("");
                Console.WriteLine("  The node list of the Euler circuit:");
                Console.WriteLine("");
                Console.WriteLine("	 I  Edge  Node");
                Console.WriteLine("");

                for (i = 0; i < NEDGE; i++)
                {
                    j = trail[i];

                    if (i + 1 == NEDGE)
                    {
                        jp1 = trail[0];
                    }
                    else
                    {
                        jp1 = trail[i + 1];
                    }

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
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  The digraph is not eulerian.");
                Console.WriteLine("");
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
            int nedge;

            Console.WriteLine("");
            Console.WriteLine("DIGRAPH_ARC_PRINT_TEST");
            Console.WriteLine("  DIGRAPH_ARC_PRINT prints a digraph.");

            nedge = 7;

            Digraph.digraph_arc_print(nedge, inode, jnode, "  The arc list of the digraph:");
        }

    }
}