using System;
using Burkardt;
using Burkardt.Types;

namespace DijkstraTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN runs an example of Dijkstra's minimum distance algorithm.
            //
            //  Discussion:
            //
            //    Given the distance matrix that defines a graph, we seek a list
            //    of the minimum distances between node 0 and all other nodes.
            //
            //    This program sets up a small example problem and solves it.
            //
            //    The correct minimum distances are:
            //
            //      0   35   15   45   49   41
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2010
            //
            //  Author:
            //
            //    Original C version by Norm Matloff, CS Dept, UC Davis.
            //    This C++ version by John Burkardt.
            //
        {
            int NV = 6;
            int i;
            int j;
            int[] mind;
            int[][] ohd = new int [NV][];

            Console.WriteLine("");
            Console.WriteLine("DIJKSTRA");
            Console.WriteLine("  Use Dijkstra's algorithm to determine the minimum");
            Console.WriteLine("  distance from node 0 to each node in a graph,");
            Console.WriteLine("  given the distances between each pair of nodes.");
            //
            //  Initialize the problem data.
            //
            init(ohd);
            //
            //  Print the distance matrix.
            //
            Console.WriteLine("");
            Console.WriteLine("  Distance matrix:");
            Console.WriteLine("");
            for (i = 0; i < NV; i++)
            {
                ohd[i] = new int[NV];
                string cout = "";
                for (j = 0; j < NV; j++)
                {
                    if (ohd[i][j] == typeMethods.i4_huge())
                    {
                        cout += "  Inf";
                    }
                    else
                    {
                        cout += "  " + ohd[i][j].ToString().PadLeft(3);
                    }
                }

                Console.WriteLine(cout);
            }

            //
            //  Carry out the algorithm.
            //
            mind = Dijkstra.dijkstra_distance(ohd);
            //
            //  Print the results.
            //
            Console.WriteLine("");
            Console.WriteLine("  Minimum distances from node 0:");
            Console.WriteLine("");
            for (i = 0; i < NV; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                    + "  " + mind[i].ToString().PadLeft(2) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("DIJKSTRA");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }

        static void init(int[][] ohd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INIT initializes the problem data.
            //
            //  Discussion:
            //
            //    The graph uses 6 nodes, and has the following diagram and
            //    distance matrix:
            //
            //    N0--15--N2-100--N3           0   40   15  Inf  Inf  Inf
            //      .      .     .            40    0   20   10   25    6
            //       .     .    .             15   20    0  100  Inf  Inf
            //        40  20  10             Inf   10  100    0  Inf  Inf
            //          .  .  .              Inf   25  Inf  Inf    0    8
            //           . . .               Inf    6  Inf  Inf    8    0
            //            N1
            //            . .
            //           .   .
            //          6    25
            //         .       .
            //        .         .
            //      N5----8-----N4
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2010
            //
            //  Author:
            //
            //    Original C version by Norm Matloff, CS Dept, UC Davis.
            //    This C++ version by John Burkardt.
            //
            //  Parameters:
            //
            //    Output, int OHD[NV][NV], the distance of the direct link between
            //    nodes I and J.
            //
        {
            for (int i = 0; i < ohd.Length; i++)
            {
                ohd[i] = new int[ohd.Length];
                for (int j = 0; j < ohd[i].Length; j++)
                {
                    if (i == j)
                    {
                        ohd[i][i] = 0;
                    }
                    else
                    {
                        ohd[i][j] = typeMethods.i4_huge();
                    }
                }
            }

            ohd[0][1] = ohd[1][0] = 40;
            ohd[0][2] = ohd[2][0] = 15;
            ohd[1][2] = ohd[2][1] = 20;
            ohd[1][3] = ohd[3][1] = 10;
            ohd[1][4] = ohd[4][1] = 25;
            ohd[2][3] = ohd[3][2] = 100;
            ohd[1][5] = ohd[5][1] = 6;
            ohd[4][5] = ohd[5][4] = 8;
        }
    }
}