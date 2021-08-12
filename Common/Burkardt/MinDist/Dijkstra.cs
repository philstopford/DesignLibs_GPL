using System;

namespace Burkardt.MinDist
{
    public static class Dijkstra
    {
        public static int[] dijkstra_distance(int[][] ohd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIJKSTRA_DISTANCE uses Dijkstra's minimum distance algorithm.
            //
            //  Discussion:
            //
            //    We essentially build a tree.  We start with only node 0 connected
            //    to the tree, and this is indicated by setting CONNECTED[0] = TRUE.
            //
            //    We initialize MIND[I] to the one step distance from node 0 to node I.
            //    
            //    Now we search among the unconnected nodes for the node MV whose minimum
            //    distance is smallest, and connect it to the tree.  For each remaining
            //    unconnected node I, we check to see whether the distance from 0 to MV
            //    to I is less than that recorded in MIND[I], and if so, we can reduce
            //    the distance.
            //
            //    After NV-1 steps, we have connected all the nodes to 0, and computed
            //    the correct minimum distances.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 July 2010
            //
            //  Author:
            //
            //    Original C version by Norm Matloff, CS Dept, UC Davis.
            //    This C++ version by John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int OHD[NV][NV], the distance of the direct link between
            //    nodes I and J.
            //
            //    Output, int DIJKSTRA_DISTANCE[NV], the minimum distance from node 0 
            //    to each node.
            //
        {
            bool[] connected;
            int i;
            int md = 0;
            int[] mind;
            int mv = 0;
            int step;
            //
            //  Start out with only node 0 connected to the tree.
            //
            connected = new bool[ohd.Length];
            connected[0] = true;
            for (i = 1; i < connected.Length; i++)
            {
                connected[i] = false;
            }

            //
            //  Initialize the minimum distance to the one-step distance.
            //
            mind = new int[ohd.Length];

            for (i = 0; i < mind.Length; i++)
            {
                mind[i] = ohd[0][i];
            }

            //
            //  Attach one more node on each iteration.
            //
            for (step = 1; step < ohd.Length; step++)
            {
                //
                //  Find the nearest unconnected node.
                //
                find_nearest(mind, connected, ref md, ref mv);

                if (mv == -1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("DIJKSTRA_DISTANCE - Warning!");
                    Console.WriteLine("  Search terminated early.");
                    Console.WriteLine("  Graph might not be connected.");
                    break;
                }

                //
                //  Mark this node as connected.
                //
                connected[mv] = true;
                //
                //  Having determined the minimum distance to node MV, see if
                //  that reduces the minimum distance to other nodes.
                //
                update_mind(mv, connected, ohd, ref mind);
            }

            return mind;
        }
        
        public static void find_nearest ( int[] mind, bool[] connected, ref int d, ref int v )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FIND_NEAREST finds the nearest unconnected node.
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
        //    Input, int MIND[NV], the currently computed minimum distance from
        //    node 0 to each node.
        //
        //    Input, bool CONNECTED[NV], is true for each connected node, whose 
        //    minimum distance to node 0 has been determined.
        //
        //    Output, int *D, the distance from node 0 to the nearest unconnected node.
        //
        //    Output, int *V, the index of the nearest unconnected node.
        //
        {
            int i;
            int i4_huge = 2147483647;

            d = i4_huge;
            v = -1;
            for ( i = 0; i < mind.Length; i++ )
            {
                if ( !connected[i] && mind[i] < d )
                {
                    d = mind[i];
                    v = i;
                }
            }
        }
        
        public static void update_mind ( int mv, bool[] connected, int[][] ohd, ref int[] mind )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UPDATE_MIND updates the minimum distance vector.
        //
        //  Discussion:
        //
        //    We've just determined the minimum distance to node MV.
        //
        //    For each node I which is not connected yet,
        //    check whether the route from node 0 to MV to I is shorter
        //    than the currently known minimum distance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 July 2010
        //
        //  Author:
        //
        //    Original C version by Norm Matloff, CS Dept, UC Davis.
        //    This C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int MV, the node whose minimum distance to node 0
        //    has just been determined.
        //
        //    Input, bool CONNECTED[NV], is true for each connected node, whose 
        //    minimum distance to node 0 has been determined.
        //
        //    Input, int OHD[NV][NV], the distance of the direct link between
        //    nodes I and J.
        //
        //    Input/output, int MIND[NV], the currently computed minimum distances
        //    from node 0 to each node.
        //
        {
            int i;
            const int i4_huge = 2147483647;

            for ( i = 0; i < connected.Length; i++ )
            {
                if ( !connected[i] )
                {
                    //
                    //  If we really use the maximum integer (or something close) to indicate
                    //  no link, then we'll get burned if we add it to another value;
                    //  Integer arithmetic can "wrap around", so that 17 + i4_huge becomes
                    //  a very negative number!  So first we eliminate the possiblity that
                    //  the link is infinite.
                    //
                    if ( ohd[mv][i] < i4_huge )
                    {
                        if ( mind[mv] + ohd[mv][i] < mind[i] )  
                        {
                            mind[i] = mind[mv] + ohd[mv][i];
                        }
                    }
                }
            }
        }
    }
}