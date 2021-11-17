using System;

namespace Burkardt.MinDist;

public static class BellmanFord
{
    public static void bellman_ford(int v_num, int e_num, int source, int[] e,
            double[] e_weight, ref double[] v_weight, ref int[] predecessor )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BELLMAN_FORD finds shortest paths from a given vertex of a weighted directed graph.
        //
        //  Discussion:
        //
        //    The Bellman-Ford algorithm is used.
        //
        //    Each edge of the graph has a weight, which may be negative.  However,
        //    it should not be the case that there is any negative loop, that is,
        //    a circuit whose total weight is negative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int V_NUM, the number of vertices.
        //
        //    Input, int E_NUM, the number of edges.
        //
        //    Input, int SOURCE, the vertex from which distances will be calculated.
        //
        //    Input, int E[2*E_NUM], the edges, given as pairs of vertex indices.
        //
        //    Input, double E_WEIGHT[E_NUM], the weight of each edge.
        //
        //    Output, double V_WEIGHT[V_NUM], the weight of each node, that is,
        //    its minimum distance from SOURCE.
        //
        //    Output, int PREDECESSOR[V_NUM], a list of predecessors, which can be
        //    used to recover the shortest path from any node back to SOURCE.
        //    
    {
        int i;
        int j;
        const double r8_big = 1.0E+30;
        double t;
        int u;
        int v;
        //
        //  Step 1: initialize the graph.
        //
        for (i = 0; i < v_num; i++)
        {
            if (i == source)
            {
                v_weight[i] = 0.0;
            }
            else
            {
                v_weight[i] = r8_big;
            }
        }

        for (i = 0; i < v_num; i++)
        {
            predecessor[i] = -1;
        }

        //
        //  Step 2: Relax edges repeatedly.
        //
        for (i = 1; i < v_num; i++)
        {
            for (j = 0; j < e_num; j++)
            {
                u = e[1 + j * 2];
                v = e[0 + j * 2];
                t = v_weight[u] + e_weight[j];
                if (t < v_weight[v])
                {
                    v_weight[v] = t;
                    predecessor[v] = u;
                }
            }
        }

        //
        //  Step 3: check for negative-weight cycles
        //
        for (j = 0; j < e_num; j++)
        {
            u = e[1 + j * 2];
            v = e[0 + j * 2];
            if (v_weight[u] + e_weight[j] < v_weight[v])
            {
                Console.WriteLine("");
                Console.WriteLine("BELLMAN_FORD - Fatal error!");
                Console.WriteLine("  Graph contains a cycle with negative weight.");
                return;
            }
        }
    }
}