namespace Burkardt.Graph;

public static class Adjacency
{
    public static void graph_01_adj(int node_num, int adj_num, ref int[] adj_row, ref int[] adj)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_01_ADJ returns the adjacency vector for graph 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan George, Joseph Liu,
        //    Computer Solution of Large Sparse Positive Definite Systems,
        //    Prentice Hall, 1981.
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ADJ_NUM, the number of adjacencies.
        //
        //    Output, int ADJ_ROW[NODE_NUM+1], node pointers into ADJ.
        //
        //    Output, int ADJ[ADJ_NUM], the adjacency information.
        //
    {
        const int ADJ_NUM = 28;
        const int NODE_NUM = 10;

        int[] adj_save =
        {
            4, 6,
            3, 5, 7, 10,
            2, 4, 5,
            1, 3, 6, 9,
            2, 3, 7,
            1, 4, 7, 8,
            2, 5, 6, 8,
            6, 7,
            4,
            2
        };
        int[] adj_row_save =
        {
            1, 3, 7, 10, 14, 17, 21, 25, 27, 28, 29
        };
        int i;

        for (i = 0; i < ADJ_NUM; i++)
        {
            adj[i] = adj_save[i];
        }

        for (i = 0; i < NODE_NUM + 1; i++)
        {
            adj_row[i] = adj_row_save[i];
        }
    }

    public static void graph_01_size(ref int node_num, ref int adj_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_01_SIZE returns the number of adjacencies for graph 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan George, Joseph Liu,
        //    Computer Solution of Large Sparse Positive Definite Systems,
        //    Prentice Hall, 1981.
        //
        //  Parameters:
        //
        //    Output, int *NODE_NUM, the number of items that can be adjacent.
        //
        //    Output, int *ADJ_NUM, the number of adjacencies.
        //
    {
        node_num = 10;
        adj_num = 28;
    }
}