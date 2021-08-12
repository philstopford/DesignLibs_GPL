using System;

namespace Burkardt.MatrixNS
{
    public static class AdjacencyMatrix
    {
        public static int adj_bandwidth(int node_num, int adj_num, int[] adj_row, int[] adj)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 January 2007
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
            //    Input, int ADJ_NUM, the number of adjacency entries.
            //
            //    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
            //    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
            //
            //    Input, int ADJ[ADJ_NUM], the adjacency structure.
            //    For each row, it contains the column indices of the nonzero entries.
            //
            //    Output, int ADJ_BANDWIDTH, the bandwidth of the adjacency
            //    matrix.
            //
        {
            int band_hi;
            int band_lo;
            int col;
            int i;
            int j;
            int value;

            band_lo = 0;
            band_hi = 0;

            for (i = 0; i < node_num; i++)
            {
                for (j = adj_row[i]; j <= adj_row[i + 1] - 1; j++)
                {
                    col = adj[j];
                    band_lo = Math.Max(band_lo, i - col);
                    band_hi = Math.Max(band_hi, col - i);
                }
            }

            value = band_lo + 1 + band_hi;

            return value;
        }

        public static int adj_perm_bandwidth(int node_num, int adj_num, int[] adj_row, int[] adj,
                int[] perm, int[] perm_inv)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
            //
            //  Discussion:
            //
            //    The matrix is defined by the adjacency information and a permutation.  
            //
            //    The routine also computes the bandwidth and the size of the envelope.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 January 2007
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
            //    Input, int ADJ_NUM, the number of adjacency entries.
            //
            //    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
            //    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
            //
            //    Input, int ADJ[ADJ_NUM], the adjacency structure.
            //    For each row, it contains the column indices of the nonzero entries.
            //
            //    Input, int PERM[NODE_NUM], PERM_INV(NODE_NUM), the permutation
            //    and inverse permutation.
            //
            //    Output, int ADJ_PERM_BANDWIDTH, the bandwidth of the permuted 
            //    adjacency matrix.
            //
        {
            int band_hi;
            int band_lo;
            int bandwidth;
            int col;
            int i;
            int j;

            band_lo = 0;
            band_hi = 0;

            for (i = 0; i < node_num; i++)
            {
                for (j = adj_row[perm[i]]; j <= adj_row[perm[i] + 1] - 1; j++)
                {
                    col = perm_inv[adj[j]];
                    band_lo = Math.Max(band_lo, i - col);
                    band_hi = Math.Max(band_hi, col - i);
                }
            }

            bandwidth = band_lo + 1 + band_hi;

            return bandwidth;
        }

        public static void adj_print(int node_num, int adj_num, int[] adj_row, int[] adj,
                string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ADJ_PRINT prints adjacency information.
            //
            //  Discussion:
            //
            //    The list has the form:
            //
            //    Row   Nonzeros
            //
            //    1       2   5   9
            //    2       7   8   9   15   78   79   81  86  91  99
            //          100 103
            //    3      48  49  53
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 December 2002
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int ADJ_NUM, the number of adjacency entries.
            //
            //    Input, int ADJ_ROW[NODE_NUM+1], organizes the adjacency entries
            //    into rows.  The entries for row I are in entries ADJ_ROW(I)
            //    through ADJ_ROW(I+1)-1.
            //
            //    Input, int ADJ[ADJ_NUM], the adjacency structure, which contains,
            //    for each row, the column indices of the nonzero entries.
            //
            //    Input, string TITLE, a title.
            //
        {
            adj_print_some(node_num, 0, node_num - 1, adj_num, adj_row, adj, title);

        }

        public static void adj_print_some(int node_num, int node_lo, int node_hi, int adj_num,
                int[] adj_row, int[] adj, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ADJ_PRINT_SOME prints some adjacency information.
            //
            //  Discussion:
            //
            //    The list has the form:
            //
            //    Row   Nonzeros
            //
            //    1       2   5   9
            //    2       7   8   9   15   78   79   81  86  91  99
            //          100 103
            //    3      48  49  53
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
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int NODE_LO, NODE_HI, the first and last nodes for
            //    which the adjacency information is to be printed.
            //
            //    Input, int ADJ_NUM, the number of adjacency entries.
            //
            //    Input, int ADJ_ROW[NODE_NUM+1], organizes the adjacency entries
            //    into rows.  The entries for row I are in entries ADJ_ROW(I)
            //    through ADJ_ROW(I+1)-1.
            //
            //    Input, int ADJ[ADJ_NUM], the adjacency structure, which contains,
            //    for each row, the column indices of the nonzero entries.
            //
            //    Input, string TITLE, a title to be printed.
            //
        {
            int i;
            int j;
            int jhi;
            int jlo;
            int jmax;
            int jmin;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("  Sparse adjacency structure:");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes       = " + node_num + "");
            ;
            Console.WriteLine("  Number of adjacencies = " + adj_num + "");
            Console.WriteLine("");
            Console.WriteLine("  Node   Min   Max          Nonzeros ");
            Console.WriteLine("");

            for (i = node_lo; i <= node_hi; i++)
            {
                jmin = adj_row[i];
                jmax = adj_row[i + 1] - 1;

                if (jmax < jmin)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(4)
                                           + "  " + jmin.ToString().PadLeft(4)
                                           + "  " + jmax.ToString().PadLeft(4) + "");
                }
                else
                {
                    for (jlo = jmin; jlo <= jmax; jlo = jlo + 5)
                    {
                        jhi = Math.Min(jlo + 4, jmax);

                        if (jlo == jmin)
                        {
                            string cout = "  " + i.ToString().PadLeft(4)
                                               + "  " + jmin.ToString().PadLeft(4)
                                               + "  " + jmax.ToString().PadLeft(4)
                                               + "   ";
                            for (j = jlo; j <= jhi; j++)
                            {
                                cout += adj[j].ToString().PadLeft(8);
                            }

                            Console.WriteLine(cout);
                        }
                        else
                        {
                            string cout2 = "                     ";
                            for (j = jlo; j <= jhi; j++)
                            {
                                cout2 += adj[j].ToString().PadLeft(8);
                            }

                            Console.WriteLine(cout2);
                        }
                    }
                }
            }
        }
    }
}