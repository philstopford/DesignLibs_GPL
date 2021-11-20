using System;
using System.Globalization;

namespace Burkardt.MatrixNS;

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
        int i;

        int band_lo = 0;
        int band_hi = 0;

        for (i = 0; i < node_num; i++)
        {
            int j;
            for (j = adj_row[i]; j <= adj_row[i + 1] - 1; j++)
            {
                int col = adj[j % adj.Length];
                band_lo = Math.Max(band_lo, i - col);
                band_hi = Math.Max(band_hi, col - i);
            }
        }

        int value = band_lo + 1 + band_hi;

        return value;
    }
        
    public static bool adj_contains_ij ( int node_num, int adj_num, int[] adj_row, int[] adj,
            int i, int j )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ADJ_CONTAINS_IJ determines if (I,J) is in an adjacency structure.
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
        //
        //    Input, int I, J, the two nodes, for which we want to know
        //    whether I is adjacent to J.
        //
        //    Output, bool ADJ_CONTAINS_IJ, is TRUE if I = J, or the adjacency
        //    structure contains the information that I is adjacent to J.
        //
    {
        int k;
        //
        //  Symmetric entries are not stored.
        //
        if ( i == j )
        {
            return true;
        }
        //
        //  Illegal I, J entries.
        //
        if ( node_num < i )
        {
            return false;
        }

        switch (i)
        {
            case < 1:
                return false;
        }
        if ( node_num < j )
        {
            return false;
        }
        switch (j)
        {
            case < 1:
                return false;
        }

        //
        //  Search the adjacency entries already stored for row I,
        //  to see if J has already been stored.
        //
        int klo = adj_row[i-1];
        int khi = adj_row[i]-1;

        for ( k = klo; k <= khi; k++ )
        {
            if (adj[k - 1] != j)
            {
                continue;
            }

            return true;
        }

        return false;
    }
        
    public static void adj_insert_ij ( int node_num, int adj_max, ref int adj_num, ref int[] adj_row,
            ref int[] adj, int i, int j )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ADJ_INSERT_IJ inserts (I,J) into an adjacency structure.
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
        //    Input, int ADJ_MAX, the maximum number of adjacency entries.
        //
        //    Input/output, int ADJ_NUM, the number of adjacency entries.
        //
        //    Input/output, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
        //    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
        //
        //    Input/output, int ADJ[ADJ_NUM], the adjacency structure.
        //
        //    Input, int I, J, the two nodes which are adjacent.
        //
    {
        int k;
        //
        //  A new adjacency entry must be made.
        //  Check that we're not exceeding the storage allocation for ADJ.
        //
        if ( adj_max < adj_num + 1 )
        {
            Console.WriteLine("");
            Console.WriteLine("ADJ_INSERT_IJ - Fatal error!");
            Console.WriteLine("  All available storage has been used.");
            Console.WriteLine("  No more information can be stored!");
            Console.WriteLine("  This error occurred for ");
            Console.WriteLine("  Row I =    " + i + "");
            Console.WriteLine("  Column J = " + j + "");
            return;
        }
        //
        //  The action is going to occur between ADJ_ROW(I) and ADJ_ROW(I+1)-1:
        //
        int j_spot = adj_row[i-1];

        for ( k = adj_row[i-1]; k <= adj_row[i]-1; k++ )
        {
            if ( adj[k-1] == j )
            {
                return;
            }

            if ( adj[k-1] < j )
            {
                j_spot = k + 1;
            }
            else
            {
                break;
            }
        }

        for ( k = adj_num; j_spot <= k; k-- )
        {
            adj[k] = adj[k-1];
        }
        adj[j_spot-1] = j;

        for ( k = i; k <= node_num; k++ )
        {
            adj_row[k] += 1;
        }

        adj_num += 1;

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
        int i;

        int band_lo = 0;
        int band_hi = 0;

        for (i = 0; i < node_num; i++)
        {
            int j;
            for (j = adj_row[(perm[i] + adj_row.Length) % adj_row.Length]; j <= adj_row[(perm[i] + 1 + adj_row.Length) % adj_row.Length] - 1; j++)
            {
                int col = perm_inv[(adj[(j + adj.Length) % adj.Length] + perm_inv.Length) % perm_inv.Length];
                band_lo = Math.Max(band_lo, i - col);
                band_hi = Math.Max(band_hi, col - i);
            }
        }

        int bandwidth = band_lo + 1 + band_hi;

        return bandwidth;
    }

    public static void adj_perm_show(int node_num, int adj_num, int[] adj_row, int[] adj,
            int[] perm, int[] perm_inv)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ADJ_PERM_SHOW displays a symbolic picture of a permuted adjacency matrix.
        //
        //  Discussion:
        //
        //    The matrix is defined by the adjacency information and a permutation.
        //
        //    The routine also computes the bandwidth and the size of the envelope.
        //
        //    If no permutation has been done, you must set PERM(I) = PERM_INV(I) = I
        //    before calling this routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 January 2007
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
        //    Input, int PERM[NODE_NUM], PERM_INV[NODE_NUM], the permutation
        //    and inverse permutation.
        //
    {
        int i;

        char[] band = new char[node_num];

        int band_lo = 0;
        int nonzero_num = 0;

        Console.WriteLine("");
        Console.WriteLine("  Nonzero structure of matrix:");
        Console.WriteLine("");

        for (i = 0; i < node_num; i++)
        {
            int k;
            for (k = 0; k < node_num; k++)
            {
                band[k] = '.';
            }

            band[i] = 'D';

            int j;
            for (j = adj_row[perm[i] - 1]; j <= adj_row[perm[i]] - 1; j++)
            {
                int col = perm_inv[adj[j - 1] - 1] - 1;

                if (col < i)
                {
                    nonzero_num += 1;
                }

                band_lo = Math.Max(band_lo, i - col);

                if (col != i)
                {
                    band[col] = 'X';
                }
            }

            string cout = "  " + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8) + " ";
            for (j = 0; j < node_num; j++)
            {
                cout += band[j];
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Lower bandwidth = " + band_lo + "");
        Console.WriteLine("  Lower envelope contains " + nonzero_num + " nonzeros.");

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

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("  Sparse adjacency structure:");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes       = " + node_num + "");
        Console.WriteLine("  Number of adjacencies = " + adj_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Node   Min   Max          Nonzeros ");
        Console.WriteLine("");

        for (i = node_lo; i <= node_hi; i++)
        {
            int jmin = adj_row[i];
            int jmax = adj_row[i + 1] - 1;

            if (jmax < jmin)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + jmin.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + jmax.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
            }
            else
            {
                int jlo;
                for (jlo = jmin; jlo <= jmax; jlo += 5)
                {
                    int jhi = Math.Min(jlo + 4, jmax);

                    int j;
                    if (jlo == jmin)
                    {
                        string cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                           + "  " + jmin.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                           + "  " + jmax.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                           + "   ";
                        for (j = jlo; j <= jhi; j++)
                        {
                            cout += adj[j % adj.Length].ToString(CultureInfo.InvariantCulture).PadLeft(8);
                        }

                        Console.WriteLine(cout);
                    }
                    else
                    {
                        string cout2 = "                     ";
                        for (j = jlo; j <= jhi; j++)
                        {
                            cout2 += adj[j % adj.Length].ToString(CultureInfo.InvariantCulture).PadLeft(8);
                        }

                        Console.WriteLine(cout2);
                    }
                }
            }
        }
    }

    public static void adj_set(int node_num, int adj_max, ref int adj_num, ref int[] adj_row,
            ref int[] adj, int irow, int jcol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ADJ_SET sets up the adjacency information.
        //
        //  Discussion:
        //
        //    The routine records the locations of each nonzero element,
        //    one at a time.
        //
        //    The first call for a given problem should be with IROW or ICOL
        //    negative.  This is a signal indicating the data structure should
        //    be initialized.
        //
        //    Then, for each case in which A(IROW,JCOL) is nonzero, or
        //    in which IROW is adjacent to JCOL, call this routine once
        //    to record that fact.
        //
        //    Diagonal entries are not to be stored.
        //
        //    The matrix is assumed to be symmetric, so setting I adjacent to J
        //    will also set J adjacent to I.
        //
        //    Repeated calls with the same values of IROW and JCOL do not
        //    actually hurt.  No extra storage will be allocated.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ADJ_MAX, the maximum dimension of the adjacency array.
        //
        //    Input/output, int *ADJ_NUM, the number of adjaceny entries.
        //
        //    Input/output, int ADJ_ROW[NODE_NUM+1].  Information about
        //    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
        //
        //    Input/output, int ADJ[ADJ_NUM], the adjacency structure.
        //
        //    Input, int IROW, JCOL, the row and column indices of a nonzero
        //    entry of the matrix.
        //
    {
        //
        //  Negative IROW or JCOL indicates the data structure should be initialized.
        //
        if (irow < 0 || jcol < 0)
        {
            Console.WriteLine("");
            Console.WriteLine("ADJ_SET - Note:");
            Console.WriteLine("  Initializing adjacency information.");
            Console.WriteLine("  Number of nodes NODE_NUM =  " + node_num + "");
            Console.WriteLine("  Maximum adjacency ADJ_MAX = " + adj_max + "");

            adj_num = 0;
            int i;
            for (i = 0; i < node_num + 1; i++)
            {
                adj_row[i] = 1;
            }

            for (i = 0; i < adj_max; i++)
            {
                adj[i] = 0;
            }

            return;
        }

        //
        //  Diagonal entries are not stored.
        //
        if (irow == jcol)
        {
            return;
        }

        if (node_num < irow)
        {
            Console.WriteLine("");
            Console.WriteLine("ADJ_SET - Fatal error!");
            Console.WriteLine("  NODE_NUM < IROW.");
            Console.WriteLine("  IROW =     " + irow + "");
            Console.WriteLine("  NODE_NUM = " + node_num + "");
            return;
        }

        switch (irow)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("ADJ_SET - Fatal error!");
                Console.WriteLine("  IROW < 1.");
                Console.WriteLine("  IROW = " + irow + "");
                return;
        }
        if (node_num < jcol)
        {
            Console.WriteLine("");
            Console.WriteLine("ADJ_SET - Fatal error!");
            Console.WriteLine("  NODE_NUM < JCOL.");
            Console.WriteLine("  JCOL =     " + jcol + "");
            Console.WriteLine("  NODE_NUM = " + node_num + "");
            return;
        }
        switch (jcol)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("ADJ_SET - Fatal error!");
                Console.WriteLine("  JCOL < 1.");
                Console.WriteLine("  JCOL = " + jcol + "");
                return;
        }

        if (!adj_contains_ij(node_num, adj_num, adj_row, adj, irow, jcol))
        {
            adj_insert_ij(node_num, adj_max, ref adj_num, ref adj_row, ref adj, irow, jcol);
        }

        if (!adj_contains_ij(node_num, adj_num, adj_row, adj, jcol, irow))
        {
            adj_insert_ij(node_num, adj_max, ref adj_num, ref adj_row, ref adj, jcol, irow);
        }

    }

    public static void adj_show(int node_num, int adj_num, int[] adj_row, int[] adj)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ADJ_SHOW displays a symbolic picture of an adjacency matrix.
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
        //    Input, int ADJ_NUM, the number of adjacency entries.
        //
        //    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
        //    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
        //
        //    Input, int ADJ[ADJ_NUM], the adjacency structure.
        //    For each row, it contains the column indices of the nonzero entries.
        //
    {
        int i;

        char[] band = new char[node_num];

        int band_lo = 0;
        int nonzero_num = 0;

        Console.WriteLine("");
        Console.WriteLine("  Nonzero structure of matrix:");
        Console.WriteLine("");

        for (i = 0; i < node_num; i++)
        {
            int k;
            for (k = 0; k < node_num; k++)
            {
                band[k] = '.';
            }

            band[i] = 'D';

            int j;
            for (j = adj_row[i]; j <= adj_row[i + 1] - 1; j++)
            {
                int col = adj[j - 1] - 1;
                if (col < i)
                {
                    nonzero_num += 1;
                }

                band_lo = Math.Max(band_lo, i - col);
                band[col] = 'X';
            }

            string cout = "  " + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8) + " ";
            for (j = 0; j < node_num; j++)
            {
                cout += band[j];
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Lower bandwidth = " + band_lo + "");
        Console.WriteLine("  Lower envelope contains " + nonzero_num + " nonzeros.");
    }

}