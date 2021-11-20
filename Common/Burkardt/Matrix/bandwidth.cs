using System;

namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static int bandwidth ( int nnodes, int element_num, int[] element_node,
            int node_num, int[] indx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BANDWIDTH determines the bandwidth of the coefficient matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NNODES, the number of local nodes per element.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
        //    index of local node I in element J.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int INDX[NODE_NUM], indicates how the value associated with the
        //    node is to be determined.  If INDX(I) is positive, then this is the
        //    index of the unknown in the finite element linear system.  The value
        //    at the node will be determined by solving the finite element system.
        //    If INDX(I) is negative, then the node is associated with a boundary
        //    condition; the value of the boundary condition is stored in the array
        //    UB, in location -INDX(I).
        //
        //    Output, int BANDWIDTH, the half bandwidth of the matrix.
        //
    {
        int element;

        int nhba = 0;

        for ( element = 1; element <= element_num; element++ )
        {
            int iln;
            for ( iln = 1; iln <= nnodes; iln++ )
            {
                int in_ = element_node[iln-1+(element-1)*nnodes];
                int i = indx[in_-1];
                switch (i)
                {
                    case > 0:
                    {
                        int jln;
                        for ( jln = 1; jln <= nnodes; jln++ )
                        {
                            int jn = element_node[jln-1+(element-1)*nnodes];
                            int j = indx[jn-1];
                            nhba = Math.Max ( nhba, j - i );
                        }

                        break;
                    }
                }
            }
        }

        return nhba;
    }
        
    public static int bandwidth ( int element_num, int[] element_node )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BANDWIDTH determines the bandwidth of the coefficient matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input,  ELEMENT_NODE[3*ELEMENT_NUM];
        //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
        //
        //    Output, int BANDWIDTH, the half bandwidth of the matrix.
        //
    {
        int element;

        int nhba = 0;

        for ( element = 0; element < element_num; element++ )
        {
            int local_i;
            for ( local_i = 0; local_i < 3; local_i++ )
            {
                int global_i = element_node[local_i+element*3];
                int local_j;
                for ( local_j = 1; local_j <= 3; local_j++ )
                {
                    int global_j = element_node[local_j+element*3];
                    nhba = Math.Max ( nhba, Math.Abs ( global_j - global_i ) );
                }
            }
        }

        return nhba;
    }
    public static void bandwidth ( int m, int n, double[] a, ref int b, ref int l, ref int d, ref int u )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BANDWIDTH returns the bandwidth of a matrix.
        //
        //  Discussion:
        //
        //    If the nonzeros of a matrix only occur in entries that are "close"
        //    to the main diagonal, we say the matrix is banded.
        //
        //    Roughly speaking, the bandwidth B of a matrix is the number of 
        //    diagonals containing nonzeros.  More precisely, it is the minimum number
        //    of contiguous diagonals that contain all the nonzeros.  It is presumed
        //    that the main diagonal is nonzero.
        //
        //    We can also measure U and L, the upper and lower "half-bandwidths" which
        //    count the number of contiguous diagonals above or below the main
        //    diagonal.
        //
        //    We may write
        //      B = L + D + U
        //    where D is presumably 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the matrix.
        //
        //    Output, int &B, the total bandwidth.
        //
        //    Output, int &L, &D, &U, the lower, diagonal, and upper 
        //    bandwidths.
        //
    {
        int i;

        l = 0;
        d = 0;
        u = 0;

        for ( i = 0; i < n; i++ )
        {
            int j = 0;
            while ( l < i - j )
            {
                if ( a[i+j*m] != 0.0 )
                {
                    l = i - j;
                    break;
                }
                j += 1;
            }

            if ( a[i+i*m] != 0.0 )
            {
                d = 1;
            }

            j = n - 1;
            while ( u < j - i )
            {
                if ( a[i+j*m] != 0.0 )
                {
                    u = j - i;
                    break;
                }
                j -= 1;
            }
        }

        b = l + d + u;
    }
}