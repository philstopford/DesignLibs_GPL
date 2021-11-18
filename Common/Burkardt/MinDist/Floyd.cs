using System;

namespace Burkardt.MinDist;

public static class Floyd
{
    public static void i4mat_floyd ( int n, ref int[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_FLOYD: shortest distances between pairs of nodes in a directed graph.
        //
        //  Discussion:
        //
        //    We assume we are given the adjacency matrix A of the directed graph.
        //
        //    We assume that A is an I4MAT, that is, a two-dimensional array of I4's.
        //
        //    The adjacency matrix is NOT assumed to be symmetric.
        //
        //    If there is not a direct link from node I to node J, the distance
        //    would formally be infinity.  We assume that such distances are assigned
        //    some large value.  For technical reasons, this number should be no
        //    greater than 2147483647/2 so that the "infinity+infinity" does not
        //    suddenly become a negative value!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, int A[N*N].
        //    On input, A(I,J) contains the direct distance from node I to node J.
        //    On output, A(I,J) contains the shortest distance from node I to node J.
        //    An input entry of -1 means there is no direct route from node I to node J.
        //    An output entry of -1 means there
        //
    {
        const int i4_huge = 2147483647;

        for (int k = 0; k < n; k++ )
        {
            for (int j = 0; j < n; j++ )
            {
                switch (a[k+j*n])
                {
                    case < i4_huge:
                    {
                        for (int i = 0; i < n; i++ )
                        {
                            a[i + j * n] = a[i + k * n] switch
                            {
                                < i4_huge => Math.Min(a[i + j * n], a[i + k * n] + a[k + j * n]),
                                _ => a[i + j * n]
                            };
                        }

                        break;
                    }
                }
            }
        }
    }
        
    public static void r8mat_floyd ( int n, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_FLOYD: shortest distance between pairs of nodes in a directed graph.
        //
        //  Discussion:
        //
        //    We assume we are given the adjacency matrix A of the directed graph.
        //
        //    We assume that A is an R8MAT, that is, a two-dimensional array of R8's.
        //
        //    The adjacency matrix is NOT assumed to be symmetric.
        //
        //    If there is not a direct link from node I to node J, the distance
        //    would formally be infinity.  We assume that such distances are assigned
        //    the value R8_HUGE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, double A[N*N].
        //    On input, A(I,J) contains the direct distance from node I to node J.
        //    On output, A(I,J) contains the shortest distance from node I to node J.
        //
    {
        const double r8_huge = 1.0E+30;

        for (int k = 0; k < n; k++ )
        {
            for (int j = 0; j < n; j++ )
            {
                switch (a[k+j*n])
                {
                    case < r8_huge:
                    {
                        for (int i = 0; i < n; i++ )
                        {
                            a[i + j * n] = a[i + k * n] switch
                            {
                                < r8_huge => Math.Min(a[i + j * n], a[i + k * n] + a[k + j * n]),
                                _ => a[i + j * n]
                            };
                        }

                        break;
                    }
                }
            }
        }
    }
}