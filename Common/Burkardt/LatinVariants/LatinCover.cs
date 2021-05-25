using Burkardt.Types;

namespace Burkardt.Latin
{
    public static partial class LatinVariants
    {
        public static int[] latin_cover(int n, int[] p)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LATIN_COVER returns a 2D Latin Square Covering.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, int P[N], a permutation which describes the
            //    first Latin square.
            //
            //    Output, int LATIN_COVER[N*N], the Latin cover.  A(I,J) = K
            //    means that (I,J) is one element of the K-th Latin square.
            //
        {
            int[] a = new int[n * n];

            typeMethods.perm_check(n, p);

            for (int i = 0; i < n; i++)
            {
                for (int k = 0; k < n; k++)
                {
                    int ik = typeMethods.i4_wrap(i + k, 0, n - 1);
                    int j = p[ik] - 1;
                    a[i + j * n] = k + 1;
                }
            }

            return a;
        }

        public static int[] latin_cover_2d(int n, int[] p1, int[] p2)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LATIN_COVER_2D returns a 2D Latin Square Covering.
            //
            //  Discussion:
            //
            //    This procedure has a chance of being extended to M dimensions.
            //
            //    A basic solution is computed, and the user is permitted to permute
            //    both the I and J coordinates.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, int P1[N], P2[N], permutations to be applied
            //    to the spatial dimensions.
            //
            //    Output, int LATIN_COVER_2D[N*N], the Latin cover.  A(I,J) = K
            //    means that (I,J) is one element of the K-th Latin square.
            //
        {
            typeMethods.perm_check(n, p1);
            typeMethods.perm_check(n, p2);

            int[] a = new int[n * n];
            int[] b = new int[n * n];
            //
            //  Set up the basic solution.
            //
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    a[i + j * n] = typeMethods.i4_wrap(i - j + 1, 1, n);
                }
            }

            //
            //  Apply permutation to dimension I.
            //
            for (int i = 0; i < n; i++)
            {
                int i1 = p1[i] - 1;
                for (int j = 0; j < n; j++)
                {
                    b[i1 + j * n] = a[i + j * n];
                }
            }

            //
            //  Apply permutation to dimension J.
            //
            for (int j = 0; j < n; j++)
            {
                int j1 = p2[j] - 1;
                for (int i = 0; i < n; i++)
                {
                    a[i + j1 * n] = b[i + j * n];
                }
            }

            return a;
        }

        public static int[] latin_cover_3d(int n, int[] p1, int[] p2, int[] p3)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LATIN_COVER_3D returns a 3D Latin Square Covering.
        //
        //  Discussion:
        //
        //    A basic solution is computed, and the user is permitted to permute
        //    I, J and K coordinates.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, int P1[N], P2[N], P3[N], permutations to be applied
        //    to the spatial dimensions.
        //
        //    Output, int LATIN_COVER_3D[N*N*N], the Latin cover.  A(I,J,K) = L
        //    means that (I,J,K) is one element of the L-th Latin square.
        //
        {

            typeMethods.perm_check(n, p1);
            typeMethods.perm_check(n, p2);
            typeMethods.perm_check(n, p3);

            int[] a = new int[n * n * n];
            int[] b = new int[n * n * n];
            //
            //  Set up the basic solution.
            //
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        int ik = typeMethods.i4_wrap(i + 1 - k, 1, n);
                        int jk = typeMethods.i4_wrap(j + 1 - k, 1, n);
                        b[i + j * n + k * n * n] = ik + (jk - 1) * n;
                    }
                }
            }

            //
            //  Apply permutation to dimension I.
            //
            for (int i = 0; i < n; i++)
            {
                int i1 = p1[i] - 1;
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        a[i1 + j * n + k * n * n] = b[i + j * n + k * n * n];
                    }
                }
            }

            //
            //  Apply permutation to dimension J.
            //
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    int j1 = p2[j] - 1;
                    for (int k = 0; k < n; k++)
                    {
                        b[i + j1 * n + k * n * n] = a[i + j * n + k * n * n];
                    }
                }
            }

            //
            //  Apply permutation to dimension K.
            //
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        int k1 = p3[k] - 1;
                        a[i + j * n + k1 * n * n] = b[i + j * n + k * n * n];
                    }
                }
            }

            return a;
        }
    }
}