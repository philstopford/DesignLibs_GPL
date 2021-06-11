using System;

namespace Burkardt
{
    public static class WathenMatrix
    {
        public static double[] wathen(int nx, int ny, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WATHEN returns the WATHEN matrix.
            //
            //  Discussion:
            //
            //    The Wathen matrix is a finite element matrix which is sparse.
            //
            //    The entries of the matrix depend in part on a physical quantity
            //    related to density.  That density is here assigned random values between
            //    0 and 100.
            //
            //    The matrix order N is determined by the input quantities NX and NY,
            //    which would usually be the number of elements in the X and Y directions.
            //    The value of N is
            //
            //      N = 3*NX*NY + 2*NX + 2*NY + 1,
            //
            //    and sufficient storage in A must have been set aside to hold
            //    the matrix.
            //
            //    A is the consistent mass matrix for a regular NX by NY grid
            //    of 8 node serendipity elements.  
            //
            //    Here is an illustration for NX = 3, NY = 2:
            //
            //     23-24-25-26-27-28-29
            //      |     |     |     |
            //     19    20    21    22
            //      |     |     |     |
            //     12-13-14-15-16-17-18
            //      |     |     |     |
            //      8     9    10    11
            //      |     |     |     |
            //      1--2--3--4--5--6--7
            //
            //    For this example, the total number of nodes is, as expected,
            //
            //      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
            //
            //  Properties:
            //
            //    A is symmetric positive definite for any positive values of the
            //    density RHO(NX,NY), which is here given the value 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Nicholas Higham,
            //    Algorithm 694: A Collection of Test Matrices in MATLAB,
            //    ACM Transactions on Mathematical Software,
            //    Volume 17, Number 3, September 1991, pages 289-305.
            //
            //    Andrew Wathen,
            //    Realistic eigenvalue bounds for the Galerkin mass matrix,
            //    IMA Journal of Numerical Analysis,
            //    Volume 7, 1987, pages 449-457.
            //
            //  Parameters:
            //
            //    Input, int NX, NY, values which determine the size of A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, double WATHEN[N*N], the matrix.
            //
        {
            double[] a;
            double[] em = {
                6.0, -6.0, 2.0, -8.0, 3.0, -8.0, 2.0, -6.0,
                -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0,
                2.0, -6.0, 6.0, -6.0, 2.0, -8.0, 3.0, -8.0,
                -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0,
                3.0, -8.0, 2.0, -6.0, 6.0, -6.0, 2.0, -8.0,
                -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0,
                2.0, -8.0, 3.0, -8.0, 2.0, -6.0, 6.0, -6.0,
                -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0
            }
            ;
            int i;
            int j;
            int kcol;
            int krow;
            int[] node = new int[8];
            double rho;

            a = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    a[i + j * n] = 0.0;
                }
            }

            for (j = 1; j <= ny; j++)
            {
                for (i = 1; i <= nx; i++)
                {
                    //
                    //  For the element (I,J), determine the indices of the 8 nodes.
                    //
                    node[0] = 3 * j * nx + 2 * j + 2 * i;
                    node[1] = node[0] - 1;
                    node[2] = node[0] - 2;
                    node[3] = (3 * j - 1) * nx + 2 * j + i - 2;
                    node[4] = (3 * j - 3) * nx + 2 * j + 2 * i - 4;
                    node[5] = node[4] + 1;
                    node[6] = node[4] + 2;
                    node[7] = node[3] + 1;
                    //
                    //  The density RHO can also be set to a random positive value.
                    //
                    for (krow = 0; krow < 8; krow++)
                    {
                        for (kcol = 0; kcol < 8; kcol++)
                        {
                            rho = 1.0;

                            if (node[krow] < 0 || n <= node[krow] ||
                                node[kcol] < 0 || n <= node[kcol])
                            {
                                Console.WriteLine("");
                                Console.WriteLine("WATHEN - Fatal error!");
                                Console.WriteLine("  I = " + i + "  J = " + j + "");
                                Console.WriteLine("  KROW = " + krow + "");
                                Console.WriteLine("  KCOL = " + kcol + "");
                                Console.WriteLine("  NODE[KROW] = " + node[krow] + "");
                                Console.WriteLine("  NODE[KCOL] = " + node[kcol] + "");
                                return null;
                            }

                            a[node[krow] + node[kcol] * n] = a[node[krow] + node[kcol] * n]
                                                             + 20.0 * rho * em[krow + kcol * 8] / 9.0;
                        }
                    }
                }
            }

            return a;
        }

        public static int wathen_order(int nx, int ny)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WATHEN_ORDER returns the order of the WATHEN matrix.
            //
            //  Discussion:
            //
            //    N = 3*NX*NY + 2*NX + 2*NY + 1,
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 June 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Nicholas Higham,
            //    Algorithm 694: A Collection of Test Matrices in MATLAB,
            //    ACM Transactions on Mathematical Software,
            //    Volume 17, Number 3, September 1991, pages 289-305.
            //
            //    Andrew Wathen,
            //    Realistic eigenvalue bounds for the Galerkin mass matrix,
            //    IMA Journal of Numerical Analysis,
            //    Volume 7, 1987, pages 449-457.
            //
            //  Parameters:
            //
            //    Input, int NX, NY, values which determine the size of A.
            //
            //    Output, int WATHEN_ORDER, the order of the matrix.
            //
        {
            int n;

            n = 3 * nx * ny + 2 * nx + 2 * ny + 1;

            return n;
        }
    }
}