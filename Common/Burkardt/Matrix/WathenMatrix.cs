using System;
using Burkardt.Uniform;

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

        public static void wathen_bandwidth(int nx, int ny, ref int l, ref int d, ref int u)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WATHEN_BANDWIDTH returns the bandwidth of the WATHEN matrix.
            //
            //  Discussion:
            //
            //    The bandwidth measures the minimal number of contiguous diagonals,
            //    including the central diagonal, which contain all the nonzero elements
            //    of a matrix.
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
            //    Output, int &L, &D, &U, the lower, diagonal, and upper 
            //    bandwidths of the matrix,
            //
        {
            l = 3 * nx + 4;
            d = 1;
            u = 3 * nx + 4;

            return;
        }

        public static double[] wathen_gb(int nx, int ny, int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WATHEN_GB returns the Wathen matrix, using general banded (GB) storage.
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
            //    The matrix is the consistent mass matrix for a regular NX by NY grid
            //    of 8 node serendipity elements.
            //
            //    The local element numbering is
            //
            //      3--2--1
            //      |     |
            //      4     8
            //      |     |
            //      5--6--7
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
            //    The matrix is symmetric positive definite for any positive values of the
            //    density RHO(X,Y).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 July 2014
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
            //    Volume 7, Number 4, October 1987, pages 449-457.
            //
            //  Parameters:
            //
            //    Input, int NX, NY, values which determine the size
            //    of the matrix.
            //
            //    Input, int N, the number of rows and columns.
            //
            //    Input/output, int &SEED, the random number seed.
            //
            //    Output, double WATHEN_GB[(9*NX+13)*N], the matrix.
            //
        {
            double[] a;
            double[] em =
            {
                6.0, -6.0, 2.0, -8.0, 3.0, -8.0, 2.0, -6.0,
                -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0,
                2.0, -6.0, 6.0, -6.0, 2.0, -8.0, 3.0, -8.0,
                -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0,
                3.0, -8.0, 2.0, -6.0, 6.0, -6.0, 2.0, -8.0,
                -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0,
                2.0, -8.0, 3.0, -8.0, 2.0, -6.0, 6.0, -6.0,
                -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0
            };
            int i;
            int ii;
            int j;
            int jj;
            int kcol;
            int krow;
            int lda;
            int ml;
            int mu;
            int[] node = new int[8];
            double rho;

            ml = 3 * nx + 4;
            mu = 3 * nx + 4;
            lda = 2 * ml + mu + 1;
            a = new double[lda * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < lda; i++)
                {
                    a[i + j * lda] = 0.0;
                }
            }

            for (j = 0; j < nx; j++)
            {
                for (i = 0; i < nx; i++)
                {
                    node[0] = 3 * (j + 1) * nx + 2 * (j + 1) + 2 * (i + 1);
                    node[1] = node[0] - 1;
                    node[2] = node[0] - 2;
                    node[3] = (3 * (j + 1) - 1) * nx + 2 * (j + 1) + (i + 1) - 2;
                    node[4] = (3 * (j + 1) - 3) * nx + 2 * (j + 1) + 2 * (i + 1) - 4;
                    node[5] = node[4] + 1;
                    node[6] = node[4] + 2;
                    node[7] = node[3] + 1;

                    rho = 100.0 * UniformRNG.r8_uniform_01(ref seed);

                    for (krow = 0; krow < 8; krow++)
                    {
                        ii = node[krow];
                        for (kcol = 0; kcol < 8; kcol++)
                        {
                            jj = node[kcol];
                            a[ii - jj + ml + mu + jj * lda] = a[ii - jj + ml + mu + jj * lda]
                                                              + rho * em[krow + kcol * 8];
                        }
                    }
                }
            }

            return a;
        }

        public static double[] wathen_ge(int nx, int ny, int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WATHEN_GE returns the Wathen matrix as a general storage (GE) matrix.
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
            //    The matrix is the consistent mass matrix for a regular NX by NY grid
            //    of 8 node serendipity elements.
            //
            //    The local element numbering is
            //
            //      3--2--1
            //      |     |
            //      4     8
            //      |     |
            //      5--6--7
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
            //    The matrix is symmetric positive definite for any positive values of the
            //    density RHO(X,Y).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 July 2014
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
            //    Volume 7, Number 4, October 1987, pages 449-457.
            //
            //  Parameters:
            //
            //    Input, int NX, NY, values which determine the size 
            //    of the matrix.
            //
            //    Input, int N, the number of rows and columns.
            //
            //    Input/output, int &SEED, the random number seed.
            //
            //    Output, double WATHEN_GE[N*N], the matrix.
            //
        {
            double[] a;
            double[] em =
            {
                6.0, -6.0, 2.0, -8.0, 3.0, -8.0, 2.0, -6.0,
                -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0,
                2.0, -6.0, 6.0, -6.0, 2.0, -8.0, 3.0, -8.0,
                -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0,
                3.0, -8.0, 2.0, -6.0, 6.0, -6.0, 2.0, -8.0,
                -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0,
                2.0, -8.0, 3.0, -8.0, 2.0, -6.0, 6.0, -6.0,
                -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0
            };
            int i;
            int ii;
            int j;
            int jj;
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

            for (j = 0; j < nx; j++)
            {
                for (i = 0; i < nx; i++)
                {
                    node[0] = 3 * (j + 1) * nx + 2 * (j + 1) + 2 * (i + 1);
                    node[1] = node[0] - 1;
                    node[2] = node[0] - 2;
                    node[3] = (3 * (j + 1) - 1) * nx + 2 * (j + 1) + (i + 1) - 2;
                    node[4] = (3 * (j + 1) - 3) * nx + 2 * (j + 1) + 2 * (i + 1) - 4;
                    node[5] = node[4] + 1;
                    node[6] = node[4] + 2;
                    node[7] = node[3] + 1;

                    rho = 100.0 * UniformRNG.r8_uniform_01(ref seed);

                    for (krow = 0; krow < 8; krow++)
                    {
                        ii = node[krow];
                        for (kcol = 0; kcol < 8; kcol++)
                        {
                            jj = node[kcol];
                            a[ii + jj * n] = a[ii + jj * n] + rho * em[krow + kcol * 8];
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

        public static double[] wathen_st(int nx, int ny, int nz_num, ref int seed, int[] row,
                int[] col)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WATHEN_ST: Wathen matrix stored in sparse triplet (ST) format.
            //
            //  Discussion:
            //
            //    When dealing with sparse matrices in MATLAB, it can be much more efficient
            //    to work first with a triple of I, J, and X vectors, and only once
            //    they are complete, convert to MATLAB's sparse format.
            //
            //    The Wathen matrix is a finite element matrix which is sparse.
            //
            //    The entries of the matrix depend in part on a physical quantity
            //    related to density.  That density is here assigned random values between
            //    0 and 100.
            //
            //    The matrix order N is determined by the input quantities NX and NY,
            //    which would usually be the number of elements in the X and Y directions.
            //
            //    The value of N is
            //
            //      N = 3*NX*NY + 2*NX + 2*NY + 1,
            //
            //    The matrix is the consistent mass matrix for a regular NX by NY grid
            //    of 8 node serendipity elements.
            //
            //    The local element numbering is
            //
            //      3--2--1
            //      |     |
            //      4     8
            //      |     |
            //      5--6--7
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
            //    The matrix is symmetric positive definite for any positive values of the
            //    density RHO(X,Y).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 July 2014
            //
            //  Author:
            //
            //    John Burkardt.
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
            //    Volume 7, Number 4, October 1987, pages 449-457.
            //
            //  Parameters:
            //
            //    Input, int NX, NY, values which determine the size of 
            //    the matrix.
            //
            //    Input, int NZ_NUM, the number of values used to 
            //    describe the matrix.
            //
            //    Input/output, int &SEED, the random number seed.
            //
            //    Output, int ROW[NZ_NUM], COL[NZ_NUM], the row and 
            //    column indices of the nonzero entries.
            //
            //    Output, double WATHEN_ST[NZ_NUM], the nonzero entries of the matrix.
            //
        {
            double[] a;
            double[] em =
            {
                6.0, -6.0, 2.0, -8.0, 3.0, -8.0, 2.0, -6.0,
                -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0,
                2.0, -6.0, 6.0, -6.0, 2.0, -8.0, 3.0, -8.0,
                -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0,
                3.0, -8.0, 2.0, -6.0, 6.0, -6.0, 2.0, -8.0,
                -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0,
                2.0, -8.0, 3.0, -8.0, 2.0, -6.0, 6.0, -6.0,
                -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0
            };
            int i;
            int j;
            int k;
            int kcol;
            int krow;
            int[] node = new int[8];
            double rho;

            a = new double[nz_num];

            for (k = 0; k < nz_num; k++)
            {
                row[k] = 0;
                col[k] = 0;
                a[k] = 0.0;
            }

            k = 0;

            for (j = 0; j < nx; j++)
            {
                for (i = 0; i < nx; i++)
                {
                    node[0] = 3 * (j + 1) * nx + 2 * (j + 1) + 2 * (i + 1);
                    node[1] = node[0] - 1;
                    node[2] = node[0] - 2;
                    node[3] = (3 * (j + 1) - 1) * nx + 2 * (j + 1) + (i + 1) - 2;
                    node[4] = (3 * (j + 1) - 3) * nx + 2 * (j + 1) + 2 * (i + 1) - 4;
                    node[5] = node[4] + 1;
                    node[6] = node[4] + 2;
                    node[7] = node[3] + 1;

                    rho = 100.0 * UniformRNG.r8_uniform_01(ref seed);

                    for (krow = 0; krow < 8; krow++)
                    {
                        for (kcol = 0; kcol < 8; kcol++)
                        {
                            row[k] = node[krow];
                            col[k] = node[kcol];
                            a[k] = rho * em[krow + kcol * 8];
                            k = k + 1;
                        }
                    }
                }
            }

            return a;
        }

        public static int wathen_st_size(int nx, int ny)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WATHEN_ST_SIZE: Size of Wathen matrix stored in sparse triplet format.
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
            //    John Burkardt.
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
            //    Volume 7, Number 4, October 1987, pages 449-457.
            //
            //  Parameters:
            //
            //    Input, integer NX, NY, values which determine the size of the matrix.
            //
            //    Output, integer NZ_NUM, the number of items of data used to describe
            //    the matrix.
            //
        {
            int nz_num;

            nz_num = nx * ny * 64;

            return nz_num;
        }

        public static double[] wathen_xy(int nx, int ny, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    wathen_xy: coordinates of Wathen nodes.
            //
            //  Discussion:
            //
            //    We will take the region to be the unit square.
            //
            //    The grid uses quadratic serendipity elements.
            //
            //    Here is an illustration of the node numbering for NX = 3, NY = 2:
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
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 February 2020
            //
            //  Author:
            //
            //    John Burkardt.
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
            //    Volume 7, Number 4, October 1987, pages 449-457.
            //
            //  Input:
            //
            //    int NX, NY, values which determine the size of the matrix.
            //
            //    int N, the number of variables.
            //
            //  Output:
            //
            //    double wathen_xy[N*2], the X and Y coordinates of nodes.
            //
        {
            int i;
            int j;
            int k;

            double[] xy = new double[n * 2];

            k = 0;

            for (j = 0; j <= 2 * ny; j++)
            {
                if ((j % 2) == 0)
                {
                    for (i = 0; i <= 2 * nx; i++)
                    {
                        xy[k + i] = (double) (i) / (double) (2 * nx);
                        xy[k + i + n] = (double) (j) / (double) (2 * ny);
                    }

                    k = k + 2 * nx + 1;
                }
                else
                {
                    for (i = 0; i <= nx; i++)
                    {
                        xy[k + i] = (double) (i) / (double) (nx);
                        xy[k + i + n] = (double) (j) / (double) (2 * ny);
                    }

                    k = k + nx + 1;
                }
            }

            return xy;
        }

    }
}