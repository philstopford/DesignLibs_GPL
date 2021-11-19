using System;
using System.Globalization;
using System.IO;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Storage;

public static class SparseTriplet
{
    public static void st_data_read(string input_filename, int m, int n, int nst, ref int[] ist,
            ref int[] jst, ref double[] ast)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ST_DATA_READ reads the data of an ST file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILENAME, the name of the ST file.
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, int NST, the number of nonzeros.
        //
        //    Output, int IST[NNZERO], JST[NNZERO], the row and column indices.
        //
        //    Output, double AST[NST], the nonzero values.
        //
    {
        int k = 0;
        try
        {
            string[] input = File.ReadAllLines(input_filename);

            for (k = 0; k < nst; k++)
            {
                string[] tokens = Helpers.splitStringByWhitespace(input[k]);

                i4 ti = typeMethods.s_to_i4(tokens[0]);
                i4 tj = typeMethods.s_to_i4(tokens[1]);
                r8 taij = typeMethods.s_to_r8(tokens[2]);

                ist[k] = ti.val;
                jst[k] = tj.val;
                ast[k] = taij.val;
            }
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("ST_DATA_READ - Fatal error!");
            Console.WriteLine("  I/O error reading data index " + k + "");
        }

    }

    public static void st_header_print(int i_min, int i_max, int j_min, int j_max, int m,
            int n, int nst)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ST_HEADER_PRINT prints the header of an ST file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I_MAX, the minimum and maximum row indices.
        //
        //    Input, int J_MIN, J_MAX, the minimum and maximum column indices.
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, int NST, the number of nonzeros.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("  ST header information:");
        Console.WriteLine("");
        Console.WriteLine("  Minimum row index I_MIN = " + i_min + "");
        Console.WriteLine("  Maximum row index I_MAX = " + i_max + "");
        Console.WriteLine("  Minimum col index J_MIN = " + j_min + "");
        Console.WriteLine("  Maximum col index J_MAX = " + j_max + "");
        Console.WriteLine("  Number of rows        M = " + m + "");
        Console.WriteLine("  Number of columns     N = " + n + "");
        Console.WriteLine("  Number of nonzeros  NST = " + nst + "");
    }

    public static void st_header_read(string input_filename, ref int i_min, ref int i_max, ref int j_min,
            ref int j_max, ref int m, ref int n, ref int nst)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ST_HEADER_READ reads the header of an ST file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILENAME, the name of the ST file.
        //
        //    Output, int &I_MIN, &I_MAX, the minimum and maximum row indices.
        //
        //    Output, int &J_MIN, &J_MAX, the minimum and maximum column indices.
        //
        //    Output, int &M, the number of rows.
        //
        //    Output, int &N, the number of columns .
        //
        //    Output, int &NST, the number of nonzeros.
        //
    {

        nst = 0;
        i_min = +typeMethods.i4_huge();
        i_max = -typeMethods.i4_huge();
        j_min = +typeMethods.i4_huge();
        j_max = -typeMethods.i4_huge();

        try
        {
            string[] input = File.ReadAllLines(input_filename);

            for (;;)
            {
                string[] tokens = Helpers.splitStringByWhitespace(input[nst]);

                i4 ti = typeMethods.s_to_i4(tokens[0]);
                i4 tj = typeMethods.s_to_i4(tokens[1]);

                int i = ti.val;
                int j = tj.val;

                nst += 1;
                i_min = Math.Min(i_min, i);
                i_max = Math.Max(i_max, i);
                j_min = Math.Min(j_min, j);
                j_max = Math.Max(j_max, j);

                if (nst == input.Length)
                {
                    break;
                }
            }
        }
        catch (Exception)
        {
        }

        m = i_max - i_min + 1;
        n = j_max - j_min + 1;
    }

    public static double[] st_mv(int m, int n, int nst, int[] ist, int[] jst, double[] ast,
            double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ST_MV multiplies an ST matrix by an R8VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, int NST, the number of nonzeros.
        //
        //    Input, int IST[NST], JST[NST], the row and column indices.
        //
        //    Input, double AST[NST], the nonzero values.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double B[M], the product vector A*X.
        //
    {
        double[] b = typeMethods.r8vec_zero_new(m);

        for (int k = 0; k < nst; k++)
        {
            int i = ist[k];
            int j = jst[k];
            b[i] += ast[k] * x[j];
        }

        return b;
    }

    public static void st_print(int m, int n, int nst, int[] ist, int[] jst, double[] ast,
            string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ST_PRINT prints a sparse matrix in ST format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, int NST, the number of ST elements.
        //
        //    Input, int IST[NST], JST[NST], the ST rows and columns.
        //
        //    Input, double AST[NST], the ST values.
        //
        //    Input, string TITLE, a title.
        //
    {
        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("     #     I     J       A");
        Console.WriteLine("  ----  ----  ----  --------------");
        Console.WriteLine("");
        for (int k = 0; k < nst; k++)
        {
            Console.WriteLine(k.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + ist[k].ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + jst[k].ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + ast[k].ToString(CultureInfo.InvariantCulture).PadLeft(16) + "");
        }
    }

    public static int st_to_ccs_size(int nst, int[] ist, int[] jst)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ST_TO_ccs_SIZE sizes CCS indexes based on ST data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NST, the number of ST elements.
        //
        //    Input, int IST[NST], JST[NST], the ST rows and columns.
        //
        //    Output, int ST_TO_ccs_SIZE, the number of CCS elements.
        //
    {
        //
        //  Make copies so the sorting doesn't confuse the user.
        //
        int[] ist2 = typeMethods.i4vec_copy_new(nst, ist);
        int[] jst2 = typeMethods.i4vec_copy_new(nst, jst);
        //
        //  Sort by column first, then row.
        //
        typeMethods.i4vec2_sort_a(nst, ref jst2, ref ist2);
        //
        //  Count the unique pairs.
        //
        int ncc = typeMethods.i4vec2_sorted_unique_count(nst, jst2, ist2);

        return ncc;
    }

    public static void st_to_ccs_index(int nst, int[] ist, int[] jst, int ncc, int n,
            ref int[] icc, ref int[] ccc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ST_TO_ccs_INDEX creates CCS indices from ST data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NST, the number of ST elements.
        //
        //    Input, int IST[NST], JST[NST], the ST rows and columns.
        //
        //    Input, int NCC, the number of CCS elements.
        //
        //    Input, int N, the number of columns in the matrix.
        //
        //    Output, int ICC[NCC], the CCS rows.
        //
        //    Output, int CCC[N+1], the compressed CCS columns.
        //
    {
        int j;
        int jhi;
        int k;
        //
        //  Make copies so the sorting doesn't confuse the user.
        //
        int[] ist2 = typeMethods.i4vec_copy_new(nst, ist);
        int[] jst2 = typeMethods.i4vec_copy_new(nst, jst);
        //
        //  Sort the elements.
        //
        typeMethods.i4vec2_sort_a(nst, ref jst2, ref ist2);
        //
        //  Get the unique elements.
        //
        int[] jcc = new int[ncc];
        typeMethods.i4vec2_sorted_uniquely(nst, jst2, ist2, ncc, ref jcc, ref icc);
        //
        //  Compress the column index.
        //
        ccc[0] = 0;
        int jlo = 0;
        for (k = 0; k < ncc; k++)
        {
            jhi = jcc[k];
            if (jhi == jlo)
            {
                continue;
            }

            for (j = jlo + 1; j <= jhi; j++)
            {
                ccc[j] = k;
            }

            jlo = jhi;
        }

        jhi = n;
        for (j = jlo + 1; j <= jhi; j++)
        {
            ccc[j] = ncc;
        }
    }

    public static double[] st_to_ccs_values(int nst, int[] ist, int[] jst, double[] ast, int ncc,
            int n, int[] icc, int[] ccc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ST_TO_ccs_VALUES creates CCS values from ST data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NST, the number of ST elements.
        //
        //    Input, int IST[NST], JST[NST], the ST rows and columns.
        //
        //    Input, double AST[NST], the ST values.
        //
        //    Input, int NCC, the number of CCS elements.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, int ICC[NCC], the CCS rows.
        //
        //    Input, int CCC[N+1], the CCS compressed columns.
        //
        //    Output, double ST_TO_ccs_VALUES[NCC], the CCS values.
        //
    {
        int i;
        int kst;

        double[] acc = new double[ncc];

        for (i = 0; i < ncc; i++)
        {
            acc[i] = 0.0;
        }

        for (kst = 0; kst < nst; kst++)
        {
            i = ist[kst];
            int j = jst[kst];

            int clo = ccc[j];
            int chi = ccc[j + 1];

            bool fail = true;

            int kcc;
            for (kcc = clo; kcc < chi; kcc++)
            {
                if (icc[kcc] != i)
                {
                    continue;
                }

                acc[kcc] += ast[kst];
                fail = false;
                break;
            }

            switch (fail)
            {
                case true:
                    Console.WriteLine("");
                    Console.WriteLine("ST_TO_ccs_VALUES - Fatal error!");
                    Console.WriteLine("  ST entry cannot be located in CCS array.");
                    Console.WriteLine("  ST index KST    = " + kst + "");
                    Console.WriteLine("  ST row IST(KST) = " + ist[kst] + "");
                    Console.WriteLine("  ST col JST(KST) = " + jst[kst] + "");
                    Console.WriteLine("  ST val AST(KST) = " + ast[kst] + "");
                    throw new Exception();
            }

        }

        return acc;
    }

    public static double[] wathen_st(int nx, int ny, int nz_num, ref int seed, ref int[] row,
            ref int[] col )

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
        int j;
        int k;
        int[] node = new int[8];

        double[] a = new double[nz_num];

        for (k = 0; k < nz_num; k++)
        {
            row[k] = 0;
            col[k] = 0;
            a[k] = 0.0;
        }

        k = 0;

        for (j = 0; j < nx; j++)
        {
            int i;
            for (i = 0; i < nx; i++)
            {
                node[0] = 3 * (j + 1) * nx + 2 * (j + 1) + 2 * (i + 1);
                node[1] = node[0] - 1;
                node[2] = node[0] - 2;
                node[3] = (3 * (j + 1) - 1) * nx + 2 * (j + 1) + i + 1 - 2;
                node[4] = (3 * (j + 1) - 3) * nx + 2 * (j + 1) + 2 * (i + 1) - 4;
                node[5] = node[4] + 1;
                node[6] = node[4] + 2;
                node[7] = node[3] + 1;

                double rho = 100.0 * UniformRNG.r8_uniform_01(ref seed);

                int krow;
                for (krow = 0; krow < 8; krow++)
                {
                    int kcol;
                    for (kcol = 0; kcol < 8; kcol++)
                    {
                        row[k] = node[krow];
                        col[k] = node[kcol];
                        a[k] = rho * em[krow + kcol * 8];
                        k += 1;
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
        int nz_num = nx * ny * 64;

        return nz_num;
    }

}