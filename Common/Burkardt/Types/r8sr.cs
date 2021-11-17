using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8sr_dif2(int n, ref int nz, ref int[] row, ref int[] col, ref double[] diag,
            ref double[] off)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SR_DIF2 sets up an R8SR second difference matrix.
        //
        //  Discussion:
        //
        //    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
        //    The off-diagonal entries of row I are stored in entries ROW(I)
        //    through ROW(I+1)-1 of OFF.  COL(J) records the column index
        //    of the entry in A(J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int &NZ, the number of offdiagonal nonzero elements
        //    in the matrix.  NZ = 2 * N - 2.
        //
        //    Output, int ROW[N+1].  The nonzero offdiagonal elements 
        //    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
        //
        //    Output, int COL[NZ], contains the column index of the 
        //    element in the corresponding position in A.
        //
        //    Output, double DIAG[N], the diagonal elements of A.
        //
        //    Output, double OFF[NZ], the off-diagonal elements of A.
        //
    {
        int i;
        int nz2;

        for (i = 0; i < n; i++)
        {
            diag[i] = -2.0;
        }

        row[0] = 0;
        nz2 = 0;

        for (i = 0; i < n; i++)
        {
            switch (i)
            {
                case 0:
                    col[nz2] = i + 1;
                    off[nz2] = 1.0;
                    nz2 += 1;

                    row[i + 1] = row[i] + 1;
                    break;
                default:
                {
                    if (i < n - 1)
                    {
                        col[nz2] = i - 1;
                        off[nz2] = 1.0;
                        nz2 += 1;

                        col[nz2] = i + 1;
                        off[nz2] = 1.0;
                        nz2 += 1;

                        row[i + 1] = row[i] + 2;
                    }
                    else
                    {
                        col[nz2] = i - 1;
                        off[nz2] = 1.0;
                        nz2 += 1;

                        row[i + 1] = row[i] + 1;
                    }

                    break;
                }
            }
        }

    }

    public static void r8sr_indicator(int n, int nz, int[] row, int[] col, ref double[] diag,
            ref double[] off)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SR_INDICATOR sets up an R8SR indicator matrix.
        //
        //  Discussion:
        //
        //    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
        //    The off-diagonal entries of row I are stored in entries ROW(I)
        //    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
        //    the entry stored in OFF(J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the number of offdiagonal nonzero elements in A.
        //
        //    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
        //    are contained in A(ROW(I)) through A(ROW(I+1)-1).
        //
        //    Input, int COL[NZ], contains the column index of the element
        //    in the corresponding position in A.
        //
        //    Output, double DIAG[N], the diagonal elements of A.
        //
        //    Output, double OFF[NZ], the off-diagonal elements of A.
        //
    {
        int fac;
        int i;
        int j;
        int k;

        fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (i = 0; i < n; i++)
        {
            j = i;
            diag[i] = fac * (i + 1) + j + 1;

            for (k = row[i]; k <= row[i + 1] - 1; k++)
            {
                j = col[k];
                off[k] = fac * (i + 1) + j + 1;
            }
        }
    }

    public static double[] r8sr_mtv(int n, int nz, int[] row, int[] col, double[] diag,
            double[] off, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SR_MTV multiplies a vector times an R8SR matrix.
        //
        //  Discussion:
        //
        //    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
        //    The off-diagonal entries of row I are stored in entries ROW(I)
        //    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
        //    the entry stored in OFF(J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the number of offdiagonal nonzero elements in A.
        //
        //    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
        //    are contained in A(ROW(I)) through A(ROW(I+1)-1).
        //
        //    Input, int COL[NZ], contains the column index of the element
        //    in the corresponding position in A.
        //
        //    Input, double DIAG[N], the diagonal elements of A.
        //
        //    Input, double OFF[NZ], the off-diagonal elements of A.
        //
        //    Input, double X[N], the vector to be multiplies by A.
        //
        //    Output, double R8SR_MTV[N], the product A' * X.
        //
    {
        double[] b;
        int i;
        int j;
        int k;

        b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = diag[i] * x[i];
        }

        for (i = 0; i < n; i++)
        {
            for (k = row[i]; k <= row[i + 1] - 1; k++)
            {
                j = col[k];
                b[j] += off[k] * x[i];
            }
        }

        return b;
    }

    public static double[] r8sr_mv(int n, int nz, int[] row, int[] col, double[] diag,
            double[] off, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SR_MV multiplies an R8SR matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
        //    The off-diagonal entries of row I are stored in entries ROW(I)
        //    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
        //    the entry stored in OFF(J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the number of offdiagonal nonzero elements in A.
        //
        //    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
        //    are contained in A(ROW(I)) through A(ROW(I+1)-1).
        //
        //    Input, int COL[NZ], contains the column index of the element
        //    in the corresponding position in A.
        //
        //    Input, double DIAG[N], the diagonal elements of A.
        //
        //    Input, double OFF[NZ], the off-diagonal elements of A.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8SR_MV[N], the product A * X.
        //
    {
        double[] b;
        int i;
        int j;
        int k;

        b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = diag[i] * x[i];
        }

        for (i = 0; i < n; i++)
        {
            for (k = row[i]; k <= row[i + 1] - 1; k++)
            {
                j = col[k];
                b[i] += off[k] * x[j];
            }
        }

        return b;
    }

    public static void r8sr_print(int n, int nz, int[] row, int[] col, double[] diag,
            double[] off, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SR_PRINT prints an R8SR matrix.
        //
        //  Discussion:
        //
        //    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
        //    The off-diagonal entries of row I are stored in entries ROW(I)
        //    through ROW(I+1)-1 of OFF.  COL(J) records the column index
        //    of the entry in A(J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the number of offdiagonal nonzero elements in A.
        //
        //    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
        //    are contained in A(ROW(I)) through A(ROW(I+1)-1).
        //
        //    Input, int COL[NZ], contains the column index of the element
        //    in the corresponding position in A.
        //
        //    Input, double DIAG[N], the diagonal elements of A.
        //
        //    Input, double OFF[NZ], the off-diagonal elements of A.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8sr_print_some(n, nz, row, col, diag, off, 0, 0, n - 1, n - 1, title);
    }

    public static void r8sr_print_some(int n, int nz, int[] row, int[] col, double[] diag,
            double[] off, int ilo, int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SR_PRINT_SOME prints some of an R8SR matrix.
        //
        //  Discussion:
        //
        //    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
        //    The off-diagonal entries of row I are stored in entries ROW(I)
        //    through ROW(I+1)-1 of OFF.  COL(J) records the column index
        //    of the entry in A(J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the number of offdiagonal nonzero elements in A.
        //
        //    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
        //    are contained in A(ROW(I)) through A(ROW(I+1)-1).
        //
        //    Input, int COL[NZ], contains the column index of the element
        //    in the corresponding position in A.
        //
        //    Input, double DIAG[N], the diagonal elements of A.
        //
        //    Input, double OFF[NZ], the off-diagonal elements of A.
        //
        //    Input, int ILO, JLO, IHI, JHI, the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int INCX = 5;

        double aij;
        int i;
        int i2hi;
        int i2lo;
        int j;
        int j2hi;
        int j2lo;
        int k;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n - 1);
            j2hi = Math.Min(j2hi, jhi);

            Console.WriteLine("");
            cout = "  Col:  ";
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += j.ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);

            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            i2lo = Math.Max(ilo, 0);
            i2hi = Math.Min(ihi, n - 1);

            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString().PadLeft(6) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    aij = 0.0;
                    if (j == i)
                    {
                        aij = diag[i];
                    }
                    else
                    {
                        for (k = row[i]; k <= row[i + 1] - 1; k++)
                        {
                            if (j == col[k])
                            {
                                aij = off[k];
                            }
                        }
                    }

                    cout += aij.ToString().PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void r8sr_random(int n, int nz, int[] row, int[] col, ref double[] diag,
            ref double[] off, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SR_RANDOM randomizes an R8SR matrix.
        //
        //  Discussion:
        //
        //    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
        //    The off-diagonal entries of row I are stored in entries ROW(I)
        //    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
        //    the entry stored in OFF(J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the number of offdiagonal nonzero elements in A.
        //
        //    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
        //    are contained in A(ROW(I)) through A(ROW(I+1)-1).
        //
        //    Input, int COL[NZ], contains the column index of the element
        //    in the corresponding position in A.
        //
        //    Output, double DIAG[N], the diagonal elements of A.
        //
        //    Output, double OFF[NZ], the off-diagonal elements of A.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
    {
        int i;
        int j;

        for (i = 0; i < n; i++)
        {
            diag[i] = UniformRNG.r8_uniform_01(ref seed);
            for (j = row[i]; j <= row[i + 1] - 1; j++)
            {
                off[j] = UniformRNG.r8_uniform_01(ref seed);
            }
        }
    }

    public static double[] r8sr_to_r8ge(int n, int nz, int[] row, int[] col, double[] diag,
            double[] off)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SR_TO_R8GE converts an R8SR matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
        //    The off-diagonal entries of row I are stored in entries ROW(I)
        //    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
        //    the entry stored in OFF(J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the number of offdiagonal nonzero elements in A.
        //
        //    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
        //    are contained in A(ROW(I)) through A(ROW(I+1)-1).
        //
        //    Input, int COL[NZ], contains the column index of the element
        //    in the corresponding position in A.
        //
        //    Input, double DIAG[N], the diagonal elements of A.
        //
        //    Input, double OFF[NZ], the off-diagonal elements of A.
        //
        //    Output, double R8SR_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        double[] b;
        int i;
        int j;

        b = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                b[i + j * n] = 0.0;
            }
        }

        for (i = 0; i < n; i++)
        {
            b[i + i * n] = diag[i];
        }

        for (i = 0; i < n; i++)
        {
            for (j = row[i]; j <= row[i + 1] - 1; j++)
            {
                b[i + col[j] * n] = off[j];
            }
        }

        return b;
    }

    public static void r8sr_zeros(int n, int nz, int[] row, int[] col, ref double[] diag,
            ref double[] off)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SR_ZEROS zeros an R8SR matrix.
        //
        //  Discussion:
        //
        //    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
        //    The off-diagonal entries of row I are stored in entries ROW(I)
        //    through ROW(I+1)-1 of OFF.  COL(J) records the column index of the
        //    the entry stored in OFF(J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ, the number of offdiagonal nonzero elements in A.
        //
        //    Input, int ROW[N+1].  The nonzero offdiagonal elements of row I of A
        //    are contained in A(ROW(I)) through A(ROW(I+1)-1).
        //
        //    Input, int COL[NZ], contains the column index of the element
        //    in the corresponding position in A.
        //
        //    Output, double DIAG[N], the diagonal elements of A.
        //
        //    Output, double OFF[NZ], the off-diagonal elements of A.
        //
    {
        int i;
        int j;

        for (i = 0; i < n; i++)
        {
            diag[i] = 0.0;
            for (j = row[i]; j <= row[i + 1] - 1; j++)
            {
                off[j] = 0.0;
            }
        }

    }
}