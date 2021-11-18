using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8gd_dif2(int n, int ndiag, int[] offset)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GD_DIF2 sets up an R8GD second difference matrix.
        //
        //  Discussion:
        //
        //    The R8GD storage format is suitable for matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0.
        //    Each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
        //
        //    Now, assuming that only a few of these diagonals contain nonzeros,
        //    then for the I-th diagonal to be saved, we stored its offset in
        //    OFFSET(I), and its entries in column I of the matrix.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 3.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal
        //    storage.  The values -1, 0 and +1 should be included.
        //
        //    Output, double R8GD_DIF2[N*NDIAG], the R8GD matrix.
        //
    {
        int i;

        switch (ndiag)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("R8GD_DIF2 - Fatal error!");
                Console.WriteLine("  NDIAG must be at least 3.");
                return null;
        }

        double[] a = r8vec_zeros_new(n * ndiag);

        for (i = 0; i < n; i++)
        {
            int jdiag;
            for (jdiag = 0; jdiag < ndiag; jdiag++)
            {
                int j = i + offset[jdiag];
                switch (j)
                {
                    case >= 0 when j < n:
                        switch (offset[jdiag])
                        {
                            case 0:
                                a[i + jdiag * n] = 2.0;
                                break;
                            case -1:
                            case +1:
                                a[i + jdiag * n] = -1.0;
                                break;
                        }

                        break;
                }
            }
        }

        return a;
    }

    public static double[] r8gd_indicator(int n, int ndiag, int[] offset)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GD_INDICATOR sets up an R8GD indicator matrix.
        //
        //  Discussion:
        //
        //    The R8GD storage format is used for matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0.
        //    Each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
        //
        //    Now, assuming that only a few of these diagonals contain nonzeros,
        //    then for the I-th diagonal to be saved, we stored its offset in
        //    OFFSET(I), and its entries in column I of the matrix.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 1, and no more than 2 * N - 1.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Output, double R8GD_INDICATOR[N*NDIAG], the R8GD matrix.
        //
    {
        int i;

        double[] a = r8vec_zeros_new(n * ndiag);

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (i = 1; i <= n; i++)
        {
            int diag;
            for (diag = 1; diag <= ndiag; diag++)
            {
                int j = i + offset[diag - 1];
                a[i - 1 + (diag - 1) * n] = j switch
                {
                    >= 1 when j <= n => fac * i + j,
                    _ => 0.0
                };
            }
        }

        return a;
    }

    public static double[] r8gd_mtv(int n, int ndiag, int[] offset, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GD_MTV multiplies a vector by an R8GD matrix.
        //
        //  Discussion:
        //
        //    The R8GD storage format is used for matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0.
        //    Each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
        //
        //    Now, assuming that only a few of these diagonals contain nonzeros,
        //    then for the I-th diagonal to be saved, we stored its offset in
        //    OFFSET(I), and its entries in column I of the matrix.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 1, and no more than 2 * N - 1.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the R8GD matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8GD_MTV[N], the product X*A.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            int diag;
            for (diag = 0; diag < ndiag; diag++)
            {
                int j = i + offset[diag];
                switch (j)
                {
                    case >= 0 when j < n:
                        b[j] += x[i] * a[i + diag * n];
                        break;
                }
            }
        }

        return b;
    }

    public static double[] r8gd_mv(int n, int ndiag, int[] offset, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GD_MV multiplies an R8GD matrix by a vector.
        //
        //  Discussion:
        //
        //    The R8GD storage format is used for matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0.
        //    Each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
        //
        //    Now, assuming that only a few of these diagonals contain nonzeros,
        //    then for the I-th diagonal to be saved, we stored its offset in
        //    OFFSET(I), and its entries in column I of the matrix.  
        //
        //  Example:
        //
        //    The "offset" value is printed near the first entry of each diagonal
        //    in the original matrix, and above the columns in the new matrix.
        //
        //    Original matrix               New Matrix
        //
        //      0    1   2   3   4   5        -3  -2   0   1   3   5
        //                             
        //        11  12   0  14   0  16      --  --  11  12  14  16
        //   -1 =  0  22  23   0  25   0      --  --  22  23  25  --
        //   -2 = 31   0  33  34   0  36      --  31  33  34  36  --
        //   -3 = 41  42   0  44  45   0      41  42  44  45  --  --
        //   -4 =  0  52  53   0  55  56      52  53  55  56  --  --
        //   -5 =  0   0  63  64  65  66      63  64  66  --  --  --
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 1, and no more than 2 * N - 1.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the R8GD matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8GD_MV[N], the product A * x.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            int diag;
            for (diag = 0; diag < ndiag; diag++)
            {
                int j = i + offset[diag];
                switch (j)
                {
                    case >= 0 when j < n:
                        b[i] += a[i + diag * n] * x[j];
                        break;
                }
            }
        }

        return b;
    }

    public static void r8gd_print(int n, int ndiag, int[] offset, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GD_PRINT prints an R8GD matrix.
        //
        //  Discussion:
        //
        //    The R8GD storage format is used for matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0.
        //    Each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
        //
        //    Now, assuming that only a few of these diagonals contain nonzeros,
        //    then for the I-th diagonal to be saved, we stored its offset in
        //    OFFSET(I), and its entries in column I of the matrix.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 April 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 1, and no more than 2 * N - 1.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the R8GD matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8gd_print_some(n, ndiag, offset, a, 1, 1, n, n, title);
    }

    public static void r8gd_print_some(int n, int ndiag, int[] offset, double[] a, int ilo,
            int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GD_PRINT_SOME prints some of an R8GD matrix.
        //
        //  Discussion:
        //
        //    The R8GD storage format is used for matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0.
        //    Each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
        //
        //    Now, assuming that only a few of these diagonals contain nonzeros,
        //    then for the I-th diagonal to be saved, we stored its offset in
        //    OFFSET(I), and its entries in column I of the matrix.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 April 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 1, and no more than 2 * N - 1.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the R8GD matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int j2lo;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            int j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n);
            j2hi = Math.Min(j2hi, jhi);

            Console.WriteLine("");
            string cout = "  Col:  ";
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += j.ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 1);
            int i2hi = Math.Min(ihi, n);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                int j2;
                for (j2 = j2lo; j2 <= j2hi; j2++)
                {
                    double aij = 0.0;
                    int diag;
                    for (diag = 0; diag < ndiag; diag++)
                    {
                        if (j2 - i == offset[diag])
                        {
                            aij = a[i - 1 + diag * n];
                        }
                    }

                    cout += aij.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8gd_random(int n, int ndiag, int[] offset, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GD_RANDOM randomizes an R8GD matrix.
        //
        //  Discussion:
        //
        //    The R8GD storage format is used for matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0.
        //    Each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
        //
        //    Now, assuming that only a few of these diagonals contain nonzeros,
        //    then for the I-th diagonal to be saved, we stored its offset in
        //    OFFSET(I), and its entries in column I of the matrix.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 1, and no more than 2 * N - 1.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8GD_RANDOM[N*NDIAG], the R8GD matrix.
        //
    {
        int i;

        double[] a = r8vec_zeros_new(n * ndiag);

        for (i = 1; i <= n; i++)
        {
            int diag;
            for (diag = 0; diag < ndiag; diag++)
            {
                int j = i + offset[diag];
                a[i - 1 + diag * n] = j switch
                {
                    >= 1 when j <= n => UniformRNG.r8_uniform_01(ref seed),
                    _ => a[i - 1 + diag * n]
                };
            }
        }

        return a;
    }

    public static double[] r8gd_to_r8ge(int n, int ndiag, int[] offset, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GD_TO_R8GE copies an R8GD matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GD storage format is used for matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0.
        //    Each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
        //
        //    Now, assuming that only a few of these diagonals contain nonzeros,
        //    then for the I-th diagonal to be saved, we stored its offset in
        //    OFFSET(I), and its entries in column I of the matrix.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 1, and no more than 2 * N - 1.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the R8GD matrix.
        //
        //    Output, double R8GD_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n * n);

        for (i = 0; i < n; i++)
        {
            int diag;
            for (diag = 0; diag < ndiag; diag++)
            {
                int j = i + offset[diag];
                b[i + j * n] = j switch
                {
                    >= 0 when j <= n - 1 => a[i + diag * n],
                    _ => b[i + j * n]
                };
            }
        }

        return b;
    }

    public static double[] r8gd_zeros(int n, int ndiag, int[] offset)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GD_ZEROS zeros an R8GD matrix.
        //
        //  Discussion:
        //
        //    The R8GD storage format is used for matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0.
        //    Each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
        //
        //    Now, assuming that only a few of these diagonals contain nonzeros,
        //    then for the I-th diagonal to be saved, we stored its offset in
        //    OFFSET(I), and its entries in column I of the matrix.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 1, and no more than 2 * N - 1.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Output, double R8GD_ZERO[N*NDIAG], the R8GD matrix.
        //
    {
        double[] a = r8vec_zeros_new(n * ndiag);

        return a;
    }
}