using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8ss_dif2(int n, ref int na, ref int[] diag, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SS_DIF2 sets up an R8SS second difference matrix.
        //
        //  Discussion:
        //
        //    The R8SS storage format is used for real symmetric skyline matrices.
        //    This storage is appropriate when the nonzero entries of the
        //    matrix are generally close to the diagonal, but the number
        //    of nonzeroes above each diagonal varies in an irregular fashion.
        //
        //    In this case, the strategy is essentially to assign column J
        //    its own bandwidth, and store the strips of nonzeros one after
        //    another.   Note that what's important is the location of the
        //    furthest nonzero from the diagonal.  A slot will be set up for
        //    every entry between that and the diagonal, whether or not
        //    those entries are zero.
        //
        //    A skyline matrix can be Gauss-eliminated without disrupting
        //    the storage scheme, as long as no pivoting is required.
        //
        //    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
        //    although the actual storage needed will generally be about half of
        //    that.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int &NA, the dimension of the array A, which for
        //    this special case will 2*N-1.
        //
        //    Output, int DIAG[N], the indices in A of the N diagonal
        //    elements.
        //
        //    Output, double A[2*N-1], the R8SS matrix.
        //
    {
        int j;

        na = 0;

        for (j = 0; j < n; j++)
        {
            switch (j)
            {
                case > 0:
                    a[na] = -1.0;
                    na += 1;
                    break;
            }

            a[na] = 2.0;
            diag[j] = na;
            na += 1;
        }

    }

    public static bool r8ss_error(int[] diag, int n, int na)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SS_ERROR checks dimensions for an R8SS matrix.
        //
        //  Discussion:
        //
        //    The R8SS storage format is used for real symmetric skyline matrices.
        //    This storage is appropriate when the nonzero entries of the
        //    matrix are generally close to the diagonal, but the number
        //    of nonzeroes above each diagonal varies in an irregular fashion.
        //
        //    In this case, the strategy is essentially to assign column J
        //    its own bandwidth, and store the strips of nonzeros one after
        //    another.   Note that what's important is the location of the
        //    furthest nonzero from the diagonal.  A slot will be set up for
        //    every entry between that and the diagonal, whether or not
        //    those entries are zero.
        //
        //    A skyline matrix can be Gauss-eliminated without disrupting
        //    the storage format, as long as no pivoting is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIAG[N], the indices in A of the N diagonal elements.
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int NA, the dimension of the array A.
        //    NA must be at least N.
        //
        //    Output, bool R8SS_ERROR, is TRUE if an error was detected.
        //
    {
        int i;

        switch (n)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R8SS_ERROR - Illegal N = " + n + "");
                return true;
        }

        if (na < n)
        {
            Console.WriteLine("");
            Console.WriteLine("R8SS_ERROR - Illegal NA < N = " + n + "");
            return true;
        }

        if (diag[0] != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("R8SS_ERROR - DIAG[0] != 1.");
            return true;
        }

        for (i = 0; i < n - 1; i++)
        {
            if (diag[i + 1] <= diag[i])
            {
                Console.WriteLine("");
                Console.WriteLine("R8SS_ERROR - DIAG[I+1] <= DIAG[I].");
                return true;
            }
        }

        if (na < diag[n - 1])
        {
            Console.WriteLine("");
            Console.WriteLine("R8SS_ERROR - NA < DIAG[N-1].");
            return true;
        }

        return false;
    }

    public static double[] r8ss_indicator(int n, ref int na, ref int[] diag)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SS_INDICATOR sets up an R8SS indicator matrix.
        //
        //  Discussion:
        //
        //    The R8SS storage format is used for real symmetric skyline matrices.
        //    This storage is appropriate when the nonzero entries of the
        //    matrix are generally close to the diagonal, but the number
        //    of nonzeroes above each diagonal varies in an irregular fashion.
        //
        //    In this case, the strategy is essentially to assign column J
        //    its own bandwidth, and store the strips of nonzeros one after
        //    another.   Note that what's important is the location of the
        //    furthest nonzero from the diagonal.  A slot will be set up for
        //    every entry between that and the diagonal, whether or not
        //    those entries are zero.
        //
        //    A skyline matrix can be Gauss-eliminated without disrupting
        //    the storage format, as long as no pivoting is required.
        //
        //    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
        //    although the actual storage needed will generally be about half of
        //    that.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 July 2016
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
        //    Output, int &NA, the dimension of the array A, which for this
        //    special case will be the maximum, ( N * ( N + 1 ) ) / 2.
        //
        //    Output, int DIAG[N], the indices in A of the N diagonal elements.
        //
        //    Output, double R8SS_INDICATOR[(N*(N+1))/2], the R8SS matrix.
        //
    {
        int j;

        double[] a = r8vec_zeros_new(n * (n + 1) / 2);

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        na = 0;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i <= j; i++)
            {
                a[na] = fac * (i + 1) + j + 1;
                if (i == j)
                {
                    diag[j] = na;
                }

                na += 1;
            }
        }

        return a;
    }

    public static double[] r8ss_mv(int n, int na, int[] diag, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SS_MV multiplies an R8SS matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8SS storage format is used for real symmetric skyline matrices.
        //    This storage is appropriate when the nonzero entries of the
        //    matrix are generally close to the diagonal, but the number
        //    of nonzeroes above each diagonal varies in an irregular fashion.
        //
        //    In this case, the strategy is essentially to assign column J
        //    its own bandwidth, and store the strips of nonzeros one after
        //    another.   Note that what's important is the location of the
        //    furthest nonzero from the diagonal.  A slot will be set up for
        //    every entry between that and the diagonal, whether or not
        //    those entries are zero.
        //
        //    A skyline matrix can be Gauss-eliminated without disrupting
        //    the storage format, as long as no pivoting is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 July 2016
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
        //    Input, int NA, the dimension of the array A.
        //    NA must be at least N.
        //
        //    Input, int DIAG[N], the indices in A of the N diagonal elements.
        //
        //    Input, double A[NA], the R8SS matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8SS_MV[N], the product vector A*x.
        //
    {
        int j;

        double[] b = r8vec_zeros_new(n);

        int diagold = -1;
        int k = 0;

        for (j = 0; j < n; j++)
        {
            int ilo = j + 1 - (diag[j] - diagold);

            int i;
            for (i = ilo; i < j; i++)
            {
                b[i] += a[k] * x[j];
                b[j] += a[k] * x[i];
                k += 1;
            }

            b[j] += a[k] * x[j];
            k += 1;
            diagold = diag[j];
        }

        return b;
    }

    public static void r8ss_print(int n, int na, int[] diag, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SS_PRINT prints an R8SS matrix.
        //
        //  Discussion:
        //
        //    The R8SS storage format is used for real symmetric skyline matrices.
        //    This storage is appropriate when the nonzero entries of the
        //    matrix are generally close to the diagonal, but the number
        //    of nonzeroes above each diagonal varies in an irregular fashion.
        //
        //    In this case, the strategy is essentially to assign column J
        //    its own bandwidth, and store the strips of nonzeros one after
        //    another.   Note that what's important is the location of the
        //    furthest nonzero from the diagonal.  A slot will be set up for
        //    every entry between that and the diagonal, whether or not
        //    those entries are zero.
        //
        //    A skyline matrix can be Gauss-eliminated without disrupting
        //    the storage format, as long as no pivoting is required.
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
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NA, the dimension of the array A.
        //
        //    Input, int DIAG[N], the indices in A of the N diagonal elements.
        //
        //    Input, double A[NA], the R8SS matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8ss_print_some(n, na, diag, a, 0, 0, n - 1, n - 1, title);
    }

    public static void r8ss_print_some(int n, int na, int[] diag, double[] a, int ilo, int jlo,
            int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SS_PRINT_SOME prints some of an R8SS matrix.
        //
        //  Discussion:
        //
        //    The R8SS storage format is used for real symmetric skyline matrices.
        //    This storage is appropriate when the nonzero entries of the
        //    matrix are generally close to the diagonal, but the number
        //    of nonzeroes above each diagonal varies in an irregular fashion.
        //
        //    In this case, the strategy is essentially to assign column J
        //    its own bandwidth, and store the strips of nonzeros one after
        //    another.   Note that what's important is the location of the
        //    furthest nonzero from the diagonal.  A slot will be set up for
        //    every entry between that and the diagonal, whether or not
        //    those entries are zero.
        //
        //    A skyline matrix can be Gauss-eliminated without disrupting
        //    the storage format, as long as no pivoting is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 July 2016
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
        //    Input, int NA, the dimension of the array A.
        //
        //    Input, int DIAG[N], the indices in A of the N diagonal elements.
        //
        //    Input, double A[NA], the R8SS matrix.
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
            j2hi = Math.Min(j2hi, n - 1);
            j2hi = Math.Min(j2hi, jhi);

            Console.WriteLine("");
            string cout = "  Col: ";
            int j;
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
            int i2lo = Math.Max(ilo, 0);
            int i2hi = Math.Min(ihi, n - 1);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString().PadLeft(4) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    double aij = 0.0;

                    int ij;
                    int ijm1;
                    if (j < i)
                    {
                        ijm1 = i switch
                        {
                            0 => 0,
                            _ => diag[i - 1]
                        };

                        ij = diag[i];
                        if (ijm1 < ij + j - i)
                        {
                            aij = a[ij + j - i];
                        }
                    }
                    else if (j == i)
                    {
                        ij = diag[j];
                        aij = a[ij];
                    }
                    else if (i < j)
                    {
                        ijm1 = j switch
                        {
                            0 => 0,
                            _ => diag[j - 1]
                        };

                        ij = diag[j];
                        if (ijm1 < ij + i - j)
                        {
                            aij = a[ij + i - j];
                        }
                    }

                    cout += aij.ToString().PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void r8ss_random(int n, ref int na, ref int[] diag, ref double[] a, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SS_RANDOM randomizes an R8SS matrix.
        //
        //  Discussion:
        //
        //    The R8SS storage format is used for real symmetric skyline matrices.
        //    This storage is appropriate when the nonzero entries of the
        //    matrix are generally close to the diagonal, but the number
        //    of nonzeroes above each diagonal varies in an irregular fashion.
        //
        //    In this case, the strategy is essentially to assign column J
        //    its own bandwidth, and store the strips of nonzeros one after
        //    another.   Note that what's important is the location of the
        //    furthest nonzero from the diagonal.  A slot will be set up for
        //    every entry between that and the diagonal, whether or not
        //    those entries are zero.
        //
        //    A skyline matrix can be Gauss-eliminated without disrupting
        //    the storage format, as long as no pivoting is required.
        //
        //    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
        //    although the actual storage needed will generally be about half of
        //    that.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 July 2016
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
        //    Output, int *NA, the dimension of the array A.
        //    NA will be at least N and no greater than ( N * ( N + 1 ) ) / 2.
        //
        //    Output, int DIAG[N], the indices in A of the N diagonal elements.
        //
        //    Output, double A[(N*(N+1))/2], the R8SS matrix.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
    {
        int j;
        int k;

        na = 0;
        //
        //  Set the values of DIAG.
        //
        diag[0] = 0;
        na = 1;

        for (j = 1; j < n; j++)
        {
            k = UniformRNG.i4_uniform_ab(1, j + 1, ref seed);
            diag[j] = diag[j - 1] + k;
            na += k;
        }

        //
        //  Now set the values of A.
        //
        int diagold = -1;
        k = 0;

        for (j = 0; j < n; j++)
        {
            int ilo = j + 1 - (diag[j] - diagold);

            int i;
            for (i = ilo; i <= j; i++)
            {
                a[k] = UniformRNG.r8_uniform_01(ref seed);
                k += 1;
            }

            diagold = diag[j];
        }
    }

    public static double[] r8ss_to_r8ge(int n, int na, int[] diag, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SS_TO_R8GE copies an R8SS matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8SS storage format is used for real symmetric skyline matrices.
        //    This storage is appropriate when the nonzero entries of the
        //    matrix are generally close to the diagonal, but the number
        //    of nonzeroes above each diagonal varies in an irregular fashion.
        //
        //    In this case, the strategy is essentially to assign column J
        //    its own bandwidth, and store the strips of nonzeros one after
        //    another.   Note that what's important is the location of the
        //    furthest nonzero from the diagonal.  A slot will be set up for
        //    every entry between that and the diagonal, whether or not
        //    those entries are zero.
        //
        //    A skyline matrix can be Gauss-eliminated without disrupting
        //    the storage format, as long as no pivoting is required.
        //
        //  Example:
        //
        //    11   0  13  0 15
        //     0  22  23  0  0
        //    31  32  33 34  0
        //     0   0  43 44  0
        //    51   0   0  0 55
        //
        //    A = ( 11 | 22 | 13, 23, 33 | 34, 44 | 15, 0, 0, 0, 55 )
        //    NA = 12
        //    DIAG = ( 0, 1, 4, 6, 11 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 July 2016
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
        //    Input, int NA, the dimension of the array A.
        //    NA must be at least N.
        //
        //    Input, int DIAG[N], the indices in A of the N diagonal elements.
        //
        //    Input, double A[NA], the R8SS matrix.
        //
        //    Output, double R8SS_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int j;

        double[] b = r8vec_zeros_new(n * n);

        int diagold = -1;
        int k = 0;

        for (j = 0; j < n; j++)
        {
            int ilo = j + 1 - (diag[j] - diagold);

            int i;
            for (i = ilo; i < j; i++)
            {
                b[i + j * n] = a[k];
                b[j + i * n] = a[k];
                k += 1;
            }

            b[j + j * n] = a[k];
            k += 1;

            diagold = diag[j];
        }

        return b;
    }

    public static double[] r8ss_zeros(int n, int na, ref int[] diag)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SS_ZEROS zeros an R8SS matrix.
        //
        //  Discussion:
        //
        //    The R8SS storage format is used for real symmetric skyline matrices.
        //    This storage is appropriate when the nonzero entries of the
        //    matrix are generally close to the diagonal, but the number
        //    of nonzeroes above each diagonal varies in an irregular fashion.
        //
        //    In this case, the strategy is essentially to assign column J
        //    its own bandwidth, and store the strips of nonzeros one after
        //    another.   Note that what's important is the location of the
        //    furthest nonzero from the diagonal.  A slot will be set up for
        //    every entry between that and the diagonal, whether or not
        //    those entries are zero.
        //
        //    A skyline matrix can be Gauss-eliminated without disrupting
        //    the storage format, as long as no pivoting is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NA, the dimension of the array A.
        //
        //    Output, int DIAG[N], the indices in A of the N diagonal elements.
        //
        //    Output, double R8SS_ZERO[NA], the R8SS matrix.
        //
    {
        int i;

        double[] a = r8vec_zeros_new(na);

        int k = -1;
        for (i = 0; i < n; i++)
        {
            k = k + i + 1;
            diag[i] = k;
        }

        return a;
    }

}