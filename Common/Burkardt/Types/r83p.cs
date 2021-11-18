using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{

    public static double r83p_det(int n, double[] a_lu, double work4)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_DET computes the determinant of a matrix factored by R83P_FA.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored 
        //    as a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Input, double A_LU[3*N], the LU factors from R83P_FA.
        //
        //    Input, double WORK4, factorization information from R83P_FA.
        //
        //    Output, double R83P_DET, the determinant of the matrix.
        //
    {
        int i;

        double det = work4;
        for (i = 0; i <= n - 2; i++)
        {
            det *= a_lu[1 + i * 3];
        }

        return det;
    }

    public static int r83p_fa(int n, ref double[] a, ref double[] work2, ref double[] work3, ref double work4)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_FA factors an R83P matrix.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored as 
        //    a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //    Once the matrix has been factored by R83P_FA, R83P_SL may be called
        //    to solve linear systems involving the matrix.
        //
        //    The logical matrix has a form which is suggested by this diagram:
        //
        //      D1 U1          L1
        //      L2 D2 U2
        //         L3 R83 U3
        //            L4 D4 U4
        //               L5 R85 U5
        //      U6          L6 D6
        //
        //    The algorithm treats the matrix as a border banded matrix:
        //
        //      ( A1  A2 )
        //      ( A3  A4 )
        //
        //    where:
        //
        //      D1 U1          | L1
        //      L2 D2 U2       |  0
        //         L3 R83 U3    |  0
        //            L4 D4 U4 |  0
        //               L5 R85 | U5
        //      ---------------+---
        //      U6  0  0  0 L6 | D6
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
        //
        //  Method:
        //
        //    The algorithm rewrites the system as:
        //
        //         X1 + inverse(A1) A2 X2 = inverse(A1) B1
        //
        //      A3 X1 +             A4 X2 = B2
        //
        //    The first equation can be "solved" for X1 in terms of X2:
        //
        //         X1 = - inverse(A1) A2 X2 + inverse(A1) B1
        //
        //    allowing us to rewrite the second equation for X2 explicitly:
        //
        //      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Input/output, double A[3*N].
        //    On input, the periodic tridiagonal matrix.  
        //    On output, the arrays have been modified to hold information
        //    defining the border-banded factorization of submatrices A1
        //    and A3.
        //
        //    Output, int R83P_FA, singularity flag.
        //    0, no singularity detected.
        //    nonzero, the factorization failed on the INFO-th step.
        //
        //    Output, double WORK2[N-1], WORK3[N-1], *WORK4, factorization information.
        //
    {
        int i;

        double[] work1 = new double[n - 1];
        //
        //  Compute inverse(A1):
        //
        int info = r83_np_fa(n - 1, ref a);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("R83P_FA - Fatal error!");
            Console.WriteLine("  R83_NP_FA returned INFO = " + info + "");
            Console.WriteLine("  Factoring failed for column INFO.");
            Console.WriteLine("  The tridiagonal matrix A1 is singular.");
            Console.WriteLine("  This algorithm cannot continue!");
            return 1;
        }

        //
        //  WORK2 := inverse(A1) * A2.
        //
        work2[0] = a[2 + (n - 1) * 3];
        for (i = 1; i < n - 2; i++)
        {
            work2[i] = 0.0;
        }

        work2[n - 2] = a[0 + (n - 1) * 3];

        int job = 0;
        work1 = r83_np_sl(n - 1, a, work2, job);
        for (i = 0; i < n - 1; i++)
        {
            work2[i] = work1[i];
        }

        //
        //  WORK3 := inverse ( A1' ) * A3'.
        //
        work3[0] = a[0 + 0 * 3];
        for (i = 1; i < n - 2; i++)
        {
            work3[i] = 0.0;
        }

        work3[n - 2] = a[2 + (n - 2) * 3];

        job = 1;
        work1 = r83_np_sl(n - 1, a, work3, job);
        for (i = 0; i < n - 1; i++)
        {
            work3[i] = work1[i];
        }

        //
        //  A4 := ( A4 - A3 * inverse(A1) * A2 )
        //
        work4 = a[1 + (n - 1) * 3] - a[0 + 0 * 3] * work2[0] - a[2 + (n - 2) * 3] * work2[n - 2];

        switch (work4)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R83P_FA - Fatal error!");
                Console.WriteLine("  The factored A4 submatrix is zero.");
                Console.WriteLine("  This algorithm cannot continue!");
                return 1;
            default:
                return 0;
        }
    }

    public static double[] r83p_indicator(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_INDICATOR sets up an R83P indicator matrix.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored 
        //    as a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
        //
        //    Here are the values as stored in an indicator matrix:
        //
        //      51 12 23 34 45
        //      11 22 33 44 55
        //      21 32 43 54 15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 2.
        //
        //    Output, double R83P_INDICATOR[3*N], the R83P indicator matrix.
        //
    {
        double[] a = new double[3 * n];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        int i = n;
        int j = 1;
        a[0] = fac * i + j;
        for (j = 2; j <= n; j++)
        {
            i = j - 1;
            a[0 + (j - 1) * 3] = fac * i + j;
        }

        for (j = 1; j <= n; j++)
        {
            i = j;
            a[1 + (j - 1) * 3] = fac * i + j;
        }

        for (j = 1; j <= n - 1; j++)
        {
            i = j + 1;
            a[2 + (j - 1) * 3] = fac * i + j;
        }

        j = n;
        a[2 + (j - 1) * 3] = fac + j;

        return a;
    }

    public static double[] r83p_ml(int n, double[] a_lu, double[] x, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_ML computes A * x or x * A, where A has been factored by R83P_FA.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored 
        //    as a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Input, double A_LU[3*N], the LU factors from R83P_FA.
        //
        //    Input, double X[N], the vector to be multiplied by the matrix.
        //
        //    Input, int JOB, indicates what product should be computed.
        //    0, compute A * x.
        //    nonzero, compute A' * x.
        //
        //    Output, double R83P_ML[N], the result of the multiplication.
        //
    {
        int i;
        //
        //  Multiply A(1:N-1,1:N-1) and X(1:N-1).
        //
        double[] b_short = r83_np_ml(n - 1, a_lu, x, job);

        double[] b = new double[n];

        for (i = 0; i < n - 1; i++)
        {
            b[i] = b_short[i];
        }

        b[n - 1] = 0.0;

        switch (job)
        {
            //
            //  Add terms from the border.
            //
            case 0:
                b[0] += a_lu[2 + (n - 1) * 3] * x[n - 1];
                b[n - 2] += a_lu[0 + (n - 1) * 3] * x[n - 1];
                b[n - 1] = a_lu[0 + 0 * 3] * x[0] + a_lu[2 + (n - 2) * 3] * x[n - 2]
                                                  + a_lu[1 + (n - 1) * 3] * x[n - 1];
                break;
            default:
                b[0] += a_lu[0 + 0 * 3] * x[n - 1];
                b[n - 2] += a_lu[2 + (n - 2) * 3] * x[n - 1];
                b[n - 1] = a_lu[2 + (n - 1) * 3] * x[0] + a_lu[0 + (n - 1) * 3] * x[n - 2]
                                                        + a_lu[1 + (n - 1) * 3] * x[n - 1];
                break;
        }

        return b;
    }

    public static double[] r83p_mtv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_MTV multiplies a vector times an R83P matrix.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored as 
        //    a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Input, double A[3*N], the R83P matrix.
        //
        //    Input, double X, the vector to be multiplied by A.
        //
        //    Output, double R83P_MTV[N], the product X * A.
        //
    {
        int i;

        double[] b = new double[n];

        b[0] = a[0 + 0 * 3] * x[n - 1] + a[1 + 0 * 3] * x[0] + a[2 + 0 * 3] * x[1];

        for (i = 2; i <= n - 1; i++)
        {
            b[i - 1] = a[0 + (i - 1) * 3] * x[i - 2] + a[1 + (i - 1) * 3] * x[i - 1] + a[2 + (i - 1) * 3] * x[i];
        }

        b[n - 1] = a[0 + (n - 1) * 3] * x[n - 2] + a[1 + (n - 1) * 3] * x[n - 1] + a[2 + (n - 1) * 3] * x[0];

        return b;
    }

    public static double[] r83p_mv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_MV multiplies an R83P matrix times a vector.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored as 
        //    a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Input, double A[3*N], the R83P matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R83P_MV[N], the product A * x.
        //
    {
        int i;

        double[] b = new double[n];

        b[0] = a[2 + (n - 1) * 3] * x[n - 1] + a[1 + 0 * 3] * x[0] + a[0 + 1 * 3] * x[1];

        for (i = 1; i < n - 1; i++)
        {
            b[i] = a[2 + (i - 1) * 3] * x[i - 1] + a[1 + i * 3] * x[i] + a[0 + (i + 1) * 3] * x[i + 1];
        }

        b[n - 1] = a[2 + (n - 2) * 3] * x[n - 2] + a[1 + (n - 1) * 3] * x[n - 1] + a[0 + 0 * 3] * x[0];

        return b;
    }

    public static void r83p_print(int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_PRINT prints an R83P matrix.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored as 
        //    a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
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
        //    N must be positive.
        //
        //    Input, double A[3*N], the R83P matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r83p_print_some(n, a, 1, 1, n, n, title);

    }

    public static void r83p_print_some(int n, double[] a, int ilo, int jlo, int ihi, int jhi,
            string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_PRINT_SOME prints some of an R83P matrix.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored as 
        //    a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
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
        //    N must be positive.
        //
        //    Input, double A[3*N], the R83P matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column, to be printed.
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

            int inc = j2hi + 1 - j2lo;

            Console.WriteLine("");
            string cout = "  Col: ";
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

            if (1 < i2lo || j2hi < n)
            {
                i2lo = Math.Max(i2lo, j2lo - 1);
            }

            int i2hi = Math.Min(ihi, n);

            if (i2hi < n || 1 < j2lo)
            {
                i2hi = Math.Min(i2hi, j2hi + 1);
            }

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";

                int j2;
                for (j2 = 1; j2 <= inc; j2++)
                {
                    j = j2lo - 1 + j2;

                    if (i == n && j == 1)
                    {
                        cout += a[0 + (j - 1) * 3].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                    }
                    else
                    {
                        switch (i)
                        {
                            case 1 when j == n:
                                cout += a[2 + (j - 1) * 3].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                                break;
                            default:
                            {
                                if (1 < i - j || 1 < j - i)
                                {
                                    cout += "              ";
                                }
                                else if (j == i + 1)
                                {
                                    cout += a[0 + (j - 1) * 3].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                                }
                                else if (j == i)
                                {
                                    cout += a[1 + (j - 1) * 3].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                                }
                                else if (j == i - 1)
                                {
                                    cout += a[2 + (j - 1) * 3].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                                }

                                break;
                            }
                        }
                    }
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r83p_random(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_RANDOM randomizes an R83P matrix.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored as 
        //    a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R83P_RANDOM[3*N], the R83P matrix.
        //
    {
        double[] a = UniformRNG.r8mat_uniform_01_new(3, n, ref seed);

        return a;
    }

    public static double[] r83p_sl(int n, double[] a_lu, double[] b, int job, ref double[] work2,
            ref double[] work3, ref double work4)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_SL solves an R83P system factored by R83P_FA.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored as 
        //    a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Input, double A_LU[3*N], the LU factors from R83P_FA.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Input, int JOB, specifies the system to solve.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
        //    Input, double WORK2(N-1), WORK3(N-1), WORK4, factor data from R83P_FA.
        //
        //    Output, double R83P_SL[N], the solution to the linear system.
        //
    {
        int i;
        double[] xnm1;

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        switch (job)
        {
            case 0:
            {
                //
                //  Solve A1 * X1 = B1.
                //
                xnm1 = r83_np_sl(n - 1, a_lu, x, job);
                //
                //  X2 = B2 - A3 * X1
                //
                for (i = 0; i < n - 1; i++)
                {
                    x[i] = xnm1[i];
                }

                x[n - 1] = x[n - 1] - a_lu[0 + 0 * 3] * x[0] - a_lu[2 + (n - 2) * 3] * x[n - 2];
                //
                //  Solve A4 * X2 = X2
                //
                x[n - 1] /= work4;
                //
                //  X1 := X1 - inverse ( A1 ) * A2 * X2.
                //
                for (i = 0; i < n - 1; i++)
                {
                    x[i] -= work2[i] * x[n - 1];
                }

                break;
            }
            default:
            {
                //
                //  Solve A1' * X1 = B1.
                //
                xnm1 = r83_np_sl(n - 1, a_lu, x, job);
                //
                //  X2 := X2 - A2' * B1
                //
                for (i = 0; i < n - 1; i++)
                {
                    x[i] = xnm1[i];
                }

                x[n - 1] = x[n - 1] - a_lu[2 + (n - 1) * 3] * x[0] - a_lu[0 + (n - 1) * 3] * x[n - 2];
                //
                //  Solve A4 * X2 = X2.
                //
                x[n - 1] /= work4;
                //
                //  X1 := X1 - transpose ( inverse ( A1 ) * A3 ) * X2.
                //
                for (i = 0; i < n - 1; i++)
                {
                    x[i] -= work3[i] * x[n - 1];
                }

                break;
            }
        }

        return x;
    }

    public static double[] r83p_to_r8ge(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_TO_R8GE copies an R83P matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored as 
        //    a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Input, double A[3*N], the R83P matrix.
        //
        //    Output, double R83P_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int i;

        double[] b = new double[n * n];

        for (i = 1; i <= n; i++)
        {
            int j;
            for (j = 1; j <= n; j++)
            {
                if (i == j)
                {
                    b[i - 1 + (j - 1) * n] = a[1 + (j - 1) * 3];
                }
                else if (j == i - 1)
                {
                    b[i - 1 + (j - 1) * n] = a[2 + (j - 1) * 3];
                }
                else if (j == i + 1)
                {
                    b[i - 1 + (j - 1) * n] = a[0 + (j - 1) * 3];
                }
                else
                {
                    switch (i)
                    {
                        case 1 when j == n:
                            b[i - 1 + (j - 1) * n] = a[2 + (j - 1) * 3];
                            break;
                        default:
                        {
                            if (i == n && j == 1)
                            {
                                b[i - 1 + (j - 1) * n] = a[0 + (j - 1) * 3];
                            }
                            else
                            {
                                b[i - 1 + (j - 1) * n] = 0.0;
                            }

                            break;
                        }
                    }
                }
            }
        }

        return b;
    }

    public static double[] r83p_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83P_ZEROS zeros an R83P matrix.
        //
        //  Discussion:
        //
        //    The R83P storage format stores a periodic tridiagonal matrix is stored as 
        //    a 3 by N array, in which each row corresponds to a diagonal, and 
        //    column locations are preserved.  The matrix value 
        //    A(1,N) is stored as the array entry A(1,1), and the matrix value
        //    A(N,1) is stored as the array entry A(3,N).
        //
        //  Example:
        //
        //    Here is how an R83P matrix of order 5 would be stored:
        //
        //      A51 A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54 A15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Output, double S3P[3*N], the R83P matrix.
        //
    {
        int j;

        double[] a = new double[3 * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                a[i + j * 3] = 0.0;
            }
        }

        return a;
    }
}