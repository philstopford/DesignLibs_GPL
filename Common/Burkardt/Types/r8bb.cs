using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8bb_add(int n1, int n2, int ml, int mu, ref double[] a, int i, int j,
            double value)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_ADD adds a value to an entry in an R8BB matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, 
        //    and N2 by N2, respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N1-1.
        //
        //    Input/output, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
        //
        //    Input, int I, J, the row and column of the entry to be incremented.
        //    Some combinations of I and J are illegal.
        //
        //    Input, double VALUE, the value to be added to the (I,J)-th entry.
        //
    {
        int ij = 0;

        switch (value)
        {
            case 0.0:
                return;
        }

        if (i < 0 || n1 + n2 <= i)
        {
            Console.WriteLine("");
            Console.WriteLine("R8BB_ADD - Fatal error!");
            Console.WriteLine("  Illegal input value of row index I = " + i + "");
            return;
        }

        if (j < 0 || n1 + n2 <= j)
        {
            Console.WriteLine("");
            Console.WriteLine("R8BB_ADD - Fatal error!");
            Console.WriteLine("  Illegal input value of column index J = " + j + "");
            return;
        }

        //
        //  The A1 block of the matrix.
        //
        //  Check for out of band problems.
        //
        //  Normally, we would check the condition MU < (J-I), but the storage
        //  format requires extra entries be set aside in case of pivoting, which
        //  means that the condition becomes MU+ML < (J-I).
        //
        if (i < n1 && j < n1)
        {
            if (mu + ml < j - i || ml < i - j)
            {
                Console.WriteLine("");
                Console.WriteLine("R8BB_ADD - Warning!");
                Console.WriteLine("  Unable to add to entry (" + i + ", " + j + ").");
            }
            else
            {
                ij = i - j + ml + mu + j * (2 * ml + mu + 1);
            }
        }
        //
        //  The A2 block of the matrix.
        //
        else if (i < n1 && n1 <= j)
        {
            ij = (2 * ml + mu + 1) * n1 + (j - n1) * n1 + i;
        }
        //
        //  The A3 and A4 blocks of the matrix.
        //
        else if (n1 <= i)
        {
            ij = (2 * ml + mu + 1) * n1 + n2 * n1 + j * n2 + (i - n1);
        }

        a[ij] += value;

    }

    public static double[] r8bb_dif2(int n1, int n2, int ml, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_DIF2 sets up an R8BB second difference matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense 
        //    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    1 <= ML, 1 <= MU.
        //
        //    Output, double R8BB_DIF2[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the matrix.
        //
    {
        double[] a;
        int i;
        int j;
        double value = 0;

        a = r8vec_zeros_new((2 * ml + mu + 1) * n1 + 2 * n1 * n2 + n2 * n2);

        if (ml < 1 || mu < 1)
        {
            Console.WriteLine("");
            Console.WriteLine("R8BB_DIF2 - Fatal error!");
            Console.WriteLine("  1 <= ML and 1 <= MU required.");
            return null;
        }

        for (i = 1; i < n1 + n2; i++)
        {
            j = i - 1;
            value = -1.0;
            r8bb_set(n1, n2, ml, mu, ref a, i, j, value);
        }

        for (i = 0; i < n1 + n2; i++)
        {
            j = i;
            value = 2.0;
            r8bb_set(n1, n2, ml, mu, ref a, i, j, value);
        }

        for (i = 0; i < n1 + n2 - 1; i++)
        {
            j = i + 1;
            value = -1.0;
            r8bb_set(n1, n2, ml, mu, ref a, i, j, value);
        }

        return a;
    }

    public static int r8bb_fa(int n1, int n2, int ml, int mu, ref double[] a, ref int[] pivot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_FA factors an R8BB matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative and no greater than N1-1.
        //
        //    Input/output, double A[(2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2 ].
        //    On input, the border-banded matrix to be factored.
        //    On output, information describing a partial factorization
        //    of the original coefficient matrix.  This information is required
        //    by R8BB_SL in order to solve linear systems associated with that
        //    matrix.
        //
        //    Output, int PIVOT[N1+N2], contains pivoting information.
        //
        //    Output, int R8BB_FA, singularity flag.
        //    0, no singularity detected.
        //    nonzero, the factorization failed on the INFO-th step.
        //
    {
        double[] b;
        int i;
        int ij;
        int ik;
        int info;
        int j;
        int jk;
        int job;
        int k;
        int nband;
        double[] x;

        nband = (2 * ml + mu + 1) * n1;
        switch (n1)
        {
            //
            //  Factor the A1 band matrix, overwriting A1 by its factors.
            //
            case > 0:
            {
                info = r8gb_fa(n1, ml, mu, ref a, ref pivot);

                if (info != 0)
                {
                    return info;
                }

                break;
            }
        }

        switch (n1)
        {
            case > 0 when 0 < n2:
            {
                //
                //  Solve A1 * x = -A2 for x, and overwrite A2 by the results.
                //
                for (i = nband + 1; i <= nband + n1 * n2; i++)
                {
                    a[i - 1] = -a[i - 1];
                }

                b = r8vec_zeros_new(n1);
                x = r8vec_zeros_new(n1);

                job = 0;
                for (j = 1; j <= n2; j++)
                {
                    for (i = 0; i < n1; i++)
                    {
                        b[i] = a[nband + (j - 1) * n1 + i];
                    }

                    x = r8gb_sl(n1, ml, mu, a, pivot, b, job);
                    for (i = 0; i < n1; i++)
                    {
                        a[nband + (j - 1) * n1 + i] = x[i];
                    }
                }

                //
                //  A4 := A4 + A3 * A2.
                //
                for (i = 1; i <= n2; i++)
                {
                    for (j = 1; j <= n1; j++)
                    {
                        ij = nband + n1 * n2 + (j - 1) * n2 + i;
                        for (k = 1; k <= n2; k++)
                        {
                            ik = nband + 2 * n1 * n2 + (k - 1) * n2 + i;
                            jk = nband + (k - 1) * n1 + j;
                            a[ik - 1] += a[ij - 1] * a[jk - 1];
                        }
                    }
                }

                break;
            }
        }

        switch (n2)
        {
            //
            //  Factor A4.
            //
            case > 0:
            {
                info = r8ge_fa(n2, ref a, ref pivot, aIndex: + (nband + 2 * n1 * n2), pivotIndex:n1);

                if (info != 0)
                {
                    return info;
                }

                break;
            }
        }

        return 0;
    }

    public static double r8bb_get(int n1, int n2, int ml, int mu, ref double[] a, int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_GET gets a value of an R8BB matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, 
        //    and N2 by N2, respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N1-1.
        //
        //    Input/output, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
        //
        //    Input, int I, J, the row and column of the entry to be incremented.
        //    Some combinations of I and J are illegal.
        //
        //    Output, double R8BB_GET, the value of the (I,J)-th entry.
        //
    {
        int ij = 0;

        if (i < 0 || n1 + n2 <= i)
        {
            Console.WriteLine("");
            Console.WriteLine("R8BB_GET - Fatal error!");
            Console.WriteLine("  Illegal input value of row index I = " + i + "");
            return 1;
        }

        if (j < 0 || n1 + n2 <= j)
        {
            Console.WriteLine("");
            Console.WriteLine("R8BB_GET - Fatal error!");
            Console.WriteLine("  Illegal input value of column index J = " + j + "");
            return 1;
        }

        //
        //  The A1 block of the matrix.
        //
        //  Check for out of band problems.
        //
        //  Normally, we would check the condition MU < (J-I), but the storage
        //  format requires extra entries be set aside in case of pivoting, which
        //  means that the condition becomes MU+ML < (J-I).
        //
        if (i < n1 && j < n1)
        {
            if (mu + ml < j - i || ml < i - j)
            {
                return 0.0;
            }

            ij = i - j + ml + mu + j * (2 * ml + mu + 1);
        }
        //
        //  The A2 block of the matrix.
        //
        else if (i < n1 && n1 <= j)
        {
            ij = (2 * ml + mu + 1) * n1 + (j - n1) * n1 + i;
        }
        //
        //  The A3 and A4 blocks of the matrix.
        //
        else if (n1 <= i)
        {
            ij = (2 * ml + mu + 1) * n1 + n2 * n1 + j * n2 + (i - n1);
        }

        return a[ij];
    }

    public static double[] r8bb_indicator(int n1, int n2, int ml, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_INDICATOR sets up an R8BB indicator matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //    The matrix is actually stored as a vector, and we will simply suggest
        //    the structure and values of the indicator matrix as:
        //
        //      00 00 00 00 00
        //      00 00 13 24 35     16 17     61 62 63 64 65     66 67
        //      00 12 23 34 45  +  26 27  +  71 72 73 74 75  +  76 77
        //      11 22 33 44 55     36 37     
        //      21 32 43 54 00     46 47     
        //                         56 57     
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative and no greater than N1-1.
        //
        //    Output, double R8BB_INDICATOR[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], 
        //    the matrix.
        //
    {
        double[] a;
        int base_;
        int fac;
        int i;
        int j;
        int row;

        a = r8vec_zeros_new((2 * ml + mu + 1) * n1 + 2 * n1 * n2 + n2 * n2);

        fac = (int) Math.Pow(10, (int) Math.Log10(n1 + n2) + 1);
        //
        //  Set the banded matrix A1.
        //
        for (j = 1; j <= n1; j++)
        {
            for (row = 1; row <= 2 * ml + mu + 1; row++)
            {
                i = row + j - ml - mu - 1;
                if (ml < row && 1 <= i && i <= n1)
                {
                    a[row - 1 + (j - 1) * (2 * ml + mu + 1)] = fac * i + j;
                }
            }
        }

        //
        //  Set the N1 by N2 rectangular strip A2.
        //
        base_ = (2 * ml + mu + 1) * n1;

        for (i = 1; i <= n1; i++)
        {
            for (j = n1 + 1; j <= n1 + n2; j++)
            {
                a[base_ + i - 1 + (j - n1 - 1) * n1] = fac * i + j;
            }
        }

        //
        //  Set the N2 by N1 rectangular strip A3.
        //
        base_ = (2 * ml + mu + 1) * n1 + n1 * n2;

        for (i = n1 + 1; i <= n1 + n2; i++)
        {
            for (j = 1; j <= n1; j++)
            {
                a[base_ + i - n1 - 1 + (j - 1) * n2] = fac * i + j;
            }
        }

        //
        //  Set the N2 by N2 square A4.
        //
        base_ = (2 * ml + mu + 1) * n1 + n1 * n2 + n2 * n1;

        for (i = n1 + 1; i <= n1 + n2; i++)
        {
            for (j = n1 + 1; j <= n1 + n2; j++)
            {
                a[base_ + i - n1 - 1 + (j - n1 - 1) * n2] = fac * i + j;
            }
        }

        return a;
    }

    public static double[] r8bb_mtv(int n1, int n2, int ml, int mu, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_MTV multiplies a vector by an R8BB matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, 
        //    and N2 by N2, respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative and no greater than N1-1.
        //
        //    Input, double A[(2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8BB matrix.
        //
        //    Input, double X[N1+N2], the vector to multiply A.
        //
        //    Output, double R8BB_MTV[N1+N2], the product X times A.
        //
    {
        double[] b;
        int i;
        int ihi;
        int ij;
        int ilo;
        int j;
        //
        //  Initialize B.
        //
        b = r8vec_zeros_new(n1 + n2);
        //
        //  Multiply by A1.
        //
        for (j = 1; j <= n1; j++)
        {
            ilo = Math.Max(1, j - mu - ml);
            ihi = Math.Min(n1, j + ml);
            ij = (j - 1) * (2 * ml + mu + 1) - j + ml + mu + 1;
            for (i = ilo; i <= ihi; i++)
            {
                b[j - 1] += x[i - 1] * a[ij + i - 1];
            }
        }

        //
        //  Multiply by A2.
        //
        for (j = n1 + 1; j <= n1 + n2; j++)
        {
            ij = (2 * ml + mu + 1) * n1 + (j - n1 - 1) * n1;
            for (i = 1; i <= n1; i++)
            {
                b[j - 1] += x[i - 1] * a[ij + i - 1];
            }
        }

        //
        //  Multiply by A3 and A4.
        //
        for (j = 1; j <= n1 + n2; j++)
        {
            ij = (2 * ml + mu + 1) * n1 + n1 * n2 + (j - 1) * n2 - n1;
            for (i = n1 + 1; i <= n1 + n2; i++)
            {
                b[j - 1] += x[i - 1] * a[ij + i - 1];
            }
        }

        return b;
    }

    public static double[] r8bb_mv(int n1, int n2, int ml, int mu, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_MV multiplies an R8BB matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative and no greater than N1-1.
        //
        //    Input, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
        //
        //    Input, double X[N1+N2], the vector to be multiplied by A.
        //
        //    Output, double R8BB_MV[N1+N2], the result of multiplying A by X.
        //
    {
        double[] b;
        int i;
        int ihi;
        int ij;
        int ilo;
        int j;
        //
        //  Initialize B.
        //
        b = r8vec_zeros_new(n1 + n2);
        //
        //  Multiply by A1.
        //
        for (j = 1; j <= n1; j++)
        {
            ilo = Math.Max(1, j - mu - ml);
            ihi = Math.Min(n1, j + ml);
            ij = (j - 1) * (2 * ml + mu + 1) - j + ml + mu + 1;
            for (i = ilo; i <= ihi; i++)
            {
                b[i - 1] += a[ij + i - 1] * x[j - 1];
            }
        }

        //
        //  Multiply by A2.
        //
        for (j = n1 + 1; j <= n1 + n2; j++)
        {
            ij = (2 * ml + mu + 1) * n1 + (j - n1 - 1) * n1;
            for (i = 1; i <= n1; i++)
            {
                b[i - 1] += a[ij + i - 1] * x[j - 1];
            }
        }

        //
        //  Multiply by A3 and A4.
        //
        for (j = 1; j <= n1 + n2; j++)
        {
            ij = (2 * ml + mu + 1) * n1 + n1 * n2 + (j - 1) * n2 - n1;
            for (i = n1 + 1; i <= n1 + n2; i++)
            {
                b[i - 1] += a[ij + i - 1] * x[j - 1];
            }
        }

        return b;
    }

    public static void r8bb_print(int n1, int n2, int ml, int mu, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_PRINT prints an R8BB matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, 
        //    and N2 by N2, respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
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
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N1-1.
        //
        //    Input, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8bb_print_some(n1, n2, ml, mu, a, 0, 0, n1 + n2 - 1, n1 + n2 - 1, title);
    }

    public static void r8bb_print_some(int n1, int n2, int ml, int mu, double[] a, int ilo,
            int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_PRINT_SOME prints some of an R8BB matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N1-1.
        //
        //    Input, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
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
        int ij;
        int j;
        int j2hi;
        int j2lo;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n1 + n2 - 1);
            j2hi = Math.Min(j2hi, jhi);

            Console.WriteLine("");
            cout = "  Col: ";

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
            i2hi = Math.Min(ihi, n1 + n2 - 1);

            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString().PadLeft(4) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    aij = 0.0;

                    if (i < n1 && j < n1)
                    {
                        if (j - i <= mu + ml && i - j <= ml)
                        {
                            ij = i - j + ml + mu + j * (2 * ml + mu + 1);
                            aij = a[ij];
                        }
                    }
                    else if (i < n1 && n1 <= j)
                    {
                        ij = (2 * ml + mu + 1) * n1 + (j - n1) * n1 + i;
                        aij = a[ij];
                    }
                    else if (n1 <= i)
                    {
                        ij = (2 * ml + mu + 1) * n1 + n2 * n1 + j * n2 + (i - n1);
                        aij = a[ij];
                    }

                    cout += aij.ToString().PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8bb_random(int n1, int n2, int ml, int mu, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_RANDOM randomizes an R8BB matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative and no greater than N1-1.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8BB_RANDOM[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the matrix.
        //
    {
        double[] a;
        int i;
        int j;
        int row;

        a = r8vec_zeros_new((2 * ml + mu + 1) * n1 + 2 * n1 * n2 + n2 * n2);
        //
        //  Randomize the banded matrix A1.
        //
        for (j = 1; j <= n1; j++)
        {
            for (row = 1; row <= 2 * ml + mu + 1; row++)
            {
                i = row + j - ml - mu - 1;
                if (ml < row && 1 <= i && i <= n1)
                {
                    a[row - 1 + (j - 1) * (2 * ml + mu + 1)] = UniformRNG.r8_uniform_01(ref seed);
                }
            }
        }

        //
        //  Randomize the rectangular strips A2+A3+A4.
        //
        for (i = (2 * ml + mu + 1) * n1; i < (2 * ml + mu + 1) * n1 + 2 * n1 * n2 + n2 * n2; i++)
        {
            a[i] = UniformRNG.r8_uniform_01(ref seed);
        }

        return a;
    }

    public static void r8bb_set(int n1, int n2, int ml, int mu, ref double[] a, int i, int j,
            double value)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_SET sets a value of an R8BB matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N1-1.
        //
        //    Input/output, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
        //
        //    Input, int I, J, the row and column of the entry to be incremented.
        //    Some combinations of I and J are illegal.
        //
        //    Input, double VALUE, the value to be assigned to the (I,J)-th entry.
        //
    {
        int ij = 0;

        if (i < 0 || n1 + n2 <= i)
        {
            Console.WriteLine("");
            Console.WriteLine("R8BB_SET - Fatal error!");
            Console.WriteLine("  Illegal input value of row index I = " + i + "");
            return;
        }

        if (j < 0 || n1 + n2 <= j)
        {
            Console.WriteLine("");
            Console.WriteLine("R8BB_SET - Fatal error!");
            Console.WriteLine("  Illegal input value of column index J = " + j + "");
            return;
        }

        //
        //  The A1 block of the matrix.
        //
        //  Check for out of band problems.
        //
        //  Normally, we would check the condition MU < (J-I), but the storage
        //  format requires extra entries be set aside in case of pivoting, which
        //  means that the condition becomes MU+ML < (J-I).
        //
        if (i < n1 && j < n1)
        {
            if (mu + ml < j - i || ml < i - j)
            {
                Console.WriteLine("");
                Console.WriteLine("R8BB_SET - Warning!");
                Console.WriteLine("  Unable to set entry (" + i + ", " + j + ").");
            }
            else
            {
                ij = i - j + ml + mu + j * (2 * ml + mu + 1);
            }
        }
        //
        //  The A2 block of the matrix.
        //
        else if (i < n1 && n1 <= j)
        {
            ij = (2 * ml + mu + 1) * n1 + (j - n1) * n1 + i;
        }
        //
        //  The A3 and A4 blocks of the matrix.
        //
        else if (n1 <= i)
        {
            ij = (2 * ml + mu + 1) * n1 + n2 * n1 + j * n2 + (i - n1);
        }

        a[ij] = value;

    }

    public static double[] r8bb_sl(int n1, int n2, int ml, int mu, double[] a_lu, int[] pivot,
            double[] b)

        //****************************************************************************80
        //
        //  Discussion:
        //
        //    R8BB_SL solves an R8BB system factored by SBB_FA.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 November 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative and no greater than N1-1.
        //
        //    Input, double A_LU[(2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the LU 
        //    factors from R8BB_FA.
        //
        //    Input, int PIVOT[N1+N2], the pivoting information from R8BB_FA.
        //
        //    Input, double B[N1+N2], the right hand side.
        //
        //    Output, double R8BB_SL[N1+N2], the solution.
        //
    {
        double[] b22 = null;
        int i;
        int ij;
        int j;
        int job;
        int nband;
        double[] x;
        double[] x1 = null;
        double[] x2 = null;

        nband = (2 * ml + mu + 1) * n1;
        switch (n1)
        {
            //
            //  Set X1 := inverse(A1) * B1.
            //
            case > 0:
                job = 0;
                x1 = r8gb_sl(n1, ml, mu, a_lu, pivot, b, job);
                break;
        }

        switch (n2)
        {
            //
            //  Modify the right hand side of the second linear subsystem.
            //  Set B22 := B2 - A3*X1.
            //
            case > 0:
            {
                b22 = r8vec_zeros_new(n2);

                for (i = 0; i < n2; i++)
                {
                    b22[i] = b[n1 + i];
                    for (j = 0; j < n1; j++)
                    {
                        ij = nband + n1 * n2 + j * n2 + i;
                        b22[i] -= a_lu[ij] * x1[j];
                    }
                }

                break;
            }
        }

        switch (n2)
        {
            //
            //  Set X2 := inverse(A4) * B22.
            //
            case > 0:
                job = 0;
                x2 = r8ge_sl_new(n2, a_lu, pivot, b22, job, aluIndex: +(nband + 2 * n1 * n2), pivotIndex: +n1);
                break;
        }

        //
        //  Modify the first subsolution.
        //  Set X1 := X1 + A2*X2.
        //
        for (i = 0; i < n1; i++)
        {
            for (j = 0; j < n2; j++)
            {
                ij = nband + j * n1 + i;
                x1[i] += a_lu[ij] * x2[j];
            }
        }

        //
        //  Set X = [ X1 | X2 ].
        //
        x = r8vec_zeros_new(n1 + n2);

        switch (n1)
        {
            case > 0:
            {
                for (i = 0; i < n1; i++)
                {
                    x[i] = x1[i];
                }

                break;
            }
        }

        switch (n2)
        {
            case > 0:
            {
                for (i = 0; i < n2; i++)
                {
                    x[n1 + i] = x2[i];
                }

                break;
            }
        }

        return x;
    }

    public static double[] r8bb_to_r8ge(int n1, int n2, int ml, int mu, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_TO_R8GE copies an R8BB matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N1-1.
        //
        //    Input, double A[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
        //
        //    Output, double R8BB_TO_R8GE[(N1+N2)*(N1+N2)], the R8GE matrix.
        //
    {
        double[] b;
        int i;
        int ij;
        int j;

        b = r8vec_zeros_new((n1 + n2) * (n1 + n2));

        for (i = 1; i <= n1; i++)
        {
            for (j = 1; j <= n1; j++)
            {
                if (mu + ml < j - i || ml < i - j)
                {
                    b[i - 1 + (j - 1) * (n1 + n2)] = 0.0;
                }
                else
                {
                    ij = i - j + ml + mu + 1 + (j - 1) * (2 * ml + mu + 1);
                    b[i - 1 + (j - 1) * (n1 + n2)] = a[ij - 1];
                }
            }
        }

        for (i = 1; i <= n1; i++)
        {
            for (j = n1 + 1; j <= n1 + n2; j++)
            {
                ij = (2 * ml + mu + 1) * n1 + (j - n1 - 1) * n1 + i;
                b[i - 1 + (j - 1) * (n1 + n2)] = a[ij - 1];
            }
        }

        for (i = n1 + 1; i <= n1 + n2; i++)
        {
            for (j = 1; j <= n1 + n2; j++)
            {
                ij = (2 * ml + mu + 1) * n1 + n2 * n1 + (j - 1) * n2 + (i - n1);
                b[i - 1 + (j - 1) * (n1 + n2)] = a[ij - 1];
            }
        }

        return b;
    }

    public static double[] r8bb_zeros(int n1, int n2, int ml, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BB_ZEROS zeros an R8BB matrix.
        //
        //  Discussion:
        //
        //    The R8BB storage format is for a border banded matrix.  Such a
        //    matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
        //    general band format.  The reason for the factor of 2 in front of
        //    ML is to allocate space that may be required if pivoting occurs.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Example:
        //
        //    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
        //
        //       00
        //       00  00
        //       00  00  00 --- ---
        //      A11 A12 A13  00 ---  A16 A17
        //      A21 A22 A23 A24  00  A26 A27
        //      --- A32 A33 A34 A35  A36 A37
        //      --- --- A43 A44 A45  A46 A47
        //      --- --- --- A54 A55  A56 A57
        //                       00
        //
        //      A61 A62 A63 A64 A65  A66 A67
        //      A71 A72 A73 A74 A75  A76 A77
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative and no greater than N1-1.
        //
        //    Output, double R8BB_ZERO[(2*ML+MU+1)*N1+2*N1*N2+N2*N2], the R8BB matrix.
        //
    {
        double[] a;

        a = r8vec_zeros_new((2 * ml + mu + 1) * n1 + 2 * n1 * n2 + n2 * n2);

        return a;
    }

}