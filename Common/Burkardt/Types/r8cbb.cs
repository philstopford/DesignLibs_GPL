using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8cbb_add(int n1, int n2, int ml, int mu, ref double[] a, int i, int j,
            double value)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_ADD adds a value to an entry of an R8CBB matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
        //    Input/output, double A[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
        //
        //    Input, int I, J, the indices of the entry to be incremented.
        //
        //    Input, double VALUE, the value to be added to the (I,J) entry.
        //
    {
        int ij = 0;

        switch (value)
        {
            case 0.0:
                return;
        }

        //
        //  Check for I or J out of bounds.
        //
        if (i < 0 || n1 + n2 <= i)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CBB_ADD - Fatal error!");
            Console.WriteLine("  Illegal input value of row index I = " + i + "");
            return;
        }

        if (j < 0 || n1 + n2 <= j)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CBB_ADD - Fatal error!");
            Console.WriteLine("  Illegal input value of column index J = " + j + "");
            return;
        }

        //
        //  The A1 block of the matrix.
        //
        //  Check for out of band problems.
        //
        if (i < n1 && j < n1)
        {
            if (mu < j - i || ml < i - j)
            {
                Console.WriteLine("");
                Console.WriteLine("R8CBB_ADD - Warning!");
                Console.WriteLine("  Unable to add to entry (" + i + ", " + j + ").");
                return;
            }

            ij = i - j + mu + j * (ml + mu + 1);
        }
        //
        //  The A2 block of the matrix:
        //
        else if (i < n1 && n1 <= j)
        {
            ij = (ml + mu + 1) * n1 + (j - n1) * n1 + i;
        }
        //
        //  The A3 and A4 blocks of the matrix.
        //
        else if (n1 <= i)
        {
            ij = (ml + mu + 1) * n1 + n2 * n1 + j * n2 + (i - n1);
        }

        a[ij] += value;
    }

    public static double[] r8cbb_dif2(int n1, int n2, int ml, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_DIF2 sets up an R8CBB second difference matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 July 2016
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
        //    Output, double R8CBB_DIF2[(ML+MU+1)*N1+2*N1*N2+N2*N2], the matrix.
        //
    {
        double[] a;
        int i;
        int j;
        double value = 0;

        a = r8vec_zeros_new((ml + mu + 1) * n1 + 2 * n1 * n2 + n2 * n2);

        for (i = 1; i < n1 + n2; i++)
        {
            j = i - 1;
            value = -1.0;
            r8cbb_set(n1, n2, ml, mu, ref a, i, j, value);
        }

        for (i = 0; i < n1 + n2; i++)
        {
            j = i;
            value = 2.0;
            r8cbb_set(n1, n2, ml, mu, ref a, i, j, value);
        }

        for (i = 0; i < n1 + n2 - 1; i++)
        {
            j = i + 1;
            value = -1.0;
            r8cbb_set(n1, n2, ml, mu, ref a, i, j, value);
        }

        return a;
    }

    public static int r8cbb_fa(int n1, int n2, int ml, int mu, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_FA factors an R8CBB matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //
        //    Once the matrix has been factored by SCCB_FA, SCCB_SL may be called
        //    to solve linear systems involving the matrix.
        //
        //    SCCB_FA uses special non-pivoting versions of LINPACK routines to
        //    carry out the factorization.  The special version of the banded
        //    LINPACK solver also results in a space saving, since no entries
        //    need be set aside for fill in due to pivoting.
        //
        //    The linear system must be border banded, of the form:
        //
        //      ( A1 A2 ) (X1) = (B1)
        //      ( A3 A4 ) (X2)   (B2)
        //
        //    where A1 is a (usually big) banded square matrix, A2 and A3 are
        //    column and row strips which may be nonzero, and A4 is a dense
        //    square matrix.
        //
        //    The algorithm rewrites the system as:
        //
        //         X1 + inverse(A1) A2 X2 = inverse(A1) B1
        //
        //      A3 X1 +             A4 X2 = B2
        //
        //    and then rewrites the second equation as
        //
        //      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
        //
        //    The algorithm will certainly fail if the matrix A1 is singular,
        //    or requires pivoting.  The algorithm will also fail if the A4 matrix,
        //    as modified during the process, is singular, or requires pivoting.
        //    All these possibilities are in addition to the failure that will
        //    if the total matrix A is singular.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2004
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
        //    Input/output, double A[ (ML+MU+1)*N1 + 2*N1*N2 + N2*N2].
        //    On input, A contains the compact border-banded coefficient matrix.
        //    On output, A contains information describing a partial factorization
        //    of the original coefficient matrix.  
        //
        //    Output, int R8CBB_FA, singularity flag.
        //    0, no singularity detected.
        //    nonzero, the factorization failed on the INFO-th step.
        //
    {
        double[] b1;
        int i;
        int ij;
        int ik;
        int info;
        int j;
        int jk;
        int job;
        int k;
        int nband;
        double[] x1;

        nband = (ml + mu + 1) * n1;
        switch (n1)
        {
            //
            //  Factor the A1 band matrix, overwriting A1 by its factors.
            //
            case > 0:
            {
                info = r8cb_np_fa(n1, ml, mu, ref a);
                if (info != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8CBB_FA - Fatal error!");
                    Console.WriteLine("  R8CB_NP_FA returned INFO = " + info + "");
                    Console.WriteLine("  Factoring failed for column INFO.");
                    Console.WriteLine("  The band matrix A1 is singular.");
                    Console.WriteLine("  This algorithm cannot continue!");
                    return 1;
                }

                break;
            }
        }

        switch (n1)
        {
            case > 0 when 0 < n2:
            {
                //
                //  Set A2 := -inverse(A1) * A2.
                //
                for (j = 0; j < n2; j++)
                {
                    for (i = 0; i < n1; i++)
                    {
                        a[nband + i + j * n1] = -a[nband + i + j * n1];
                    }
                }

                b1 = r8vec_zeros_new(n1);
                x1 = r8vec_zeros_new(n1);
                job = 0;

                for (j = 0; j < n2; j++)
                {
                    for (i = 0; i < n1; i++)
                    {
                        b1[i] = a[nband + i + j * n1];
                    }

                    x1 = r8cb_np_sl(n1, ml, mu, a, b1, job);
                    for (i = 0; i < n1; i++)
                    {
                        a[nband + i + j * n1] = x1[i];
                    }
                }

                //
                //  Set A4 := A4 + A3*A2
                //
                for (i = 1; i <= n2; i++)
                {
                    for (j = 1; j <= n1; j++)
                    {
                        ij = nband + n1 * n2 + (j - 1) * n2 + i - 1;
                        for (k = 1; k <= n2; k++)
                        {
                            ik = nband + 2 * n1 * n2 + (k - 1) * n2 + i - 1;
                            jk = nband + (k - 1) * n1 + j - 1;
                            a[ik] += a[ij] * a[jk];
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
                info = r8ge_np_fa(n2, ref a,  + (nband + 2 * n1 * n2));

                if (info != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8CBB_FA - Fatal error!");
                    Console.WriteLine("  R8GE_NP_FA returned INFO = " + info + "");
                    Console.WriteLine("  This indicates singularity in column " + n1 + info + ".");
                    Console.WriteLine("  The dense matrix A4 is singular.");
                    Console.WriteLine("  This algorithm cannot continue!");
                    return 1;
                }

                break;
            }
        }

        return 0;
    }

    public static double r8cbb_get(int n1, int n2, int ml, int mu, ref double[] a, int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_GET gets the value of an entry of an R8CBB matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
        //    Input/output, double A[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
        //
        //    Input, int I, J, the indices of the entry to be incremented.
        //
        //    Output, double R8CBB_GET, the value of the (I,J) entry.
        //
    {
        int ij = 0;
        //
        //  Check for I or J out of bounds.
        //
        if (i < 0 || n1 + n2 <= i)
        {
            return 0.0;
        }

        if (j < 0 || n1 + n2 <= j)
        {
            return 0.0;
        }

        //
        //  The A1 block of the matrix.
        //
        //  Check for out of band problems.
        //
        if (i < n1 && j < n1)
        {
            if (mu < j - i || ml < i - j)
            {
                return 0.0;
            }

            ij = i - j + mu + j * (ml + mu + 1);
        }
        //
        //  The A2 block of the matrix:
        //
        else if (i < n1 && n1 <= j)
        {
            ij = (ml + mu + 1) * n1 + (j - n1) * n1 + i;
        }
        //
        //  The A3 and A4 blocks of the matrix.
        //
        else if (n1 <= i)
        {
            ij = (ml + mu + 1) * n1 + n2 * n1 + j * n2 + (i - n1);
        }

        return a[ij];
    }

    public static double[] r8cbb_indicator(int n1, int n2, int ml, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_INDICATOR sets up an R8CBB indicator matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 January 2004
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
        //    Output, double R8CBB_INDICATOR[(ML+MU+1)*N1+2*N1*N2+N2*N2], the R8CBB
        //    indicator matrix.
        //
    {
        double[] a;
        int base_;
        int fac;
        int i;
        int j;
        int row;

        a = r8vec_zeros_new((ml + mu + 1) * n1 + 2 * n1 * n2 + n2 * n2);

        fac = (int) Math.Pow(10, (int) Math.Log10(n1 + n2) + 1);
        //
        //  Set the banded matrix A1.
        //
        for (j = 1; j <= n1; j++)
        {
            for (row = 1; row <= ml + mu + 1; row++)
            {
                i = row + j - mu - 1;
                a[row - 1 + (j - 1) * (ml + mu + 1)] = i switch
                {
                    >= 1 when i <= n1 => fac * i + j,
                    _ => a[row - 1 + (j - 1) * (ml + mu + 1)]
                };
            }
        }

        //
        //  Set the N1 by N2 rectangular strip A2.
        //
        base_ = (ml + mu + 1) * n1;

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
        base_ = (ml + mu + 1) * n1 + n1 * n2;

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
        base_ = (ml + mu + 1) * n1 + n1 * n2 + n2 * n1;

        for (i = n1 + 1; i <= n1 + n2; i++)
        {
            for (j = n1 + 1; j <= n1 + n2; j++)
            {
                a[base_ + i - n1 - 1 + (j - n1 - 1) * n2] = fac * i + j;
            }
        }

        return a;
    }

    public static double[] r8cbb_mtv(int n1, int n2, int ml, int mu, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_MTV multiplies a vector by an R8CBB matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
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
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N1-1.
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, double A[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
        //
        //    Input, double X[N1+N2], the vector to multiply the matrix.
        //
        //    Output, double R8CBB_MTV[N1+N2], the product X * A.
        //
    {
        double[] b;
        int i;
        int ihi;
        int ij;
        int ilo;
        int j;
        //
        //  Set B to zero.
        //
        b = r8vec_zeros_new(n1 + n2);
        //
        //  Multiply by A1.
        //
        for (j = 1; j <= n1; j++)
        {
            ilo = Math.Max(1, j - mu);
            ihi = Math.Min(n1, j + ml);
            ij = (j - 1) * (ml + mu + 1) - j + mu + 1;
            for (i = ilo; i <= ihi; i++)
            {
                b[j] += x[i - 1] * a[ij + i - 1];
            }
        }

        //
        //  Multiply by A2.
        //
        for (j = n1 + 1; j <= n1 + n2; j++)
        {
            ij = (ml + mu + 1) * n1 + (j - n1 - 1) * n1;
            for (i = 1; i <= n1; i++)
            {
                b[j] += x[i - 1] * a[ij + i - 1];
            }
        }

        //
        //  Multiply by A3 and A4.
        //
        for (j = 1; j <= n1 + n2; j++)
        {
            ij = (ml + mu + 1) * n1 + n1 * n2 + (j - 1) * n2 - n1;
            for (i = n1 + 1; i <= n1 + n2; i++)
            {
                b[j - 1] += x[i - 1] * a[ij + i - 1];
            }
        }

        return b;
    }

    public static double[] r8cbb_mv(int n1, int n2, int ml, int mu, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_MV multiplies an R8CBB matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 October 1998
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N1-1.
        //
        //    Input, int N1, N2, the order of the banded and dense blocks.
        //    N1 and N2 must be nonnegative, and at least one must be positive.
        //
        //    Input, double A[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
        //
        //    Input, double X[N1+N2], the vector to be multiplied by A.
        //
        //    Output, double R8CBB_MV[N1+N2], the result of multiplying A by X.
        //
    {
        double[] b;
        int i;
        int ihi;
        int ij;
        int ilo;
        int j;
        //
        //  Set B to zero.
        //
        b = r8vec_zeros_new(n1 + n2);
        //
        //  Multiply by A1.
        //
        for (j = 1; j <= n1; j++)
        {
            ilo = Math.Max(1, j - mu);
            ihi = Math.Min(n1, j + ml);
            ij = (j - 1) * (ml + mu + 1) - j + mu + 1;
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
            ij = (ml + mu + 1) * n1 + (j - n1 - 1) * n1;
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
            ij = (ml + mu + 1) * n1 + n1 * n2 + (j - 1) * n2 - n1;
            for (i = n1 + 1; i <= n1 + n2; i++)
            {
                b[i - 1] += a[ij + i - 1] * x[j - 1];
            }
        }

        return b;
    }

    public static void r8cbb_print(int n1, int n2, int ml, int mu, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_PRINT prints an R8CBB matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
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
        //    Input, double A[(ML+MU+1)*N1+2*N1*N2+N2*N2], the R8CBB matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8cbb_print_some(n1, n2, ml, mu, a, 0, 0, n1 + n2 - 1, n1 + n2 - 1, title);
    }

    public static void r8cbb_print_some(int n1, int n2, int ml, int mu, double[] a, int ilo,
            int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_PRINT_SOME prints some of an R8CBB matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
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
        //    Input, double A[(ML+MU+1)*N1+2*N1*N2+N2*N2], the R8CBB matrix.
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
                    aij = r8cbb_get(n1, n2, ml, mu, ref a, i, j);
                    cout += aij.ToString().PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8cbb_random(int n1, int n2, int ml, int mu, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_RANDOM randomizes an R8CBB matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
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
        //    Output, double R8CBB_RANDOM[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
        //
    {
        double[] a;
        int i;
        int j;
        int row;

        a = r8vec_zeros_new((ml + mu + 1) * n1 + 2 * n1 * n2 + n2 * n2);
        //
        //  Randomize the banded matrix A1.
        //
        for (j = 1; j <= n1; j++)
        {
            for (row = 1; row <= ml + mu + 1; row++)
            {
                i = row + j - mu - 1;
                a[row - 1 + (j - 1) * (ml + mu + 1)] = i switch
                {
                    >= 1 when i <= n1 => UniformRNG.r8_uniform_01(ref seed),
                    _ => a[row - 1 + (j - 1) * (ml + mu + 1)]
                };
            }
        }

        //
        //  Randomize the rectangular strips A2+A3+A4.
        //
        for (i = (ml + mu + 1) * n1 + 1; i <= (ml + mu + 1) * n1 + 2 * n1 * n2 + n2 * n2; i++)
        {
            a[i - 1] = UniformRNG.r8_uniform_01(ref seed);
        }

        return a;
    }

    public static void r8cbb_set(int n1, int n2, int ml, int mu, ref double[] a, int i, int j,
            double value)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_SET sets an entry of an R8CBB matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
        //    Input/output, double A[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the R8CBB matrix.
        //
        //    Input, int I, J, the indices of the entry to be incremented.
        //
        //    Input, double VALUE, the value to be assigned to the (I,J) entry.
        //
    {
        int ij = 0;
        //
        //  Check for I or J out of bounds.
        //
        if (i < 0 || n1 + n2 <= i)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CBB_SET - Fatal error!");
            Console.WriteLine("  Illegal input value of row index I = " + i + "");
            return;
        }

        if (j < 0 || n1 + n2 <= j)
        {
            Console.WriteLine("");
            Console.WriteLine("R8CBB_SET - Fatal error!");
            Console.WriteLine("  Illegal input value of column index J = " + j + "");
            return;
        }

        //
        //  The A1 block of the matrix.
        //
        //  Check for out of band problems.
        //
        if (i < n1 && j < n1)
        {
            if (mu < j - i || ml < i - j)
            {
                Console.WriteLine("");
                Console.WriteLine("R8CBB_SET - Warning!");
                Console.WriteLine("  Unable to set entry (" + i + ", " + j + ").");
                return;
            }

            ij = i - j + mu + j * (ml + mu + 1);
        }
        //
        //  The A2 block of the matrix:
        //
        else if (i < n1 && n1 <= j)
        {
            ij = (ml + mu + 1) * n1 + (j - n1) * n1 + i;
        }
        //
        //  The A3 and A4 blocks of the matrix.
        //
        else if (n1 <= i)
        {
            ij = (ml + mu + 1) * n1 + n2 * n1 + j * n2 + (i - n1);
        }

        a[ij] = value;
    }

    public static double[] r8cbb_sl(int n1, int n2, int ml, int mu, double[] a_lu, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_SL solves an R8CBB system factored by R8CBB_FA.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //
        //    The linear system A * x = b is decomposable into the block system:
        //
        //      ( A1 A2 ) * (X1) = (B1)
        //      ( A3 A4 )   (X2)   (B2)
        //
        //    where A1 is a (usually big) banded square matrix, A2 and A3 are
        //    column and row strips which may be nonzero, and A4 is a dense
        //    square matrix.
        //
        //    All the arguments except B are input quantities only, which are
        //    not changed by the routine.  They should have exactly the same values
        //    they had on exit from R8CBB_FA.
        //
        //    If more than one right hand side is to be solved, with the same
        //    matrix, R8CBB_SL should be called repeatedly.  However, R8CBB_FA only
        //    needs to be called once to create the factorization.
        //
        //    See the documentation of R8CBB_FA for details on the matrix storage.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 January 2004
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
        //    Input, double A_LU[ (ML+MU+1)*N1 + 2*N1*N2 + N2*N2].
        //    the LU factors from R8CBB_FA.
        //
        //    Input, double B[N1+N2], the right hand side of the linear system.
        //
        //    Output, double R8CBB_SL[N1+N2], the solution.
        //
    {
        double[] b2;
        int i;
        int ij;
        int j;
        int job;
        int nband;
        double[] x;
        double[] x1 = null;
        double[] x2 = null;

        nband = (ml + mu + 1) * n1;
        switch (n1)
        {
            //
            //  Set X1 := inverse(A1) * B1.
            //
            case > 0:
                job = 0;
                x1 = r8cb_np_sl(n1, ml, mu, a_lu, b, job);
                break;
        }

        //
        //  Modify the right hand side of the second linear subsystem.
        //  Set B2 = B2-A3*X1.
        //
        b2 = r8vec_zeros_new(n2);

        for (i = 0; i < n2; i++)
        {
            b2[i] = b[n1 + i];
        }

        for (j = 0; j < n1; j++)
        {
            for (i = 0; i < n2; i++)
            {
                ij = nband + n1 * n2 + j * n2 + i;
                b2[i] -= a_lu[ij] * x1[j];
            }
        }

        switch (n2)
        {
            //
            //  Solve A4*X2 = B2.
            //
            case > 0:
                job = 0;
                x2 = r8ge_np_sl(n2, a_lu, b2, job, aluIndex: + (nband + 2 * n1 * n2));
                break;
        }

        //
        //  Modify the first subsolution.
        //  Set X1 = X1+A2*X2.
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
        //  Collect X1 and X2 into X.
        //
        x = r8vec_zeros_new(n1 + n2);

        for (i = 0; i < n1; i++)
        {
            x[i] = x1[i];
        }

        for (i = 0; i < n2; i++)
        {
            x[n1 + i] = x2[i];
        }

        return x;
    }

    public static double[] r8cbb_to_r8ge(int n1, int n2, int ml, int mu, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_TO_R8GE copies an R8CBB matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 October 2003
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
        //    Input, double A[(ML+MU+1)*N1+2*N1*N2+N2*N2], the R8CBB matrix.
        //
        //    Output, double R8CBB_TO_R8GE[(N1+N2)*(N1+N2)], the R8GE matrix.
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
                    ij = i - j + mu + 1 + (j - 1) * (ml + mu + 1);
                    b[i - 1 + (j - 1) * (n1 + n2)] = a[ij - 1];
                }
            }
        }

        for (i = 1; i <= n1; i++)
        {
            for (j = n1 + 1; j <= n1 + n2; j++)
            {
                ij = (ml + mu + 1) * n1 + (j - n1 - 1) * n1 + i;
                b[i - 1 + (j - 1) * (n1 + n2)] = a[ij - 1];
            }
        }

        for (i = n1 + 1; i <= n1 + n2; i++)
        {
            for (j = 1; j <= n1 + n2; j++)
            {
                ij = (ml + mu + 1) * n1 + n2 * n1 + (j - 1) * n2 + (i - n1);
                b[i - 1 + (j - 1) * (n1 + n2)] = a[ij - 1];
            }
        }

        return b;
    }

    public static double[] r8cbb_zeros(int n1, int n2, int ml, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8CBB_ZEROS zeros an R8CBB matrix.
        //
        //  Discussion:
        //
        //    The R8CBB storage format is for a compressed border banded matrix.  
        //    Such a matrix has the logical form:
        //
        //      A1 | A2
        //      ---+---
        //      A3 | A4
        //
        //    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
        //    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
        //    respectively.  
        //
        //    The R8CBB format is the same as the R8BB format, except that the banded
        //    matrix A1 is stored in compressed band form rather than standard
        //    banded form.  In other words, we do not include the extra room
        //    set aside for fill in during pivoting.
        //
        //    A should be defined as a vector.  The user must then store
        //    the entries of the four blocks of the matrix into the vector A.
        //    Each block is stored by columns.
        //
        //    A1, the banded portion of the matrix, is stored in
        //    the first (ML+MU+1)*N1 entries of A, using the obvious variant
        //    of the LINPACK general band format.
        //
        //    The following formulas should be used to determine how to store
        //    the entry corresponding to row I and column J in the original matrix:
        //
        //    Entries of A1:
        //
        //      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
        //
        //      Store the I, J entry into location
        //      (I-J+MU+1)+(J-1)*(ML+MU+1).
        //
        //    Entries of A2:
        //
        //      1 <= I <= N1, N1+1 <= J <= N1+N2.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+(J-N1-1)*N1+I.
        //
        //    Entries of A3:
        //
        //      N1+1 <= I <= N1+N2, 1 <= J <= N1.
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //
        //    Entries of A4:
        //
        //      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
        //
        //      Store the I, J entry into location
        //      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
        //      (same formula used for A3).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 September 2003
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
        //    Output, double R8CBB_ZERO[(ML+MU+1)*N1 + 2*N1*N2 + N2*N2], the matrix.
        //
    {
        double[] a;

        a = r8vec_zeros_new((ml + mu + 1) * n1 + 2 * n1 * n2 + n2 * n2);

        return a;
    }

}