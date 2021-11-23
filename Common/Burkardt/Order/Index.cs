using System;
using Burkardt.Types;

namespace Burkardt.IndexNS;

public static class Index
{
    public static int index0(int i_min, int i, int i_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX0 indexes a 1D vector using a zero base.
        //
        //  Discussion:
        //
        //    Index       Element
        //    ---------   --------
        //    0           I_MIN
        //    INDEX0      I
        //   (INDEX_MAX)  I_MAX
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I, I_MAX, for the first index,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX0, the index of element I.
        //
    {
        int value = i - i_min;

        return value;
    }

    public static int index01(int i_min, int i, int i_max, int j_min, int j, int j_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX01 indexes a 2D array by columns, with a zero base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
        //    and increasing the row index first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I, I_MAX, for row indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int J_MIN, J, J_MAX, for column indices,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX01, the index of element (I,J).
        //
    {
        int value = i - i_min
                    + (i_max + 1 - i_min) * (j - j_min);

        return value;
    }

    public static int index012(int i_min, int i, int i_max, int j_min, int j, int j_max,
            int k_min, int k, int k_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX012 indexes a 3D array by columns with zero base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at entry (I_MIN,J_MIN,K_MIN), 
        //    and increasing the row index first, then the column index.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I, I_MAX, for row indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int J_MIN, J, J_MAX, for column indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int K_MIN, K, K_MAX, for plane indices,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX012, the index of element (I,J,K).
        //
    {
        int value = + (i - i_min)
                    + (i_max + 1 - i_min) * (j - j_min) *
                    +(i_max + 1 - i_min) * (j_max + 1 - j_min) * (k - k_min);

        return value;
    }

    public static int index0123(int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max,
            int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX0123 indexes a 4D array by columns, with a zero base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
        //    and increasing the initial index first, then the second, third and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I1_MIN, I1, I1_MAX, for index 1,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I2_MIN, I2, I2_MAX, for index 2,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I3_MIN, I3, I3_MAX, for index 3,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I4_MIN, I4, I4_MAX, for index 4,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX0123, the index of (I1,I2,I3,I4).
        //
    {
        int value = + (i1 - i1_min)
                    + (i1_max + 1 - i1_min) * (i2 - i2_min)
                    + (i1_max + 1 - i1_min) * (i2_max + 1 - i2_min)
                                            * (i3 - i3_min)
                    + (i1_max + 1 - i1_min) * (i2_max + 1 - i2_min)
                                            * (i3_max + 1 - i3_min) * (i4 - i4_min);

        return value;
    }

    public static int index0n(int n, int[] i_min, int[] i, int[] i_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX0N indexes an N-dimensional array by columns, with zero base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at entry 
        //      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
        //    and increasing the first index up to I_MAX(1), 
        //    then the second and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of indices.
        //
        //    Input, int I_MIN[N], the minimum indices.
        //
        //    Input, int I[N], the indices.
        //
        //    Input, int I_MAX[N], for maximum indices.
        //
        //    Output, int INDEX0N, the index of element I.
        //
    {
        int j;

        int value = i[n - 1] - i_min[n - 1];

        for (j = n - 2; 0 <= j; j--)
        {
            value = value * (i_max[j] + 1 - i_min[j]) + (i[j] - i_min[j]);
        }
        
        return value;
    }

    public static int index1(int i_min, int i, int i_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX1 indexes a 1D vector using a unit base.
        //
        //  Discussion:
        //
        //    Index       Element
        //    ---------   --------
        //    1           I_MIN
        //    INDEX1      I
        //   (INDEX_MAX)  I_MAX
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I, I_MAX, for the first index,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX1, the index of element I.
        //
    {
        const int index_min = 1;

        int value = index_min + (i - i_min);

        return value;
    }

    public static int index10(int i_min, int i, int i_max, int j_min, int j, int j_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX10 indexes a 2D array by rows, with a zero base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
        //    and increasing the column index first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I, I_MAX, for row indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int J_MIN, J, J_MAX, for column indices,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX10, the index of element (I,J).
        //
    {
        int value = j - j_min
                    + (i - i_min) * (j_max + 1 - j_min);

        return value;
    }

    public static int index12(int i_min, int i, int i_max, int j_min, int j, int j_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX12 indexes a 2D array by columns, with a unit base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
        //    and increasing the row index first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I, I_MAX, for row indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int J_MIN, J, J_MAX, for column indices,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX12, the index of element (I,J).
        //
    {
        const int index_min = 1;

        int value = index_min
                    + (i - i_min)
                    + (i_max + 1 - i_min) * (j - j_min);

        return value;
    }

    public static int index123(int i_min, int i, int i_max, int j_min, int j, int j_max,
            int k_min, int k, int k_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX123 indexes a 3D array by columns with unit base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at entry (I_MIN,J_MIN,K_MIN), 
        //    and increasing the row index first, then the column index.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I, I_MAX, for row indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int J_MIN, J, J_MAX, for column indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int K_MIN, K, K_MAX, for plane indices,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX123, the index of element (I,J,K).
        //
    {
        int index_min = 1;

        int value = index_min
                    + (i - i_min)
                    + (i_max + 1 - i_min) * (j - j_min) *
                    +(i_max + 1 - i_min) * (j_max + 1 - j_min) * (k - k_min);

        return value;
    }

    public static int index1234(int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max,
            int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX1234 indexes a 4D array by columns, with a unit base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
        //    and increasing the initial index first, then the second, third and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I1_MIN, I1, I1_MAX, for index 1,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I2_MIN, I2, I2_MAX, for index 2,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I3_MIN, I3, I3_MAX, for index 3,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I4_MIN, I4, I4_MAX, for index 4,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX1234, the index of (I1,I2,I3,I4).
        //
    {
        const int index_min = 1;

        int value = index_min
                    + (i1 - i1_min)
                    + (i1_max + 1 - i1_min) * (i2 - i2_min)
                    + (i1_max + 1 - i1_min) * (i2_max + 1 - i2_min)
                                            * (i3 - i3_min)
                    + (i1_max + 1 - i1_min) * (i2_max + 1 - i2_min)
                                            * (i3_max + 1 - i3_min) * (i4 - i4_min);

        return value;
    }

    public static int index1n(int n, int[] i_min, int[] i, int[] i_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX1N indexes an N-dimensional array by columns, with unit base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at entry 
        //      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
        //    and increasing the first index up to I_MAX(1), 
        //    then the second and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of indices.
        //
        //    Input, int I_MIN[N], the minimum indices.
        //
        //    Input, int I[N], the indices.
        //
        //    Input, int I_MAX[N], for maximum indices.
        //
        //    Output, int INDEX1N, the index of element I.
        //
    {
        const int index_min = 1;
        int j;

        int value = i[n - 1] - i_min[n - 1];

        for (j = n - 2; 0 <= j; j--)
        {
            value = value * (i_max[j] + 1 - i_min[j]) + (i[j] - i_min[j]);
        }

        value += index_min;

        return value;
    }

    public static int index21(int i_min, int i, int i_max, int j_min, int j, int j_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX21 indexes a 2D array by rows, with a unit base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
        //    and increasing the column index first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I, I_MAX, for row indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int J_MIN, J, J_MAX, for column indices,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX21, the index of element (I,J).
        //
    {
        const int index_min = 1;

        int value = index_min
                    + (j - j_min)
                    + (i - i_min) * (j_max + 1 - j_min);

        return value;
    }

    public static int index210(int i_min, int i, int i_max, int j_min, int j, int j_max,
            int k_min, int k, int k_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX210 indexes a 3D array by rows, with zero base.
        //
        //  Discussion:
        //
        //    When we say "by rows", we really just mean that entries of the array are 
        //    indexed starting at entry (I_MIN,J_MIN,K_MIN), and the increasing the LAST
        //    index first, then the next-to-the-last, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I, I_MAX, for row indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int J_MIN, J, J_MAX, for column indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int K_MIN, K, K_MAX, for plane indices,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX210, the index of element (I,J,K).
        //
    {
        int value = + (k - k_min)
                    + (j - j_min) * (k_max + 1 - k_min)
                    + (i - i_min) * (j_max + 1 - j_min) * (k_max + 1 - k_min);

        return value;
    }

    public static int index321(int i_min, int i, int i_max, int j_min, int j, int j_max,
            int k_min, int k, int k_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX321 indexes a 3D array by rows, with zero base.
        //
        //  Discussion:
        //
        //    When we say "by rows", we really just mean that entries of the array are 
        //    indexed starting at entry (I_MIN,J_MIN,K_MIN), and the increasing the LAST
        //    index first, then the next-to-the-last, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I_MIN, I, I_MAX, for row indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int J_MIN, J, J_MAX, for column indices,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int K_MIN, K, K_MAX, for plane indices,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX321, the index of element (I,J,K).
        //
    {
        const int index_min = 1;

        int value = index_min
                    + (k - k_min)
                    + (j - j_min) * (k_max + 1 - k_min)
                    + (i - i_min) * (j_max + 1 - j_min) * (k_max + 1 - k_min);

        return value;
    }

    public static int index3210(int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max,
            int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX3210 indexes a 4D array by rows, with zero base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
        //    and increasing the last index, then the next to last, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I1_MIN, I1, I1_MAX, for index 1,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I2_MIN, I2, I2_MAX, for index 2,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I3_MIN, I3, I3_MAX, for index 3,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I4_MIN, I4, I4_MAX, for index 4,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX3210, the index of (I1,I2,I3,I4).
        //
    {
        int value = + (i4 - i4_min)
                    + (i3 - i3_min)
                    * (i4_max + 1 - i4_min)
                    + (i2 - i2_min) * (i3_max + 1 - i3_min)
                                    * (i4_max + 1 - i4_min)
                    + (i1 - i1_min) * (i2_max + 1 - i2_min) * (i3_max + 1 - i3_min)
                    * (i4_max + 1 - i4_min);

        return value;
    }

    public static int index4321(int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max,
            int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX4321 indexes a 4D array by rows, with unit base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
        //    and increasing the last index, then the next to last, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I1_MIN, I1, I1_MAX, for index 1,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I2_MIN, I2, I2_MAX, for index 2,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I3_MIN, I3, I3_MAX, for index 3,
        //    the minimum, the index, and the maximum.
        //
        //    Input, int I4_MIN, I4, I4_MAX, for index 4,
        //    the minimum, the index, and the maximum.
        //
        //    Output, int INDEX4321, the index of (I1,I2,I3,I4).
        //
    {
        int index_min = 1;

        int value = index_min
                    + (i4 - i4_min)
                    + (i3 - i3_min)
                    * (i4_max + 1 - i4_min)
                    + (i2 - i2_min) * (i3_max + 1 - i3_min)
                                    * (i4_max + 1 - i4_min)
                    + (i1 - i1_min) * (i2_max + 1 - i2_min) * (i3_max + 1 - i3_min)
                    * (i4_max + 1 - i4_min);

        return value;
    }

    public static int indexn0(int n, int[] i_min, int[] i, int[] i_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEXN0 indexes an N-dimensional array by rows, with zero base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at entry 
        //      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
        //    and increasing the last index up to I_MAX(N), 
        //    then the next-to-last and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of indices.
        //
        //    Input, int I_MIN[N], the minimum indices.
        //
        //    Input, int I[N], the indices.
        //
        //    Input, int I_MAX[N], for maximum indices.
        //
        //    Output, int INDEXN0, the index of element I.
        //
    {
        int j;

        int value = i[0] - i_min[0];

        for (j = 1; j < n; j++)
        {
            value = value * (i_max[j] + 1 - i_min[j]) + (i[j] - i_min[j]);
        }
        
        return value;
    }

    public static int indexn1(int n, int[] i_min, int[] i, int[] i_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEXN1 indexes an N-dimensional array by rows, with unit base.
        //
        //  Discussion:
        //
        //    Entries of the array are indexed starting at entry 
        //      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
        //    and increasing the last index up to I_MAX(N), 
        //    then the next-to-last and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of indices.
        //
        //    Input, int I_MIN[N], the minimum indices.
        //
        //    Input, int I[N], the indices.
        //
        //    Input, int I_MAX[N], for maximum indices.
        //
        //    Output, int INDEXN1, the index of element I.
        //
    {
        const int index_min = 1;
        int j;

        int value = i[0] - i_min[0];

        for (j = 1; j < n; j++)
        {
            value = value * (i_max[j] + 1 - i_min[j]) + (i[j] - i_min[j]);
        }

        value += index_min;

        return value;
    }

    public static void index_box2_next_2d(int n1, int n2, int ic, int jc, ref int i, ref int j,
            ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_BOX2_NEXT_2D produces index vectors on the surface of a box in 2D.
        //
        //  Discussion:
        //
        //    The box has center at (IC,JC), and half-widths N1 and N2.
        //    The index vectors are exactly those which are between (IC-N1,JC-N1) and
        //    (IC+N1,JC+N2) with the property that at least one of I and J
        //    is an "extreme" value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the half-widths of the box, that is, the
        //    maximum distance allowed between (IC,JC) and (I,J).
        //
        //    Input, int IC, JC, the central cell of the box.
        //
        //    Input/output, int &I, &J.  On input, the previous index set.
        //    On output, the next index set.  On the first call, MORE should
        //    be set to FALSE, and the input values of I and J are ignored.
        //
        //    Input/output, bool &MORE.
        //    On the first call for a given box, the user should set MORE to FALSE.
        //    On return, the routine sets MORE to TRUE.
        //    When there are no more indices, the routine sets MORE to FALSE.
        //
    {
        switch (more)
        {
            case false:
                more = true;
                i = ic - n1;
                j = jc - n2;
                return;
        }

        if (i == ic + n1 && j == jc + n2)
        {
            more = false;
            return;
        }

        //
        //  Increment J.
        //
        j += 1;
        //
        //  Check J.
        //
        if (jc + n2 < j)
        {
            j = jc - n2;
            i += 1;
        }
        else if (j < jc + n2 && (i == ic - n1 || i == ic + n1))
        {
        }
        else
        {
            j = jc + n2;
        }

    }

    public static void index_box2_next_3d(int n1, int n2, int n3, int ic, int jc, int kc,
            ref int i, ref int j, ref int k, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_BOX2_NEXT_3D produces index vectors on the surface of a box in 3D.
        //
        //  Discussion:
        //
        //    The box has a central cell of (IC,JC,KC), with half widths of
        //    (N1,N2,N3).  The index vectors are exactly those between
        //    (IC-N1,JC-N2,KC-N3) and (IC+N1,JC+N2,KC+N3) with the property that 
        //    at least one of I, J, and K is an "extreme" value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the "half widths" of the box, that is, the
        //    maximum distances from the central cell allowed for I, J and K.
        //
        //    Input, int IC, JC, KC, the central cell of the box.
        //
        //    Input/output, int &I, &J, &K.  On input, the previous index set.
        //    On output, the next index set.  On the first call, MORE should
        //    be set to FALSE, and the input values of I, J, and K are ignored.
        //
        //    Input/output, bool &MORE.
        //    On the first call for a given box, the user should set MORE to FALSE.
        //    On return, the routine sets MORE to TRUE.
        //    When there are no more indices, the routine sets MORE to FALSE.
        //
    {
        switch (more)
        {
            case false:
                more = true;
                i = ic - n1;
                j = jc - n2;
                k = kc - n3;
                return;
        }

        if (i == ic + n1 && j == jc + n2 && k == kc + n3)
        {
            more = false;
            return;
        }

        //
        //  Increment K.
        //
        k += 1;
        //
        //  Check K.
        //
        if (kc + n3 < k)
        {
            k = kc - n3;
            j += 1;
        }
        else if (k < kc + n3 &&
                 (i == ic - n1 || i == ic + n1 ||
                  j == jc - n2 || j == jc + n2))
        {
            return;
        }
        else
        {
            k = kc + n3;
            return;
        }

        //
        //  Check J.
        //
        if (jc + n2 < j)
        {
            j = jc - n2;
            i += 1;
        }
        else if (j < jc + n2 &&
                 (i == ic - n1 || i == ic + n1 ||
                  k == kc - n3 || k == kc + n3))
        {
        }
        else
        {
            j = jc + n2;
        }

    }

    public static void index_box_next_2d(int n1, int n2, ref int i, ref int j, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_BOX_NEXT_2D produces index vectors on the surface of a box in 2D.
        //
        //  Discussion:
        //
        //    The index vectors are exactly those which are between (1,1) and
        //    (N1,N2) with the property that at least one of I and J
        //    is an "extreme" value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the "dimensions" of the box, that is, the
        //    maximum values allowed for I and J.  The minimum values are
        //    assumed to be 1.
        //
        //    Input/output, int &I, &J.  On input, the previous index set.
        //    On output, the next index set.  On the first call, MORE should
        //    be set to FALSE, and the input values of I and J are ignored.
        //
        //    Input/output, bool &MORE.
        //    On the first call for a given box, the user should set MORE to FALSE.
        //    On return, the routine sets MORE to TRUE.
        //    When there are no more indices, the routine sets MORE to FALSE.
        //
    {
        switch (more)
        {
            case false:
                more = true;
                i = 1;
                j = 1;
                return;
        }

        if (i == n1 && j == n2)
        {
            more = false;
            return;
        }

        //
        //  Increment J.
        //
        j += 1;
        //
        //  Check J.
        //
        if (n2 < j)
        {
            j = 1;
            i += 1;
        }
        else if (j < n2 && (i == 1 || i == n1))
        {
        }
        else
        {
            j = n2;
        }

    }

    public static void index_box_next_3d(int n1, int n2, int n3, ref int i, ref int j, ref int k,
            ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_BOX_NEXT_3D produces index vectors on the surface of a box in 3D.
        //
        //  Discussion:
        //
        //    The index vectors are exactly those which are between (1,1,1) and
        //    (N1,N2,N3) with the property that at least one of I, J, and K
        //    is an "extreme" value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the "dimensions" of the box, that is, the
        //    maximum values allowed for I, J and K.  The minimum values are
        //    assumed to be 1.
        //
        //    Input/output, int &I, &J, &K.  On input, the previous index set.
        //    On output, the next index set.  On the first call, MORE should
        //    be set to FALSE, and the input values of I, J, and K are ignored.
        //
        //    Input/output, bool &MORE.
        //    On the first call for a given box, the user should set MORE to FALSE.
        //    On return, the routine sets MORE to TRUE.
        //    When there are no more indices, the routine sets MORE to FALSE.
        //
    {
        switch (more)
        {
            case false:
                more = true;
                i = 1;
                j = 1;
                k = 1;
                return;
        }

        if (i == n1 && j == n2 && k == n3)
        {
            more = false;
            return;
        }

        //
        //  Increment K.
        //
        k += 1;
        //
        //  Check K.
        //
        if (n3 < k)
        {
            k = 1;
            j += 1;
        }
        else if (k < n3 && (i == 1 || i == n1 || j == 1 || j == n2))
        {
            return;
        }
        else
        {
            k = n3;
            return;
        }

        //
        //  Check J.
        //
        if (n2 < j)
        {
            j = 1;
            i += 1;
        }
        else if (j < n2 && (i == 1 || i == n1 || k == 1 || k == n3))
        {
        }
        else
        {
            j = n2;
        }

    }

    public static void index_next0(int n, int hi, ref int[] a, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_NEXT0 generates all index vectors within given upper limits.
        //
        //  Discussion:
        //
        //    The index vectors are generated in such a way that the reversed
        //    sequences are produced in lexicographic order.
        //
        //  Example:
        //
        //    N = 3,
        //    HI = 3
        //
        //    1   2   3
        //    ---------
        //    1   1   1
        //    2   1   1
        //    3   1   1
        //    1   2   1
        //    2   2   1
        //    3   2   1
        //    1   3   1
        //    2   3   1
        //    3   3   1
        //    1   1   2
        //    2   1   2
        //    3   1   2
        //    1   2   2
        //    2   2   2
        //    3   2   2
        //    1   3   2
        //    2   3   2
        //    3   3   2
        //    1   1   3
        //    2   1   3
        //    3   1   3
        //    1   2   3
        //    2   2   3
        //    3   2   3
        //    1   3   3
        //    2   3   3
        //    3   3   3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, int HI, the upper limit for the array indices.
        //    The lower limit is implicitly 1 and HI must be at least 1.
        //
        //    Input/output, int A[N].
        //    On startup calls, with MORE = FALSE, the input value of A
        //    doesn't matter, because the routine initializes it.
        //    On calls with MORE = TRUE, the input value of A must be
        //    the output value of A from the previous call.  (In other words,
        //    just leave it alone!).
        //    On output, A contains the successor set of indices to the input
        //    value.
        //
        //    Input/output, bool &MORE.  Set this variable FALSE before
        //    the first call.  Normally, MORE will be returned TRUE but
        //    once all the vectors have been generated, MORE will be
        //    reset FALSE and you should stop calling the program.
        //
    {
        int i;

        switch (more)
        {
            case false:
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = 1;
                }

                switch (hi)
                {
                    case < 1:
                        more = false;
                        Console.WriteLine("");
                        Console.WriteLine("INDEX_NEXT0 - Fatal error!");
                        Console.WriteLine("  HI is " + hi + "");
                        Console.WriteLine("  but HI must be at least 1.");
                        return;
                }

                break;
            }
            default:
            {
                int inc = 0;

                while (hi <= a[inc])
                {
                    a[inc] = 1;
                    inc += 1;
                }

                a[inc] += 1;
                break;
            }
        }

        //
        //  See if there are more entries to compute.
        //
        more = false;

        for (i = 0; i < n; i++)
        {
            if (a[i] < hi)
            {
                more = true;
            }
        }

    }

    public static void index_next1(int n, int[] hi, ref int[] a, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_NEXT1 generates all index vectors within given upper limits.
        //
        //  Discussion:
        //
        //    The index vectors are generated in such a way that the reversed
        //    sequences are produced in lexicographic order.
        //
        //  Example:
        //
        //    N = 3,
        //    HI(1) = 4, HI(2) = 2, HI(3) = 3
        //
        //    1   2   3
        //    ---------
        //    1   1   1
        //    2   1   1
        //    3   1   1
        //    4   1   1
        //    1   2   1
        //    2   2   1
        //    3   2   1
        //    4   2   1
        //    1   1   2
        //    2   1   2
        //    3   1   2
        //    4   1   2
        //    1   2   2
        //    2   2   2
        //    3   2   2
        //    4   2   2
        //    1   1   3
        //    2   1   3
        //    3   1   3
        //    4   1   3
        //    1   2   3
        //    2   2   3
        //    3   2   3
        //    4   2   3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, int HI[N], the upper limits for the array indices.
        //    The lower limit is implicitly 1, and each HI(I) should be at least 1.
        //
        //    Input/output, int A[N].
        //    On startup calls, with MORE = FALSE, the input value of A
        //    doesn't matter, because the routine initializes it.
        //    On calls with MORE = TRUE, the input value of A must be
        //    the output value of A from the previous call.  (In other words,
        //    just leave it alone!).
        //    On output, A contains the successor set of indices to the input
        //    value.
        //
        //    Input/output, bool &MORE.  Set this variable FALSE before
        //    the first call.  Normally, MORE will be returned TRUE but
        //    once all the vectors have been generated, MORE will be
        //    reset FALSE and you should stop calling the program.
        //
    {
        int i;

        switch (more)
        {
            case false:
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = 1;
                }

                for (i = 0; i < n; i++)
                {
                    switch (hi[i])
                    {
                        case < 1:
                            more = false;
                            Console.WriteLine("");
                            Console.WriteLine("INDEX_NEXT1 - Fatal error!");
                            Console.WriteLine("  Entry " + i + " of HI is " + hi[i] + "");
                            Console.WriteLine("  but all entries must be at least 1.");
                            return;
                    }
                }

                break;
            }
            default:
            {
                int inc = 0;

                while (hi[inc] <= a[inc])
                {
                    a[inc] = 1;
                    inc += 1;
                }

                a[inc] += 1;
                break;
            }
        }

        //
        //  See if there are more entries to compute.
        //
        more = false;

        for (i = 0; i < n; i++)
        {
            if (a[i] < hi[i])
            {
                more = true;
            }
        }

    }

    public static void index_next2(int n, int[] lo, int[] hi, ref int[] a, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_NEXT2 generates all index vectors within given lower and upper limits.
        //
        //  Example:
        //
        //    N = 3,
        //    LO(1) = 1, LO(2) = 10, LO(3) = 4
        //    HI(1) = 2, HI(2) = 11, HI(3) = 6
        //
        //    1   2   3
        //    ---------
        //    1  10   4
        //    2  10   4
        //    1  11   4
        //    2  11   4
        //    1  10   5
        //    2  10   5
        //    1  11   5
        //    2  11   5
        //    1  10   6
        //    2  10   6
        //    1  11   6
        //    2  11   6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.  The rank of
        //    the object being indexed.
        //
        //    Input, int LO[N], HI[N], the lower and upper limits for the array
        //    indices.  LO(I) should be less than or equal to HI(I), for each I.
        //
        //    Input/output, int A[N].
        //    On startup calls, with MORE = FALSE, the input value of A
        //    doesn't matter, because the routine initializes it.
        //    On calls with MORE = TRUE, the input value of A must be
        //    the output value of A from the previous call.  (In other words,
        //    just leave it alone!).
        //    On output, A contains the successor set of indices to the input
        //    value.
        //
        //    Input/output, bool &MORE.  Set this variable FALSE before
        //    the first call.  Normally, MORE will be returned TRUE but
        //    once all the vectors have been generated, MORE will be
        //    reset FALSE and you should stop calling the program.
        //
    {
        int i;

        switch (more)
        {
            case false:
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = lo[i];
                }

                for (i = 0; i < n; i++)
                {
                    if (hi[i] < lo[i])
                    {
                        more = false;
                        Console.WriteLine("");
                        Console.WriteLine("INDEX_NEXT2 - Fatal error!");
                        Console.WriteLine("  Entry " + i + " of HI is " + hi[i] + "");
                        Console.WriteLine("  Entry " + i + " of LO is " + lo[i] + "");
                        Console.WriteLine("  but LO(I) <= HI(I) is required.");
                        return;
                    }
                }

                break;
            }
            default:
            {
                int inc = 0;

                while (hi[inc] <= a[inc])
                {
                    a[inc] = lo[inc];
                    inc += 1;
                }

                a[inc] += 1;
                break;
            }
        }

        //
        //  See if there are more entries to compute.
        //
        more = false;

        for (i = 0; i < n; i++)
        {
            if (a[i] < hi[i])
            {
                more = true;
            }
        }
    }

    public static int index_rank0(int n, int hi, int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_RANK0 ranks an index vector within given upper limits.
        //
        //  Example:
        //
        //    N = 3,
        //    HI = 3
        //    A = ( 3, 1, 2 )
        //
        //    RANK = 12
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, int HI, the upper limit for the array indices.
        //    The lower limit is implicitly 1, and HI should be at least 1.
        //
        //    Input, int A[N], the index vector to be ranked.
        //
        //    Output, int INDEX_RANK0, the rank of the index vector, or -1 if A
        //    is not a legal index.
        //
    {
        int i;

        int rank = -1;
        for (i = 0; i < n; i++)
        {
            if (a[i] < 1 || hi < a[i])
            {
                return rank;
            }
        }

        rank = 0;
        for (i = n - 1; 0 <= i; i--)
        {
            rank = hi * rank + a[i];
        }

        rank = 1;
        int range = 1;
        for (i = 0; i < n; i++)
        {
            rank += (a[i] - 1) * range;
            range *= hi;
        }

        return rank;
    }

    public static int index_rank1(int n, int[] hi, int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_RANK1 ranks an index vector within given upper limits.
        //
        //  Example:
        //
        //    N = 3,
        //    HI(1) = 4, HI(2) = 2, HI(3) = 3
        //    A = ( 4, 1, 2 )
        //
        //    RANK = 12
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, int HI[N], the upper limits for the array indices.
        //    The lower limit is implicitly 1, and each HI(I) should be at least 1.
        //
        //    Input, int A[N], the index to be ranked.
        //
        //    Output, int INDEX_RANK1, the rank of the index vector, or -1 if A
        //    is not a legal index.
        //
    {
        int i;

        int rank = -1;
        for (i = 0; i < n; i++)
        {
            if (a[i] < 1 || hi[i] < a[i])
            {
                return rank;
            }
        }

        rank = 0;
        for (i = n - 1; 0 <= i; i--)
        {
            rank = hi[i] * rank + a[i];
        }

        rank = 1;
        int range = 1;
        for (i = 0; i < n; i++)
        {
            rank += (a[i] - 1) * range;
            range *= hi[i];
        }

        return rank;
    }

    public static int index_rank2(int n, int[] lo, int[] hi, int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_RANK2 ranks an index vector within given lower and upper limits.
        //
        //  Example:
        //
        //    N = 3,
        //    LO(1) = 1, LO(2) = 10, LO(3) = 4
        //    HI(1) = 2, HI(2) = 11, HI(3) = 6
        //    A = ( 1, 11, 5 )
        //
        //    RANK = 7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, int LO[N], HI[N], the lower and upper limits for the array
        //    indices.  LO(I) should be less than or equal to HI(I), for each I.
        //
        //    Input, int A[N], the index vector to be ranked.
        //
        //    Output, int INDEX_RANK2, the rank of the index vector, or -1 if A
        //    is not a legal index vector.
        //
    {
        int i;
        int rank;

        for (i = 0; i < n; i++)
        {
            if (a[i] >= lo[i] && hi[i] >= a[i])
            {
                continue;
            }

            rank = -1;
            return rank;
        }

        rank = 1;
        int range = 1;
        for (i = 0; i < n; i++)
        {
            rank += (a[i] - lo[i]) * range;
            range *= hi[i] + 1 - lo[i];
        }

        return rank;
    }

    public static void index_unrank0(int n, int hi, int rank, ref int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_UNRANK0 unranks an index vector within given upper limits.
        //
        //  Example:
        //
        //    N = 3,
        //    HI = 3
        //    RANK = 12
        //
        //    A = ( 3, 1, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, int HI, the upper limit for the array indices.
        //    The lower limit is implicitly 1, and HI should be at least 1.
        //
        //    Input, int RANK, the rank of the desired index vector.
        //
        //    Output, int A[N], the index vector of the given rank.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i] = 0;
        }

        switch (rank)
        {
            //
            //  The rank might be too small.
            //
            case < 1:
                return;
        }

        int range = (int)Math.Pow(hi, n);
        //
        //  The rank might be too large.
        //
        if (range < rank)
        {
            return;
        }

        int k = rank - 1;

        for (i = n - 1; 0 <= i; i--)
        {
            range /= hi;
            int j = k / range;
            a[i] = j + 1;
            k -= j * range;
        }

    }

    public static void index_unrank1(int n, int[] hi, int rank, ref int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_UNRANK1 unranks an index vector within given upper limits.
        //
        //  Example:
        //
        //    N = 3,
        //    HI(1) = 4, HI(2) = 2, HI(3) = 3
        //    RANK = 11
        //
        //    A = ( 3, 1, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, int HI[N], the upper limits for the array indices.
        //    The lower limit is implicitly 1, and each HI(I) should be at least 1.
        //
        //    Input, int RANK, the rank of the desired index vector.
        //
        //    Output, int A[N], the index vector of the given rank.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i] = 0;
        }

        switch (rank)
        {
            //
            //  The rank might be too small.
            //
            case < 1:
                return;
        }

        int range = typeMethods.i4vec_product(n, hi);
        //
        //  The rank might be too large.
        //
        if (range < rank)
        {
            return;
        }

        int k = rank - 1;

        for (i = n - 1; 0 <= i; i--)
        {
            range /= hi[i];
            int j = k / range;
            a[i] = j + 1;
            k -= j * range;
        }

    }

    public static void index_unrank2(int n, int[] lo, int[] hi, int rank, ref int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_UNRANK2 unranks an index vector within given lower and upper limits.
        //
        //  Example:
        //
        //    N = 3,
        //    LO(1) = 1, LO(2) = 10, LO(3) = 4
        //    HI(1) = 2, HI(2) = 11, HI(3) = 6
        //    RANK = 7
        //
        //    A = ( 1, 11, 5 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, int LO[N], HI[N], the lower and upper limits for the array
        //    indices.  It should be the case that LO(I) <= HI(I) for each I.
        //
        //    Input, int RANK, the rank of the desired index.
        //
        //    Output, int A[N], the index vector of the given rank.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i] = 0;
        }

        switch (rank)
        {
            //
            //  The rank might be too small.
            //
            case < 1:
                return;
        }

        int range = 1;
        for (i = 0; i < n; i++)
        {
            range *= hi[i] + 1 - lo[i];
        }

        //
        //  The rank might be too large.
        //
        if (range < rank)
        {
            return;
        }

        int k = rank - 1;
        for (i = n - 1; 0 <= i; i--)
        {
            range /= hi[i] + 1 - lo[i];
            int j = k / range;
            a[i] = j + lo[i];
            k -= j * range;
        }

    }

}