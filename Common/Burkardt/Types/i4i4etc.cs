using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void i4i4_sort_a(int i1, int i2, ref int j1, ref int j2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4I4_SORT_A ascending sorts a pair of I4's.
        //
        //  Discussion:
        //
        //    The program allows the reasonable call:
        //
        //      i4i4_sort_a ( i1, i2, &i1, &i2 );
        //
        //    and this will return the reasonable result.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I1, I2, the values to sort.
        //
        //    Output, int J1, J2, the sorted values.
        //
    {
        int k1;
        int k2;
        //
        //  Copy arguments, so that the user can make "reasonable" calls like:
        //
        //    i4i4_sort_a ( i1, i2, &i1, &i2 );
        //
        k1 = i1;
        k2 = i2;

        j1 = Math.Min(k1, k2);
        j2 = Math.Max(k1, k2);
    }

    public static void i4i4i4_sort_a(int i1, int i2, int i3, ref int j1, ref int j2, ref int j3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4I4I4_SORT_A ascending sorts a triple of I4's.
        //
        //  Discussion:
        //
        //    The program allows the reasonable call:
        //
        //      i4i4i4_sort_a ( i1, i2, i3, &i1, &i2, &i3 );
        //
        //    and this will return the reasonable result.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I1, I2, I3, the values to sort.
        //
        //    Output, int *J1, *J2, *J3, the sorted values.
        //
    {
        int k1;
        int k2;
        int k3;
        //
        //  Copy arguments, so that the user can make "reasonable" calls like:
        //
        //    i4i4i4_sort_a ( i1, i2, i3, &i1, &i2, &i3 );
        //
        k1 = i1;
        k2 = i2;
        k3 = i3;

        j1 = Math.Min(Math.Min(k1, k2), Math.Min(k2, k3));
        j2 = Math.Min(Math.Max(k1, k2),
            Math.Min(Math.Max(k2, k3), Math.Max(k3, k1)));
        j3 = Math.Max(Math.Max(k1, k2), Math.Max(k2, k3));
    }
}