namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_copy(int n, double[] a1, ref double[] a2, int a1index = 0, int a2index = 0)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_COPY copies an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input, double A1[N], the vector to be copied.
        //
        //    Output, double A2[N], the copy of A1.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a2[i + a2index] = a1[i + a1index];
        }
    }


    public static double[] r8vec_copy_new(int n, double[] a1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_COPY_NEW copies an R8VEC to a new R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input, double A1[N], the vector to be copied.
        //
        //    Output, double R8VEC_COPY_NEW[N], the copy of A1.
        //
    {
        int i;

        double[] a2 = new double[n];

        for (i = 0; i < n; i++)
        {
            a2[i] = a1[i];
        }

        return a2;
    }
        
}