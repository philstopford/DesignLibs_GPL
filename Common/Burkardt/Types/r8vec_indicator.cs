namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_indicator ( int n, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDICATOR sets an R8VEC to the indicator vector {1,2,3...}.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Output, double A[N], the array to be initialized.
        //
    {
        int i;

        for ( i = 0; i <= n - 1; i++ )
        {
            a[i] = i + 1;
        }
    }

    public static double[] r8vec_indicator_new ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDICATOR_NEW sets an R8VEC to the indicator vector {1,2,3...}.
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
        //    20 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Output, double R8VEC_INDICATOR_NEW[N], the indicator array.
        //
    {
        int i;

        double[] a = new double[n];

        for ( i = 0; i <= n-1; i++ ) 
        {
            a[i] = i + 1;
        }

        return a;
    }

    public static void r8vec_indicator0(int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDICATOR0 sets an R8VEC to the indicator vector (0,1,2,...)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Output, double A[N], the array.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i] = i;
        }
    }

    public static double[] r8vec_indicator0_new(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDICATOR0_NEW sets an R8VEC to the indicator vector {0,1,2,...}.
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
        //    27 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Output, double R8VEC_INDICATOR0_NEW[N], the indicator array.
        //
    {
        int i;

        double[] a = new double[n];

        for (i = 0; i < n; i++)
        {
            a[i] = i;
        }

        return a;
    }

    public static void r8vec_indicator1(int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDICATOR1 sets an R8VEC to the indicator vector (1,2,3,...)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Output, double A[N], the array.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i] = i + 1;
        }

    }

    public static double[] r8vec_indicator1_new(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDICATOR1_NEW sets an R8VEC to the indicator vector {1,2,3,...}.
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
        //    27 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Output, double R8VEC_INDICATOR1_NEW[N], the indicator array.
        //
    {
        int i;

        double[] a = new double[n];

        for (i = 0; i < n; i++)
        {
            a[i] = i + 1;
        }

        return a;
    }

}