using System;
using System.Linq;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static int r8vec_min_index ( int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MIN_INDEX returns the index of the minimum value in an R8VEC.
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
        //    02 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], the array.
        //
        //    Output, int R8VEC_MIN_INDEX, the index of the smallest entry.
        //
    {
        int min_index;

        switch (n)
        {
            case <= 0:
                min_index = -1;
                break;
            default:
            {
                min_index = 0;

                int i;
                for ( i = 1; i < n; i++ )
                {
                    if ( a[i] < a[min_index] )
                    {
                        min_index = i;
                    }
                }

                break;
            }
        }

        return min_index;
    }
        
    public static double r8vec_max(int n, double[] dvec)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MAX returns the value of the maximum element in an R8VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 1998
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double *RVEC, a pointer to the first entry of the array.
        //
        //    Output, double R8VEC_MAX, the value of the maximum element.  This
        //    is set to 0.0 if N <= 0.
        //
    {
        switch (dvec.Length)
        {
            case <= 0:
                return 0;
        }

        // Limit to the number of items in the array as a maximum
        n = Math.Min(n, dvec.Length);

        return n == dvec.Length ? dvec.Max() : dvec.Take(n).Max();
    }
        
    public static int r8vec_max_index ( int n, double[] a )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MAX_INDEX returns the index of the maximum value in an R8VEC.
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
        //    02 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], the array.
        //
        //    Output, int R8VEC_MAX_INDEX, the index of the largest entry.
        //
    {
        int max_index;

        switch (n)
        {
            case <= 0:
                max_index = -1;
                break;
            default:
            {
                max_index = 0;

                int i;
                for ( i = 1; i < n; i++ )
                {
                    if ( a[max_index] < a[i] )
                    {
                        max_index = i;
                    }
                }

                break;
            }
        }

        return max_index;
    }
        
    public static double r8vec_min(int n, double[] dvec)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MIN returns the value of the minimum element in an R8VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 1998
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double *RVEC, a pointer to the first entry of the array.
        //
        //    Output, double R8VEC_MIN, the value of the minimum element.  This
        //    is set to 0.0 if N <= 0.
        //
    {
        switch (dvec.Length)
        {
            case <= 0:
                return 0;
        }

        // Limit to the number of items in the array as a maximum
        n = Math.Min(n, dvec.Length);

        return n == dvec.Length ? dvec.Min() : dvec.Take(n).Min();
    }

    public static double r8vec_amax(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_AMAX returns the maximum absolute value in an R8VEC.
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
        //    18 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], the array.
        //
        //    Output, double AMAX, the value of the entry
        //    of largest magnitude.
        //
    {
        int i;

        double amax = 0.0;
        for (i = 0; i < n; i++)
        {
            if (amax < Math.Abs(a[i]))
            {
                amax = Math.Abs(a[i]);
            }
        }

        return amax;
    }

    public static int r8vec_amax_index(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_AMAX_INDEX returns the index of the maximum absolute value in an R8VEC.
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
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], the array.
        //
        //    Output, int R8VEC_AMAX_INDEX, the index of the entry of largest magnitude.
        //
    {
        int amax_index;

        switch (n)
        {
            case <= 0:
                amax_index = -1;
                break;
            default:
            {
                amax_index = 1;
                double amax = Math.Abs(a[0]);

                int i;
                for (i = 2; i <= n; i++)
                {
                    if (amax < Math.Abs(a[i - 1]))
                    {
                        amax_index = i;
                        amax = Math.Abs(a[i - 1]);
                    }
                }

                break;
            }
        }

        return amax_index;
    }

    public static double r8vec_amin(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_AMIN returns the minimum absolute value in an R8VEC.
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
        //    18 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], the array.
        //
        //    Output, double R8VEC_AMIN, the value of the entry
        //    of smallest magnitude.
        //
    {
        int i;

        double value = r8_huge();
        for (i = 0; i < n; i++)
        {
            if (Math.Abs(a[i]) < value)
            {
                value = Math.Abs(a[i]);
            }
        }

        return value;
    }

    public static int r8vec_amin_index(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_AMIN_INDEX returns the index of the minimum absolute value in an R8VEC.
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
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], the array.
        //
        //    Output, int R8VEC_AMIN_INDEX, the index of the entry of smallest magnitude.
        //
    {
        int amin_index;

        switch (n)
        {
            case <= 0:
                amin_index = -1;
                break;
            default:
            {
                amin_index = 1;
                double amin = Math.Abs(a[0]);

                int i;
                for (i = 2; i <= n; i++)
                {
                    if (Math.Abs(a[i - 1]) < amin)
                    {
                        amin_index = i;
                        amin = Math.Abs(a[i - 1]);
                    }
                }

                break;
            }
        }

        return amin_index;
    }
        
    public static int r8vec_max_abs_index ( int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MAX_ABS_INDEX returns the index of the maximum absolute value in an R8VEC.
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
        //    08 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], the array.
        //
        //    Output, int R8VEC_MAX_ABS_INDEX, the index of the entry of 
        //    largest absolute value.
        //
    {
        int max_index;

        switch (n)
        {
            case <= 0:
                max_index = -1;
                break;
            default:
            {
                max_index = 0;

                int i;
                for ( i = 1; i < n; i++ )
                {
                    if ( Math.Abs ( a[max_index] ) < Math.Abs ( a[i] ) )
                    {
                        max_index = i;
                    }
                }

                break;
            }
        }

        return max_index;
    }

    public static double r8vec_min_pos ( int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MIN_POS returns the minimum positive value of an R8VEC.
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
        //    08 November 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries.
        //
        //    Input, double A[N], the array.
        //
        //    Output, double R8VEC_MIN_POS, the smallest positive entry.
        //
    {
        int i;

        double value = r8_huge();

        for ( i = 0; i < n; i++ )
        {
            switch (a[i])
            {
                case > 0.0:
                {
                    if ( a[i] < value )
                    {
                        value = a[i];
                    }

                    break;
                }
            }
        }
        return value;
    }

}