using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static bool r8vec_is_ascending(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_ASCENDING determines if an R8VEC is (weakly) ascending.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    For example, if:
        //
        //      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
        //
        //    then
        //
        //      R8VEC_IS_ASCENDING = TRUE
        //
        //    The sequence is not required to be strictly ascending.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the array.
        //
        //    Input, double X[N], the array to be examined.
        //
        //    Output, bool R8VEC_IS_ASCENDING, is TRUE if the
        //    entries ascend.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n - 1; i++)
        {
            if (!(x[i + 1] < x[i]))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_ascending_strictly(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_ASCENDING_STRICTLY determines if an R8VEC is strictly ascending.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    Notice the effect of entry number 6 in the following results:
        //
        //      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.4, 9.8 )
        //      Y = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
        //      Z = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.6, 9.8 )
        //
        //      R8VEC_IS_ASCENDING_STRICTLY ( X ) = FALSE
        //      R8VEC_IS_ASCENDING_STRICTLY ( Y ) = FALSE
        //      R8VEC_IS_ASCENDING_STRICTLY ( Z ) = TRUE
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the array.
        //
        //    Input, double X[N], the array to be examined.
        //
        //    Output, bool R8VEC_IS_ASCENDING_STRICTLY, is TRUE if the
        //    entries strictly ascend.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n - 1; i++)
        {
            if (!(x[i + 1] <= x[i]))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_binary(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_BINARY is true if the entries in an R8VEC are all 0 or 1.
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
        //    31 March 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double X[N], the vector to be checked.
        //
        //    Output, bool R8VEC_IS_BINARY is true if are entries are 0 or 1.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n; i++)
        {
            if (x[i] == 0.0 || Math.Abs(x[i] - 1.0) <= typeMethods.r8_epsilon())
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_distinct(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_DISTINCT is true if the entries in an R8VEC are distinct.
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
        //    05 April 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double X[N], the vector to be checked.
        //
        //    Output, bool R8VEC_IS_DISTINCT is true if all entries are distinct.
        //
    {
        int i;

        bool value = true;

        for (i = 1; i < n; i++)
        {
            int j;
            for (j = 0; j < i; j++)
            {
                if (!(Math.Abs(x[i] - x[j]) <= typeMethods.r8_epsilon()))
                {
                    continue;
                }

                value = false;
                break;
            }
        }

        return value;
    }

    public static bool r8vec_is_in_01(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_IN_01 is TRUE if the entries of an R8VEC are in the range [0,1].
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
        //    06 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries.
        //
        //    Input, double X[N], the vector
        //
        //    Output, bool R8VEC_IS_IN_01, is TRUE if every entry is
        //    between 0 and 1.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n; i++)
        {
            if (!(x[i] < 0.0) && !(1.0 < x[i]))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_in_ab(int n, double[] x, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_IN_AB is TRUE if the entries of an R8VEC are in the range [A,B].
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
        //    15 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries.
        //
        //    Input, double X[N], the vector
        //
        //    Input, double A, B, the limits of the range.
        //
        //    Output, bool R8VEC_IS_IN_AB, is TRUE if every entry is
        //    between A and B.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n; i++)
        {
            if (!(x[i] < a) && !(b < x[i]))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_insignificant(int n, double[] r, double[] s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_INSIGNIFICANT determines if an R8VEC is insignificant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vectors.
        //
        //    Input, double R[N], the vector to be compared against.
        //
        //    Input, double S[N], the vector to be compared.
        //
        //    Output, bool R8VEC_IS_INSIGNIFICANT, is TRUE if S is insignificant
        //    compared to R.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n; i++)
        {
            double t = r[i] + s[i];
            double tol = r8_epsilon() * Math.Abs(r[i]);

            if (!(tol < Math.Abs(r[i] - t)))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_integer(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_INTEGER is TRUE if an R8VEC is integral.
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
        //    03 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], the vector
        //
        //    Output, bool R8VEC_IS_INTEGER, is TRUE if every entry is an integer.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n; i++)
        {
            if (!(Math.Abs(a[i] - (int) a[i]) > typeMethods.r8_epsilon()))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_negative(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_NEGATIVE: all entries of R8VEC are strictly negative.
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
        //    24 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vector.
        //
        //    Input, double A[N], the vector.
        //
        //    Output, bool R8VEC_IS_NEGATIVE, is TRUE if all entries are strictly 
        //    negative.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n; i++)
        {
            if (!(0 <= a[i]))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_negative_any(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_NEGATIVE_ANY: ( any ( A < 0 ) ) for R8VEC's.
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
        //    09 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries.
        //
        //    Input, double A[N], the vector to check.
        //
        //    Output, bool R8VEC_IS_NEGATIVE_ANY is TRUE if any entry
        //    is less than zero.
        //
    {
        int i;

        bool value = false;

        for (i = 0; i < n; i++)
        {
            if (!(a[i] < 0.0))
            {
                continue;
            }

            value = true;
            break;
        }

        return value;
    }

    public static bool r8vec_is_nonnegative(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_NONNEGATIVE is true if all entries in an R8VEC are nonnegative.
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
        //    04 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double X[N], the vector to be checked.
        //
        //    Output, bool R8VEC_IS_NONNEGATIVE is true if all entries are nonnegative.
        //
    {
        int i;

        bool value = true;
        for (i = 0; i < n; i++)
        {
            if (!(x[i] < 0.0))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_nonpositive(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_NONPOSITIVE: ( all ( A <= 0 ) ) for R8VEC's.
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
        //    08 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries.
        //
        //    Input, double A[N], the vector to check.
        //
        //    Output, bool R8VEC_IS_NONPOSITIVE is TRUE if all entries
        //    are less than or equal to zero.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n; i++)
        {
            if (!(0.0 < a[i]))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_nonzero_any(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_NONZERO_ANY: ( any A nonzero ) for R8VEC's.
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
        //    25 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries.
        //
        //    Input, double A[N], the vector to check.
        //
        //    Output, bool R8VEC_IS_NONZERO_ANY is TRUE if any entry is nonzero.
        //
    {
        int i;

        bool value = false;

        for (i = 0; i < n; i++)
        {
            if (a[i] == 0.0)
            {
                continue;
            }

            value = true;
            break;
        }

        return value;
    }

    public static bool r8vec_is_one(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_ONE is true if the entries in an R8VEC are all 1.
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
        //    30 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double X[N], the vector to be checked.
        //
        //    Output, bool R8VEC_IS_ONE is true if all entries are 1.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n; i++)
        {
            if (!(Math.Abs(x[i] - 1.0) > typeMethods.r8_epsilon()))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_positive(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_POSITIVE: all entries of R8VEC are strictly positive.
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
        //    24 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vector.
        //
        //    Input, double A[N], the vector.
        //
        //    Output, bool R8VEC_IS_POSITIVE, is TRUE if all entries are strictly 
        //    positive.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n; i++)
        {
            if (!(a[i] <= 0.0))
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }

    public static bool r8vec_is_zero(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IS_ZERO is true if the entries in an R8VEC are all zero.
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
        //    30 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double X[N], the vector to be checked.
        //
        //    Output, bool R8VEC_IS_ZERO is true if all entries are zero.
        //
    {
        int i;

        bool value = true;

        for (i = 0; i < n; i++)
        {
            if (x[i] == 0.0)
            {
                continue;
            }

            value = false;
            break;
        }

        return value;
    }
}