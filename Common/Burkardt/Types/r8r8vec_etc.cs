using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8r8r8vec_index_insert_unique(int maxn, ref int n, ref double[] x, ref double[] y,
            ref double[] z, ref int[] indx, double xval, double yval, double zval, ref int ival,
            ref int ierror )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8R8 value in an indexed sorted list.
        //
        //  Discussion:
        //
        //    If the input value does not occur in the current list, it is added,
        //    and N, X, Y, Z and INDX are updated.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MAXN, the maximum size of the list.
        //
        //    Input/output, int &N, the size of the list.
        //
        //    Input/output, double X[N], Y[N], Z[N], the R8R8R8 vector.
        //
        //    Input/output, int INDX[N], the sort index of the list.
        //
        //    Input, double XVAL, YVAL, ZVAL, the value to be inserted
        //    if it is not already in the list.
        //
        //    Output, int &IVAL, the index in X, Y, Z corresponding to the
        //    value XVAL, YVAL, ZVAL.
        //
        //    Output, int &IERROR, 0 for no error, 1 if an error occurred.
        //
    {
        int equal = 0;
        int less = 0;
        int more = 0;

        ierror = 0;

        switch (n)
        {
            case <= 0 when maxn <= 0:
                ierror = 1;
                Console.WriteLine("");
                Console.WriteLine("R8R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!");
                Console.WriteLine("  Not enough space to store new data.");
                return;
            case <= 0:
                n = 1;
                x[0] = xval;
                y[0] = yval;
                z[0] = zval;
                indx[0] = 1;
                ival = 1;
                return;
        }

        //
        //  Does ( XVAL, YVAL, ZVAL ) already occur in ( X, Y, Z)?
        //
        r8r8r8vec_index_search(n, x, y, z, indx, xval, yval, zval,
            ref less, ref equal, ref more);

        switch (equal)
        {
            case 0 when maxn <= n:
                ierror = 1;
                Console.WriteLine("");
                Console.WriteLine("R8R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!");
                Console.WriteLine("  Not enough space to store new data.");
                return;
            case 0:
            {
                x[n] = xval;
                y[n] = yval;
                z[n] = zval;
                ival = n + 1;
                for (int i = n - 1; more - 1 <= i; i--)
                {
                    indx[i + 1] = indx[i];
                }

                indx[more - 1] = n + 1;
                n += 1;
                break;
            }
            default:
                ival = indx[equal - 1];
                break;
        }
    }

    public static void r8r8r8vec_index_search(int n, double[] x, double[] y, double[] z,
            int[] indx, double xval, double yval, double zval, ref int less, ref int equal,
            ref int more )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8R8R8VEC_INDEX_SEARCH searches for an R8R8R8 value in an indexed sorted list.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the list.
        //
        //    Input, double X[N], Y[N], Z[N], the list.
        //
        //    Input, int INDX[N], the sort index of the list.
        //
        //    Input, double XVAL, YVAL, ZVAL, the value to be sought.
        //
        //    Output, int &LESS, &EQUAL, &MORE, the indexes in INDX of the
        //    entries of X that are just less than, equal to, and just greater
        //    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
        //    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
        //    is the greatest entry of X, then MORE is N+1.
        //
    {
        switch (n)
        {
            case <= 0:
                less = 0;
                equal = 0;
                more = 0;
                return;
        }

        int lo = 1;
        int hi = n;

        double xlo = x[indx[lo - 1] - 1];
        double ylo = y[indx[lo - 1] - 1];
        double zlo = z[indx[lo - 1] - 1];

        double xhi = x[indx[hi - 1] - 1];
        double yhi = y[indx[hi - 1] - 1];
        double zhi = z[indx[hi - 1] - 1];

        int compare = r8r8r8_compare(xval, yval, zval, xlo, ylo, zlo);

        switch (compare)
        {
            case -1:
                less = 0;
                equal = 0;
                more = 1;
                return;
            case 0:
                less = 0;
                equal = 1;
                more = 2;
                return;
        }

        compare = r8r8r8_compare(xval, yval, zval, xhi, yhi, zhi);

        switch (compare)
        {
            case 1:
                less = n;
                equal = 0;
                more = n + 1;
                return;
            case 0:
                less = n - 1;
                equal = n;
                more = n + 1;
                return;
        }

        for (;;)
        {
            if (lo + 1 == hi)
            {
                less = lo;
                equal = 0;
                more = hi;
                return;
            }

            int mid = (lo + hi) / 2;
            double xmid = x[indx[mid - 1] - 1];
            double ymid = y[indx[mid - 1] - 1];
            double zmid = z[indx[mid - 1] - 1];

            compare = r8r8r8_compare(xval, yval, zval, xmid, ymid, zmid);

            switch (compare)
            {
                case 0:
                    equal = mid;
                    less = mid - 1;
                    more = mid + 1;
                    return;
                case -1:
                    hi = mid;
                    break;
                case +1:
                    lo = mid;
                    break;
            }
        }
    }

    public static void r8r8vec_index_insert_unique(int maxn, ref int n, ref double[] x, ref double[] y,
            ref int[] indx, double xval, double yval, ref int ival, ref int ierror )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8 value in an indexed sorted list.
        //
        //  Discussion:
        //
        //    If the input value does not occur in the current list, it is added,
        //    and N, X, Y and INDX are updated.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MAXN, the maximum size of the list.
        //
        //    Input/output, int &N, the size of the list.
        //
        //    Input/output, double X[N], Y[N], the list of R8R8 vectors.
        //
        //    Input/output, int INDX[N], the sort index of the list.
        //
        //    Input, double XVAL, YVAL, the value to be inserted if it is
        //    not already in the list.
        //
        //    Output, int &IVAL, the index in X, Y corresponding to the
        //    value XVAL, YVAL.
        //
        //    Output, int &IERROR, 0 for no error, 1 if an error occurred.
        //
    {
        int equal = 0;
        int less = 0;
        int more = 0;

        ierror = 0;

        switch (n)
        {
            case <= 0 when maxn <= 0:
                Console.WriteLine("");
                Console.WriteLine("R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!");
                Console.WriteLine("  Not enough space to store new data.");
                return;
            case <= 0:
                n = 1;
                x[0] = xval;
                y[0] = yval;
                indx[0] = 1;
                ival = 1;
                return;
        }

        //
        //  Does ( XVAL, YVAL ) already occur in ( X, Y )?
        //
        r8r8vec_index_search(n, x, y, indx, xval, yval, ref less, ref equal, ref more);

        switch (equal)
        {
            case 0 when maxn <= n:
                Console.WriteLine("");
                Console.WriteLine("R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!");
                Console.WriteLine("  Not enough space to store new data.");
                return;
            case 0:
            {
                x[n] = xval;
                y[n] = yval;
                ival = n + 1;
                for (int i = n - 1; more - 1 <= i; i--)
                {
                    indx[i + 1] = indx[i];
                }

                indx[more - 1] = n + 1;
                n += 1;
                break;
            }
            default:
                ival = indx[equal - 1];
                break;
        }
    }

    public static void r8r8vec_index_search(int n, double[] x, double[] y, int[] indx,
            double xval, double yval, ref int less, ref int equal, ref int more )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8R8VEC_INDEX_SEARCH searches for an R8R8 value in an indexed sorted list.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the current list.
        //
        //    Input, double X[N], Y[N], the list.
        //
        //    Input, int INDX[N], the sort index of the list.
        //
        //    Input, double XVAL, YVAL, the value to be sought.
        //
        //    Output, int &LESS, &EQUAL, &MORE, the indexes in INDX of the
        //    entries of X that are just less than, equal to, and just greater
        //    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
        //    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
        //    is the greatest entry of X, then MORE is N+1.
        //
    {
        switch (n)
        {
            case <= 0:
                less = 0;
                equal = 0;
                more = 0;
                return;
        }

        int lo = 1;
        int hi = n;

        double xlo = x[indx[lo - 1] - 1];
        double ylo = y[indx[lo - 1] - 1];

        double xhi = x[indx[hi - 1] - 1];
        double yhi = y[indx[hi - 1] - 1];

        int compare = r8r8_compare(xval, yval, xlo, ylo);

        switch (compare)
        {
            case -1:
                less = 0;
                equal = 0;
                more = 1;
                return;
            case 0:
                less = 0;
                equal = 1;
                more = 2;
                return;
        }

        compare = r8r8_compare(xval, yval, xhi, yhi);

        switch (compare)
        {
            case 1:
                less = n;
                equal = 0;
                more = n + 1;
                return;
            case 0:
                less = n - 1;
                equal = n;
                more = n + 1;
                return;
        }

        for (;;)
        {
            if (lo + 1 == hi)
            {
                less = lo;
                equal = 0;
                more = hi;
                return;
            }

            int mid = (lo + hi) / 2;
            double xmid = x[indx[mid - 1] - 1];
            double ymid = y[indx[mid - 1] - 1];

            compare = r8r8_compare(xval, yval, xmid, ymid);

            switch (compare)
            {
                case 0:
                    equal = mid;
                    less = mid - 1;
                    more = mid + 1;
                    return;
                case -1:
                    hi = mid;
                    break;
                case +1:
                    lo = mid;
                    break;
            }
        }
    }

}