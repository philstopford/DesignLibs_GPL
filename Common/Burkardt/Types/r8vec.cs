using System;
using System.IO;
using System.Linq;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_reverse ( int n, ref double[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_REVERSE reverses the elements of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of double precision values.
            //
            //    Input:
            //
            //      N = 5,
            //      A = ( 11.0, 12.0, 13.0, 14.0, 15.0 )
            //
            //    Output:
            //
            //      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 April 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input/output, double A[N], the array to be reversed.
            //
        {
            int i;
            double temp;

            for ( i = 0; i < n / 2; i++ )
            {
                temp     = a[i];
                a[i]     = a[n-1-i];
                a[n-1-i] = temp;
            }
        }

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
                a[i] = ( double ) ( i + 1 );
            }
        }
        
        public static bool r8vec_distinct ( int n, double[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2003
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
            //    Output, bool R8VEC_DISTINCT is true if all N elements of X are distinct.
            //
        {
            int i;
            int j;

            for ( i = 1; i <= n - 1; i++ )
            {
                for ( j = 1; j <= i - 1; j++ )
                {
                    if ( x[i] == x[j] )
                    {
                        return false;
                    }
                }
            }

            return true;
        }
        
        public static void r8vec_01_to_ab(int n, ref double[] a, double amax, double amin)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_01_TO_AB shifts and rescales data to lie within given bounds.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    On input, A contains the original data, which is presumed to lie
            //    between 0 and 1.  However, it is not necessary that this be so.
            //
            //    On output, A has been shifted and rescaled so that all entries which
            //    on input lay in [0,1] now lie between AMIN and AMAX.  Other entries will
            //    be mapped in a corresponding way.
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
            //    Input, int N, the number of data values.
            //
            //    Input/output, double A[N], the vector to be rescaled.
            //
            //    Input, double AMAX, AMIN, the maximum and minimum values
            //    allowed for A.
            //
        {
            double amax2;
            double amax3;
            double amin2;
            double amin3;
            int i;

            if (amax == amin)
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = amin;
                }

                return;
            }

            amax2 = Math.Max(amax, amin);
            amin2 = Math.Min(amax, amin);

            amin3 = r8vec_min(n, a);
            amax3 = r8vec_max(n, a);

            if (amax3 != amin3)
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = ((amax3 - a[i]) * amin2
                            + (a[i] - amin3) * amax2)
                           / (amax3 - amin3);
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = 0.5 * (amax2 + amin2);
                }
            }
        }

        public static void r8vec_add(int n, double[] a1, ref double[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_ADD adds one R8VEC to another.
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
            //    22 September 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, double A1[N], the vector to be added.
            //
            //    Input/output, double A2[N], the vector to be increased.
            //    On output, A2 = A2 + A1.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                a2[i] = a2[i] + a1[i];
            }
        }

        public static int r8vec_compare(int n, double[] a, double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_COMPARE compares two R8VEC's.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    The lexicographic ordering is used.
            //
            //  Example:
            //
            //    Input:
            //
            //      A1 = ( 2.0, 6.0, 2.0 )
            //      A2 = ( 2.0, 8.0, 12.0 )
            //
            //    Output:
            //
            //      ISGN = -1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, double A[N], B[N], the vectors to be compared.
            //
            //    Output, int R8VEC_COMPARE, the results of the comparison:
            //    -1, A is lexicographically less than B,
            //     0, A is equal to B,
            //    +1, A is lexicographically greater than B.
            //
        {
            int isgn;
            int k;

            isgn = 0;

            for (k = 0; k < n; k++)
            {
                if (a[k] < b[k])
                {
                    isgn = -1;
                    return isgn;
                }
                else if (b[k] < a[k])
                {
                    isgn = +1;
                    return isgn;
                }
            }

            return isgn;
        }








        public static double r8vec_length(int dim_num, double[] x)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_LENGTH returns the Euclidean length of an R8VEC
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double X[DIM_NUM], the vector.
            //
            //    Output, double R8VEC_LENGTH, the Euclidean length of the vector.
            //
        {
            double value = 0.0;
            for (int i = 0; i < dim_num; i++)
            {
                value = value + Math.Pow(x[i], 2);
            }

            value = Math.Sqrt(value);

            return value;
        }

        public static double r8vec_distance(int dim_num, double[] v1, double[] v2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_DISTANCE returns the Euclidean distance between two R8VEC's.
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
            //    11 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double V1[DIM_NUM], V2[DIM_NUM], the vectors.
            //
            //    Output, double R8VEC_DISTANCE, the Euclidean distance
            //    between the vectors.
            //
        {
            int i;
            double value;

            value = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                value = Math.Pow(v1[i] - v2[i], 2);
            }

            value = Math.Sqrt(value);

            return value;
        }

        public static void r8vec_divide(int n, ref double[] a, double s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_DIVIDE divides an R8VEC by a nonzero scalar.
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
            //    30 August 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input/output, double A[N].  On input, the vector to be scaled.
            //    On output, each entry has been divided by S.
            //
            //    Input, double S, the divisor.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                a[i] = a[i] / s;
            }

            return;
        }

        public static double r8vec_entropy(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_ENTROPY computes the entropy of an R8VEC.
            //
            //  Discussion:
            //
            //    Typically, the entries represent probabilities, and must sum to 1.
            //    For this function, the only requirement is that the entries be nonnegative.
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 August 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries.
            //
            //    Input, double X[N], the vector.
            //    Each entry must be nonnegative.
            //
            //    Output, double R8VEC_ENTROPY, the entropy of the
            //    normalized vector.
            //
        {
            int i;
            double value;
            double x_sum;
            double xi;

            for (i = 0; i < n; i++)
            {
                if (x[i] < 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8VEC_ENTROPY - Fatal error!");
                    Console.WriteLine("  Some entries are negative.");
                    return (1);
                }
            }

            x_sum = 0.0;
            for (i = 0; i < n; i++)
            {
                x_sum = x_sum + x[i];
            }

            if (x_sum == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_ENTROPY - Fatal error!");
                Console.WriteLine("  Entries sum to 0.");
                return (1);
            }

            value = 0.0;
            for (i = 0; i < n; i++)
            {
                if (0.0 < x[i])
                {
                    xi = x[i] / x_sum;
                    value = value - r8_log_2(xi) * xi;
                }
            }

            return value;
        }

        public static bool r8vec_is_integer ( int n, double[] a )

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
            bool value;

            value = true;

            for ( i = 0; i < n; i++ )
            {
                if ( a[i] != ( double ) ( int ) a[i] )
                {
                    value = false;
                    break;
                }
            }
            return value;
        }

        public static double[] r8vec_linspace_new(int n, double a_first, double a_last)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
            //
            //    In other words, the interval is divided into N-1 even subintervals,
            //    and the endpoints of intervals are used as the points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 March 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double A_FIRST, A_LAST, the first and last entries.
            //
            //    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
            //
        {
            double[] a = new double[n];

            if (n == 1)
            {
                a[0] = (a_first + a_last) / 2.0;
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    a[i] = ((double) (n - 1 - i) * a_first
                            + (double) (i) * a_last)
                           / (double) (n - 1);
                }
            }

            return a;
        }




        public static void r8vec_mesh_2d(int nx, int ny, double[] xvec, double[] yvec,
                ref double[] xmat, ref double[] ymat)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_MESH_2D creates a 2D mesh from X and Y vectors.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    NX = 2
            //    XVEC = ( 1, 2, 3 )
            //    NY = 3
            //    YVEC = ( 4, 5 )
            //
            //    XMAT = (
            //      1, 2, 3
            //      1, 2, 3 )
            //
            //    YMAT = (
            //      4, 4, 4
            //      5, 5, 5 ) 
            // 
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 July 2013
            //
            //  Parameters:
            //
            //    Input, int NX, NY, the number of X and Y values.
            //
            //    Input, double XVEC[NX], YVEC[NY], the X and Y coordinate
            //    values.
            //
            //    Output, double XMAT[NX*NY], YMAT[NX*NY], the coordinate
            //    values of points on an NX by NY mesh.
            //
        {
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    xmat[i + j * nx] = xvec[i];
                }
            }

            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    ymat[i + j * nx] = yvec[j];
                }
            }
        }

        public static void r8vec_permute(int n, int[] p, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_PERMUTE permutes an R8VEC in place.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    This routine permutes an array of real "objects", but the same
            //    logic can be used to permute an array of objects of any arithmetic
            //    type, or an array of objects of any complexity.  The only temporary
            //    storage required is enough to store a single object.  The number
            //    of data movements made is N + the number of cycles of order 2 or more,
            //    which is never more than N + N/2.
            //
            //  Example:
            //
            //    Input:
            //
            //      N = 5
            //      P = (   1,   3,   4,   0,   2 )
            //      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
            //
            //    Output:
            //
            //      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects.
            //
            //    Input, int P[N], the permutation.
            //
            //    Input/output, double A[N], the array to be permuted.
            //
        {
            double a_temp;
            int i;
            int iget;
            int iput;
            int istart;

            perm_check0(n, p);
            //
            //  In order for the sign negation trick to work, we need to assume that the
            //  entries of P are strictly positive.  Presumably, the lowest number is 0.
            //  So temporarily add 1 to each entry to force positivity.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] + 1;
            }

            //
            //  Search for the next element of the permutation that has not been used.
            //
            for (istart = 1; istart <= n; istart++)
            {
                if (p[istart - 1] < 0)
                {
                    continue;
                }
                else if (p[istart - 1] == istart)
                {
                    p[istart - 1] = -p[istart - 1];
                    continue;
                }
                else
                {
                    a_temp = a[istart - 1];
                    iget = istart;
                    //
                    //  Copy the new value into the vacated entry.
                    //
                    for (;;)
                    {
                        iput = iget;
                        iget = p[iget - 1];

                        p[iput - 1] = -p[iput - 1];

                        if (iget < 1 || n < iget)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("R8VEC_PERMUTE - Fatal error!");
                            Console.WriteLine("  A permutation index is out of range.");
                            Console.WriteLine("  P(" + iput + ") = " + iget + "");
                            return;
                        }

                        if (iget == istart)
                        {
                            a[iput - 1] = a_temp;
                            break;
                        }

                        a[iput - 1] = a[iget - 1];
                    }
                }
            }

            //
            //  Restore the signs of the entries.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = -p[i];
            }

            //
            //  Restore the entries.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] - 1;
            }

            return;
        }
        
        public static void r8vec_transpose_print ( int n, double[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Example:
        //
        //    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
        //    TITLE = 'My vector:  '
        //
        //    My vector:
        //        1.0    2.1    3.2    4.3    5.4
        //        6.5    7.6    8.7    9.8   10.9
        //       11.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input, double A[N], the vector to be printed.
        //
        //    Input, string TITLE, a title.
        //
        {
            int i;
            int ihi;
            int ilo;

            Console.WriteLine("");
            Console.WriteLine(title + "");

            if ( n <= 0 )
            {
                Console.WriteLine("  (Empty)");
                return;
            }

            for ( ilo = 0; ilo < n; ilo = ilo + 5 )
            {
                ihi = Math.Min ( ilo + 5, n );
                string cout = "";
                for ( i = ilo; i < ihi; i++ )
                {
                    cout += "  " + a[i].ToString().PadLeft(12);
                }
                Console.WriteLine(cout);
            }
        }

        public static double r8vec_product(int n, double[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_PRODUCT returns the product of the entries of an R8VEC.
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
            //    17 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double A[N], the vector.
            //
            //    Output, double R8VEC_PRODUCT, the product of the vector.
            //
        {
            int i;
            double product;

            product = 1.0;
            for (i = 0; i < n; i++)
            {
                product = product * a[i];
            }

            return product;
        }
        
        public static void r8vec_sort_heap_a ( int n, ref double[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SORT_HEAP_A ascending sorts an R8VEC using heap sort.
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
            //    30 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input/output, double A[N].
            //    On input, the array to be sorted;
            //    On output, the array has been sorted.
            //
        {
            int n1;
            double temp;

            if ( n <= 1 )
            {
                return;
            }
            //
            //  1: Put A into descending heap form.
            //
            r8vec_heap_d ( n, ref a );
            //
            //  2: Sort A.
            //
            //  The largest object in the heap is in A[0].
            //  Move it to position A[N-1].
            //
            temp = a[0];
            a[0] = a[n-1];
            a[n-1] = temp;
            //
            //  Consider the diminished heap of size N1.
            //
            for ( n1 = n-1; 2 <= n1; n1-- )
            {
                //
                //  Restore the heap structure of the initial N1 entries of A.
                //
                r8vec_heap_d ( n1, ref a );
                //
                //  Take the largest object from A[0] and move it to A[N1-1].
                //
                temp = a[0];
                a[0] = a[n1-1];
                a[n1-1] = temp;
            }
        }

        public static double[] r8vec_midspace_new(int n, double a, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_MIDSPACE_NEW creates a vector of linearly spaced values.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    This function divides the interval [a,b] into n subintervals, and then
            //    returns the midpoints of those subintervals.
            //
            //  Example:
            //
            //    N = 5, A = 10, B = 20
            //    X = [ 11, 13, 15, 17, 19 ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double A, B, the endpoints of the interval.
            //
            //    Output, double R8VEC_MIDSPACE_NEW[N], a vector of linearly spaced data.
            //
        {
            double[] x = new double[n];

            for (int i = 0; i < n; i++)
            {
                x[i] = ((double) (2 * n - 2 * i - 1) * a
                        + (double) (2 * i + 1) * b)
                       / (double) (2 * n);
            }

            return x;
        }



        public static bool r8vec_eq(int n, double[] a1, int startIndexA1, double[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_EQ is true if every pair of entries in two vectors is equal.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, double A1[N], A2[N], two vectors to compare.
            //
            //    Output, bool R8VEC_EQ.
            //    R8VEC_EQ is TRUE if every pair of elements A1(I) and A2(I) are equal,
            //    and FALSE otherwise.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                if (a1[startIndexA1 + i] != a2[i])
                {
                    return false;
                }
            }

            return true;
        }

        public static bool r8vec_gt(int n, double[] a1, int startIndexA1, double[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_GT == ( A1 > A2 ) for real vectors.
            //
            //  Discussion:
            //
            //    The comparison is lexicographic.
            //
            //    A1 > A2  <=>                              A1(1) > A2(1) or
            //                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
            //                 ...
            //                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the vectors.
            //
            //    Input, double A1[N], A2[N], the vectors to be compared.
            //
            //    Output, bool R8VEC_GT, is TRUE if and only if A1 > A2.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                if (a2[i] < a1[startIndexA1 + i])
                {
                    return true;
                }
                else if (a1[startIndexA1 + i] < a2[i])
                {
                    return false;
                }
            }

            return false;
        }

        public static bool r8vec_lt(int n, double[] a1, int startIndexA1, double[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_LT == ( A1 < A2 ) for real vectors.
            //
            //  Discussion:
            //
            //    The comparison is lexicographic.
            //
            //    A1 < A2  <=>                              A1(1) < A2(1) or
            //                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
            //                 ...
            //                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the vectors.
            //
            //    Input, double A1[N], A2[N], the vectors to be compared.
            //
            //    Output, bool R8VEC_LT, is TRUE if and only if A1 < A2.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                if (a1[startIndexA1 + i] < a2[i])
                {
                    return true;
                }
                else if (a2[i] < a1[startIndexA1 + i])
                {
                    return false;
                }

            }

            return false;
        }

        public static void r8vec_part_quick_a(int n, ref double[] a, int startIndexA, ref int l, ref int r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_PART_QUICK_A reorders an R8VEC as part of a quick sort.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    The routine reorders the entries of A.  Using A[0] as a
            //    key, all entries of A that are less than or equal to A[0] will
            //    precede A[0] which precedes all entries that are greater than A[0].
            //
            //  Example:
            //
            //    Input:
            //
            //  N = 8
            //
            //  A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
            //
            //    Output:
            //
            //  L = 3, R = 6
            //
            //  A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
            //        -------        -------
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries of A.
            //
            //    Input/output, double A[N].  On input, the array to be checked.
            //    On output, A has been reordered as described above.
            //
            //    Output, int L, R, the indices of A that define the three segments.
            //    Let KEY = the input value of A[0].  Then
            //    I <= L             A(I) < KEY;
            //     L < I < R         A(I) = KEY;
            //             R <= I    A(I) > KEY.
            //
        {
            int i;
            double key;
            int m;
            double temp;

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_PART_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            }
            else if (n == 1)
            {
                l = 0;
                r = 2;
                return;
            }

            key = a[startIndexA + 0];
            m = 1;
            //
            //  The elements of unknown size have indices between L+1 and R-1.
            //
            l = 1;
            r = n + 1;

            for (i = 2; i <= n; i++)
            {

                if (key < a[startIndexA + l])
                {
                    r = r - 1;
                    temp = a[startIndexA + (r - 1)];
                    a[startIndexA + (r - 1)] = a[startIndexA + l];
                    a[startIndexA + l] = temp;
                }
                else if (a[startIndexA + (l)] == key)
                {
                    m = m + 1;
                    temp = a[startIndexA + (m - 1)];
                    a[startIndexA + (m - 1)] = a[startIndexA + (l)];
                    a[startIndexA + (l)] = temp;
                    l = l + 1;
                }
                else if (a[startIndexA + (l)] < key)
                {
                    l = l + 1;
                }

            }

            //
            //  Now shift small elements to the left, and KEY elements to center.
            //
            for (i = 1; i <= l - m; i++)
            {
                a[startIndexA + (i - 1)] = a[startIndexA + (i + m - 1)];
            }

            l = l - m;

            for (i = l + 1; i <= l + m; i++)
            {
                a[startIndexA + (i - 1)] = key;
            }
        }


        public static void r8vec_swap(int n, ref double[] a1, int startIndexA1, ref double[] a2, int startIndexA2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SWAP swaps the entries of two R8VEC's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the arrays.
            //
            //    Input/output, double A1[N], A2[N], the vectors to swap.
            //
        {
            for (int i = 0; i < n; i++)
            {
                double temp = a1[i + startIndexA1];
                a1[i + startIndexA1] = a2[i + startIndexA2];
                a2[i + startIndexA2] = temp;
            }
        }

        public static int[] r8vec_first_index(int n, double[] a, double tol )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_FIRST_INDEX indexes the first occurrence of values in an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    For element A(I) of the vector, FIRST_INDEX(I) is the index in A of
        //    the first occurrence of the value A(I).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Input, double A[N], the unsorted array to examine.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int R8VEC_FIRST_INDEX[N], the first occurrence index.
        //
        {
            int[] first_index;
            int i;
            int j;

            first_index = new int[n];

            for (i = 0; i < n; i++)
            {
                first_index[i] = -1;
            }

            for (i = 0; i < n; i++)
            {
                if (first_index[i] == -1)
                {
                    first_index[i] = i;
                    for (j = i + 1; j < n; j++)
                    {
                        if (Math.Abs(a[i] - a[j]) <= tol)
                        {
                            first_index[j] = i;
                        }
                    }
                }
            }

            return first_index;
        }

        public static double r8vec_frac(int n, ref double[] a, int k )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_FRAC searches for the K-th smallest entry in an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    Hoare's algorithm is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 August 2004
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Input/output, double A[N].
        //    On input, A is the array to search.
        //    On output, the elements of A have been somewhat rearranged.
        //
        //    Input, int K, the fractile to be sought.  If K = 1, the minimum
        //    entry is sought.  If K = N, the maximum is sought.  Other values
        //    of K search for the entry which is K-th in size.  K must be at
        //    least 1, and no greater than N.
        //
        //    Output, double R8VEC_FRAC, the value of the K-th fractile of A.
        //
        {
            double frac;
            int i;
            int iryt;
            int j;
            int left;
            double temp;
            double x;

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_FRAC - Fatal error!");
                Console.WriteLine("  Illegal nonpositive value of N = " + n + "");
                return(1);
            }

            if (k <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_FRAC - Fatal error!");
                Console.WriteLine("  Illegal nonpositive value of K = " + k + "");
                return(1);
            }

            if (n < k)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_FRAC - Fatal error!");
                Console.WriteLine("  Illegal N < K, K = " + k + "");
                return(1);
            }

            left = 1;
            iryt = n;

            for (;;)
            {
                if (iryt <= left)
                {
                    frac = a[k - 1];
                    break;
                }

                x = a[k - 1];
                i = left;
                j = iryt;

                for (;;)
                {
                    if (j < i)
                    {
                        if (j < k)
                        {
                            left = i;
                        }

                        if (k < i)
                        {
                            iryt = j;
                        }

                        break;
                    }

                    //
                    //  Find I so that X <= A(I).
                    //
                    while (a[i - 1] < x)
                    {
                        i = i + 1;
                    }

                    //
                    //  Find J so that A(J) <= X.
                    //
                    while (x < a[j - 1])
                    {
                        j = j - 1;
                    }

                    if (i <= j)
                    {
                        temp = a[i - 1];
                        a[i - 1] = a[j - 1];
                        a[j - 1] = temp;
                        i = i + 1;
                        j = j - 1;
                    }
                }
            }

            return frac;
        }

        public static double[] r8vec_fraction(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_FRACTION returns the fraction parts of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    If we regard a real number as
            //
            //      R8 = SIGN * ( WHOLE + FRACTION )
            //
            //    where
            //
            //      SIGN is +1 or -1,
            //      WHOLE is a nonnegative integer
            //      FRACTION is a nonnegative real number strictly less than 1,
            //
            //    then this routine returns the value of FRACTION.
            //
            //  Example:
            //
            //     R8    R8_FRACTION
            //
            //    0.00      0.00
            //    1.01      0.01
            //    2.02      0.02
            //   19.73      0.73
            //   -4.34      0.34
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of arguments.
            //
            //    Input, double X[N], the arguments.
            //
            //    Output, double R8_FRACTION[N], the fraction parts.
            //
        {
            double[] fraction;
            int i;

            fraction = new double[n];

            for (i = 0; i < n; i++)
            {
                fraction[i] = Math.Abs(x[i]) - (double) ((int) (Math.Abs(x[i])));
            }

            return fraction;
        }

        public static bool r8vec_gt(int n, double[] a1, double[] a2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_GT == ( A1 > A2 ) for two R8VEC's.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The comparison is lexicographic.
        //
        //    A1 > A2  <=>                              A1(1) > A2(1) or
        //                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
        //                 ...
        //                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vectors.
        //
        //    Input, double A1[N], A2[N], the vectors to be compared.
        //
        //    Output, bool R8VEC_GT, is TRUE if and only if A1 > A2.
        //
        {
            int i;

            for (i = 0; i < n; i++)
            {

                if (a2[i] < a1[i])
                {
                    return true;
                }
                else if (a1[i] < a2[i])
                {
                    return false;
                }

            }

            return false;
        }

        public static double[] r8vec_house_column ( int n, double[] a_vec, int k )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The routine returns a vector V that defines a Householder
        //    premultiplier matrix H(V) that zeros out the subdiagonal entries of
        //    column K of the matrix A.
        //
        //       H(V) = I - 2 * v * v'
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Input, double A_VEC[N], a row or column of the matrix A.
        //
        //    Input, int K, the index of the row or column.
        //
        //    Output, double R8VEC_HOUSE_COLUMN[N], a vector of unit L2 norm which
        //    defines an orthogonal Householder premultiplier matrix H with the property
        //    that the K-th column of H*A is zero below the diagonal.
        //
        {
            int i;
            double s;
            double[] v;

            v = r8vec_zeros_new ( n );

            if ( k < 1 || n <= k )
            {
                return v;
            }

            s = r8vec_norm_l2 ( n+1-k, a_vec, +k-1 );

            if ( s == 0.0 )
            {
                return v;
            }

            v[k-1] = a_vec[k-1] + Math.Abs ( s ) * r8_sign ( a_vec[k-1] );

            r8vec_copy ( n-k, a_vec, ref v, k, k );
            //
            //  Normalize.
            //
            s = r8vec_norm_l2 ( n-k+1, v, +k-1 );

            for ( i = k - 1; i < n; i++ )
            {
                v[i] = v[i] / s;
            }

            return v;
        }

        public static double[] r8vec_identity_row_new ( int n, int i )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_IDENTITY_ROW_NEW sets an R8VEC to the I-th row of the identity.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 March 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Input, int I, indicates the row.  0 <= I < N.
            //
            //    Output, double R8VEC_IDENTITY_ROW_NEW[N], the array.
            //
        {
            double[] a;
            int j;

            a = new double[n];

            for ( j = 0; j < n; j++ )
            {
                a[j] = 0.0;
            }

            if ( 0 <= i && i < n )
            {
                a[i] = 1.0;
            }

            return a;
        }
        
        public static double[] r8vec_ones_new ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_ONES_NEW creates a vector of 1's.
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
            //    14 March 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Output, double R8VEC_ONES_NEW[N], a vector of 1's.
            //
        {
            double[] a;
            int i;

            a = new double[n];

            for ( i = 0; i < n; i++ )
            {
                a[i] = 1.0;
            }
            return a;
        }

        public static bool r8vec_in_ab ( int n, double[] x, double a, double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_IN_AB is TRUE if the entries of an R8VEC are in the range [A,B].
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
        //    Output, bool R8VEC_IN_AB, is TRUE if every entry is
        //    between A and B.
        //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                if ( x[i] < a || b < x[i] )
                {
                    return false;
                }
            }
            return true;
        }
        
        public static void r8vec_backtrack ( int n, int maxstack, double[] stack, ref double[] x, 
        ref int indx, ref int k, ref int nstack, ref int[] ncan )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BACKTRACK supervises a backtrack search for a real vector.
        //
        //  Discussion:
        //
        //    The routine tries to construct a real vector one index at a time,
        //    using possible candidates as supplied by the user.
        //
        //    At any time, the partially constructed vector may be discovered to be
        //    unsatisfactory, but the routine records information about where the
        //    last arbitrary choice was made, so that the search can be
        //    carried out efficiently, rather than starting out all over again.
        //
        //    First, call the routine with INDX = 0 so it can initialize itself.
        //
        //    Now, on each return from the routine, if INDX is:
        //      1, you've just been handed a complete candidate vector;
        //         Admire it, analyze it, do what you like.
        //      2, please determine suitable candidates for position X(K).
        //         Return the number of candidates in NCAN(K), adding each
        //         candidate to the end of STACK, and increasing NSTACK.
        //      3, you're done.  Stop calling the routine;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 July 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms for Computers and Calculators,
        //    Second Edition,
        //    Academic Press, 1978,
        //    ISBN: 0-12-519260-6,
        //    LC: QA164.N54.
        //
        //  Parameters:
        //
        //    Input, int N, the number of positions to be filled in the vector.
        //
        //    Input, int MAXSTACK, the maximum length of the stack.
        //
        //    Input, double STACK[MAXSTACK], a list of all current candidates for
        //    all positions 1 through K.
        //
        //    Input/output, double X[N], the partial or complete candidate vector.
        //
        //    Input/output, int &INDX, a communication flag.
        //    On input,
        //      0 to start a search.
        //    On output:
        //      1, a complete output vector has been determined and returned in X(1:N);
        //      2, candidates are needed for position X(K);
        //      3, no more possible vectors exist.
        //
        //    Inout/output, int &K, if INDX=2, the current vector index being considered.
        //
        //    Input/output, int &NSTACK, the current length of the stack.
        //
        //    Input/output, int NCAN[N], lists the current number of candidates for
        //    positions 1 through K.
        //
        {
            //
            //  If this is the first call, request a candidate for position 1.
            //
            if (indx == 0)
            {
                k = 1;
                nstack = 0;
                indx = 2;
                return;
            }

            //
            //  Examine the stack.
            //
            for (;;)
            {
                //
                //  If there are candidates for position K, take the first available
                //  one off the stack, and increment K.
                //
                //  This may cause K to reach the desired value of N, in which case
                //  we need to signal the user that a complete set of candidates
                //  is being returned.
                //
                if (0 < ncan[k - 1])
                {
                    x[k - 1] = stack[nstack - 1];
                    nstack = nstack - 1;

                    ncan[k - 1] = ncan[k - 1] - 1;

                    if (k != n)
                    {
                        k = k + 1;
                        indx = 2;
                    }
                    else
                    {
                        indx = 1;
                    }

                    break;
                }
                //
                //  If there are no candidates for position K, then decrement K.
                //  If K is still positive, repeat the examination of the stack.
                //
                else
                {
                    k = k - 1;

                    if (k <= 0)
                    {
                        indx = 3;
                        break;
                    }
                }
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
            double[] a;
            int i;

            a = new double[n];

            for ( i = 0; i <= n-1; i++ ) 
            {
                a[i] = ( double ) ( i + 1 );
            }

            return a;
        }

    }
}