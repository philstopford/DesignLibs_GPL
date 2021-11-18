using System;
using System.Globalization;
using System.Numerics;
using Burkardt.SortNS;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void c8vec_copy(int n, Complex[] a1, ref Complex[] a2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_COPY copies a C8VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input, Complex A1[N], the vector to be copied.
        //
        //    Output, Complex A2[N], the copy of A1.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a2[i] = a1[i];
        }
    }

    public static Complex[] c8vec_copy_new(int n, Complex[] a1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_COPY_NEW copies a C8VEC to a "new" C8VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input, Complex A1[N], the vector to be copied.
        //
        //    Output, Complex C8VEC_COPY_NEW[N], the copy of A1.
        //
    {
        int i;

        Complex[] a2 = new Complex[n];

        for (i = 0; i < n; i++)
        {
            a2[i] = a1[i];
        }

        return a2;
    }

    public static Complex[] c8vec_indicator_new(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_INDICATOR_NEW sets a C8VEC to the indicator vector.
        //
        //  Discussion:
        //
        //    A C8VEC is a vector of Complex values.
        //
        //    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Output, Complex C8VEC_INDICATOR_NEW[N], the array.
        //
    {
        int i;

        Complex[] a = new Complex [n];

        for (i = 0; i < n; i++)
        {
            a[i] = new Complex(i + 1, -i - 1);
        }

        return a;
    }

    public static void c8vec_nint(int n, ref Complex[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_NINT rounds the entries of a C8VEC.
        //
        //  Discussion:
        //
        //    A C8VEC is a vector of Complex values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input/output, Complex A[N], the vector to be nint'ed.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i] = c8_nint(a[i]);
        }
    }

    public static double c8vec_norm_l1(int n, Complex[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_NORM_L1 returns the L1 norm of a C8VEC.
        //
        //  Discussion:
        //
        //    The vector L1 norm is defined as:
        //
        //      value = sum ( 1 <= I <= N ) abs ( A(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 February 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries.
        //
        //    Input, Complex A[N], the vector.
        //
        //    Output, double C8VEC_NORM_L1, the norm.
        //
    {
        int i;

        double value = 0.0;
        for (i = 0; i < n; i++)
        {
            value += c8_abs(a[i]);
        }

        return value;
    }

    public static double c8vec_norm_l2(int n, Complex[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_NORM_L2 returns the L2 norm of a C8VEC.
        //
        //  Discussion:
        //
        //    The vector L2 norm is defined as:
        //
        //      value = sqrt ( sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries.
        //
        //    Input, Complex A[N], the vector.
        //
        //    Output, double C8VEC_NORM_L2, the norm.
        //
    {
        int i;

        double value = 0.0;
        for (i = 0; i < n; i++)
        {
            value = value
                    + a[i].Real * a[i].Real
                    + a[i].Imaginary * a[i].Imaginary;
        }

        value = Math.Sqrt(value);

        return value;
    }

    public static double c8vec_norm_li(int n, Complex[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_NORM_LI returns the Loo norm of a C8VEC.
        //
        //  Discussion:
        //
        //    The vector Loo norm is defined as:
        //
        //      value = max ( 1 <= I <= N ) abs ( A(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 February 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries.
        //
        //    Input, Complex A[N], the vector.
        //
        //    Output, double C8VEC_NORM_L2, the norm.
        //
    {
        int i;

        double value = 0.0;
        for (i = 0; i < n; i++)
        {
            value = Math.Max(value, c8_abs(a[i]));
        }

        return value;
    }

    public static void c8vec_print(int n, Complex[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_PRINT prints a C8VEC.
        //
        //  Discussion:
        //
        //    A C8VEC is a vector of Complex values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input, Complex A[N], the vector to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine(title + "");

        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + ": " + a[i].Real
                                   + "  " + a[i].Imaginary + "");
        }
    }

    public static void c8vec_print_part(int n, Complex[] a, int max_print,
            string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_PRINT_PART prints "part" of a C8VEC.
        //
        //  Discussion:
        //
        //    The user specifies MAX_PRINT, the maximum number of lines to print.
        //
        //    If N, the size of the vector, is no more than MAX_PRINT, then
        //    the entire vector is printed, one entry per line.
        //
        //    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
        //    followed by a line of periods suggesting an omission,
        //    and the last entry.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries of the vector.
        //
        //    Input, Complex A[N], the vector to be printed.
        //
        //    Input, int MAX_PRINT, the maximum number of lines
        //    to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        switch (max_print)
        {
            case <= 0:
                return;
        }

        switch (n)
        {
            case <= 0:
                return;
        }

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");

        if (n <= max_print)
        {
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + a[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + a[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
        else
        {
            switch (max_print)
            {
                case >= 3:
                {
                    for (i = 0; i < max_print - 2; i++)
                    {
                        Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                               + "  " + a[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                               + "  " + a[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    }

                    Console.WriteLine("  ........  ..............  ..............");
                    i = n - 1;
                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + a[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + a[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
                }
                default:
                {
                    for (i = 0; i < max_print - 1; i++)
                    {
                        Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                               + "  " + a[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                               + "  " + a[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    }

                    i = max_print - 1;
                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + a[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + a[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + "...more entries...");
                    break;
                }
            }
        }
    }

    public static void c8vec_print_some(int n, Complex[] a, int i_lo, int i_hi,
            string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_PRINT_SOME prints some of a C8VEC.
        //
        //  Discussion:
        //
        //    A C8VEC is a vector of Complex values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input, Complex A[N], the vector to be printed.
        //
        //    Input, int I_LO, I_HI, the first and last entries to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");
        for (i = i_lo; i <= i_hi; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + ": " + a[i].Real
                                   + "  " + a[i].Imaginary + "");
        }

    }

    public static void c8vec_sort_a_l1(int n, ref Complex[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_SORT_A_L1 ascending sorts a C8VEC by L1 norm.
        //
        //  Discussion:
        //
        //    The L1 norm of A+Bi is sqrt ( A * A + B * B ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, length of input array.
        //
        //    Input/output, Complex X[N].
        //    On input, an unsorted array.
        //    On output, X has been sorted.
        //
    {
        int i = 0;
        int indx = 0;
        int isgn = 0;
        int j = 0;

        SortHeapExternalData data = new();
            
        for (;;)
        {
            Sort.sort_heap_external(ref data, n, ref indx, ref i, ref j, isgn);

            if (0 < indx)
            {
                (x[i - 1], x[j - 1]) = (x[j - 1], x[i - 1]);
            }
            else if (indx < 0)
            {
                if (c8_le_l1(x[i - 1], x[j - 1]))
                {
                    isgn = -1;
                }
                else
                {
                    isgn = +1;
                }
            }
            else
            {
                break;
            }
        }
    }

    public static void c8vec_sort_a_l2(int n, ref Complex[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_SORT_A_L2 ascending sorts a C8VEC by L2 norm.
        //
        //  Discussion:
        //
        //    The L2 norm of A+Bi is sqrt ( A * A + B * B ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, length of input array.
        //
        //    Input/output, Complex X[N].
        //    On input, an unsorted array.
        //    On output, X has been sorted.
        //
    {
        int i = 0;
        int indx = 0;
        int isgn = 0;
        int j = 0;
        SortHeapExternalData data = new();

        for (;;)
        {
            Sort.sort_heap_external(ref data, n, ref indx, ref i, ref j, isgn);

            if (0 < indx)
            {
                (x[i - 1], x[j - 1]) = (x[j - 1], x[i - 1]);
            }
            else if (indx < 0)
            {
                double normsq_i = Math.Pow(x[i - 1].Real, 2)
                                  + Math.Pow(x[i - 1].Imaginary, 2);

                double normsq_j = Math.Pow(x[j - 1].Real, 2)
                                  + Math.Pow(x[j - 1].Imaginary, 2);

                if (normsq_i < normsq_j)
                {
                    isgn = -1;
                }
                else
                {
                    isgn = +1;
                }
            }
            else
            {
                break;
            }
        }

    }

    public static void c8vec_sort_a_li(int n, ref Complex[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_SORT_A_LI ascending sorts a C8VEC by Loo norm.
        //
        //  Discussion:
        //
        //    The Loo norm of A+Bi is max ( abs ( A ), abs ( B ) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, length of input array.
        //
        //    Input/output, Complex X[N].
        //    On input, an unsorted array.
        //    On output, X has been sorted.
        //
    {
        int i = 0;
        int indx = 0;
        int isgn = 0;
        int j = 0;

        SortHeapExternalData data = new();
            
        for (;;)
        {
            Sort.sort_heap_external(ref data, n, ref indx, ref i, ref j, isgn);

            if (0 < indx)
            {
                (x[i - 1], x[j - 1]) = (x[j - 1], x[i - 1]);
            }
            else if (indx < 0)
            {
                if (c8_le_li(x[i - 1], x[j - 1]))
                {
                    isgn = -1;
                }
                else
                {
                    isgn = +1;
                }
            }
            else
            {
                break;
            }
        }
    }

    public static Complex[] c8vec_spiral_new(int n, int m, Complex c1,
            Complex c2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_SPIRAL_NEW returns N points on a spiral between C1 and C2.
        //
        //  Discussion:
        //
        //    A C8VEC is a vector of C8's.
        //
        //    Let the polar form of C1 be ( R1, T1 ) and the polar form of C2 
        //    be ( R2, T2 ) where, if necessary, we increase T2 by 2*PI so that T1 <= T2.
        //    
        //    Then the polar form of the I-th point C(I) is:
        //
        //      R(I) = ( ( N - I     ) * R1 
        //             + (     I - 1 ) * R2 ) 
        //              / ( N    - 1 )
        //
        //      T(I) = ( ( N - I     ) * T1 
        //             + (     I - 1 ) * ( T2 + M * 2 * PI ) ) 
        //             / ( N     - 1 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points on the spiral.
        //
        //    Input, int M, the number of full circuits the 
        //    spiral makes.
        //
        //    Input, Complex C1, C2, the first and last points 
        //    on the spiral.
        //
        //    Output, Complex C8VEC_SPIRAL_NEW[N], the points.
        //
    {
        int i;

        Complex[] c = new Complex[n];

        double r1 = c8_abs(c1);
        double r2 = c8_abs(c2);

        double t1 = c8_arg(c1);
        double t2 = c8_arg(c2);

        switch (m)
        {
            case 0:
            {
                if (t2 < t1)
                {
                    t2 += 2.0 * Math.PI;
                }

                break;
            }
            case > 0:
            {
                if (t2 < t1)
                {
                    t2 += 2.0 * Math.PI;
                }

                t2 += m * 2.0 * Math.PI;
                break;
            }
            case < 0:
            {
                if (t1 < t2)
                {
                    t2 -= 2.0 * Math.PI;
                }

                t2 -= m * 2.0 * Math.PI;
                break;
            }
        }

        for (i = 0; i < n; i++)
        {
            double ri = ((n - i - 1) * r1
                         + i * r2)
                        / (n - 1);

            double ti = ((n - i - 1) * t1
                         + i * t2)
                        / (n - 1);

            c[i] = polar_to_c8(ri, ti);
        }

        return c;
    }

    public static Complex[] c8vec_uniform_01_new(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_UNIFORM_01_NEW returns a unit pseudorandom C8VEC.
        //
        //  Discussion:
        //
        //    A C8VEC is a vector of Complex values.
        //
        //    The angles should be uniformly distributed between 0 and 2 * PI,
        //    the square roots of the radius uniformly distributed between 0 and 1.
        //
        //    This results in a uniform distribution of values in the unit circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of values to compute.
        //
        //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
        //    On output, SEED has been updated.
        //
        //    Output, Complex C8VEC_UNIFORM_01_NEW[N], the pseudorandom 
        //    complex vector.
        //
    {
        int i;

        Complex[] c = new Complex [n];

        for (i = 0; i < n; i++)
        {
            int k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            switch (seed)
            {
                case < 0:
                    seed += 2147483647;
                    break;
            }

            double r = Math.Sqrt(seed * 4.656612875E-10);

            k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            switch (seed)
            {
                case < 0:
                    seed += 2147483647;
                    break;
            }

            double theta = 2.0 * Math.PI * (seed * 4.656612875E-10);

            c[i] = r * new Complex(Math.Cos(theta), Math.Sin(theta));
        }

        return c;
    }

    public static Complex[] c8vec_unity_new(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_UNITY_NEW returns the N roots of unity in a C8VEC.
        //
        //  Discussion:
        //
        //    A C8VEC is a vector of Complex values.
        //
        //    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
        //
        //    X(1:N)^N = ( (1,0), (1,0), ..., (1,0) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Output, Complex C8VEC_UNITY_NEW[N], the N roots of unity.
        //
    {
        int i;

        Complex[] a = new Complex [n];

        for (i = 0; i < n; i++)
        {
            double theta = Math.PI * (2 * i) / n;
            a[i] = new Complex(Math.Cos(theta), Math.Sin(theta));
        }

        return a;
    }
}