using System;
using Burkardt.Probability;

namespace Burkardt.Types
{
    public class r8
    {
        public bool error { get; set; }
        public double val { get; set; }
        public int lchar { get; set; }
    }

    public class r8vec
    {
        public bool error { get; set; }
        public double[] rvec { get; set; }
    }

    public static partial class typeMethods
    {
        public static r8 s_to_r8(string s)
//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
        {
            r8 ret = new r8 {lchar = -1};
            double rexp;
            char TAB = (char) 9;

            int nchar = s_len_trim(s);
            int isgn = 1;
            double rtop = 0.0;
            double rbot = 1.0;
            int jsgn = 1;
            int jtop = 0;
            int jbot = 1;
            int ihave = 1;
            int iterm = 0;

            for (;;)
            {
                char c = s[ret.lchar + 1];
                ret.lchar = ret.lchar + 1;
//
//  Blank or TAB character.
//
                if (c == ' ' || c == TAB)
                {
                    if (ihave == 2)
                    {
                    }
                    else if (ihave == 6 || ihave == 7)
                    {
                        iterm = 1;
                    }
                    else if (1 < ihave)
                    {
                        ihave = 11;
                    }
                }
//
//  Comma.
//
                else if (c == ',' || c == ';')
                {
                    if (ihave != 1)
                    {
                        iterm = 1;
                        ihave = 12;
                        ret.lchar = ret.lchar + 1;
                    }
                }
//
//  Minus sign.
//
                else if (c == '-')
                {
                    if (ihave == 1)
                    {
                        ihave = 2;
                        isgn = -1;
                    }
                    else if (ihave == 6)
                    {
                        ihave = 7;
                        jsgn = -1;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
//
//  Plus sign.
//
                else if (c == '+')
                {
                    if (ihave == 1)
                    {
                        ihave = 2;
                    }
                    else if (ihave == 6)
                    {
                        ihave = 7;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
//
//  Decimal point.
//
                else if (c == '.')
                {
                    if (ihave < 4)
                    {
                        ihave = 4;
                    }
                    else if (6 <= ihave && ihave <= 8)
                    {
                        ihave = 9;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
//
//  Exponent marker.
//
                else if ((Char.ToUpper(c) == 'E') || (Char.ToUpper(c) == 'D'))
                {
                    if (ihave < 6)
                    {
                        ihave = 6;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
//
//  Digit.
//
                else if (ihave < 11 && '0' <= c && c <= '9')
                {
                    if (ihave <= 2)
                    {
                        ihave = 3;
                    }
                    else if (ihave == 4)
                    {
                        ihave = 5;
                    }
                    else if (ihave == 6 || ihave == 7)
                    {
                        ihave = 8;
                    }
                    else if (ihave == 9)
                    {
                        ihave = 10;
                    }

                    int ndig = ch_to_digit(c);

                    if (ihave == 3)
                    {
                        rtop = 10.0 * rtop + (double) ndig;
                    }
                    else if (ihave == 5)
                    {
                        rtop = 10.0 * rtop + (double) ndig;
                        rbot = 10.0 * rbot;
                    }
                    else if (ihave == 8)
                    {
                        jtop = 10 * jtop + ndig;
                    }
                    else if (ihave == 10)
                    {
                        jtop = 10 * jtop + ndig;
                        jbot = 10 * jbot;
                    }
                }
//
//  Anything else is regarded as a terminator.
//
                else
                {
                    iterm = 1;
                }

//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
                if (iterm == 1 || nchar <= ret.lchar + 1)
                {
                    break;
                }

            }

//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
            if (iterm != 1 && (ret.lchar) + 1 == nchar)
            {
                ret.lchar = nchar;
            }

//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
            if (ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7)
            {
                ret.error = true;
                return ret;
            }

//
//  Number seems OK.  Form it.
//
            if (jtop == 0)
            {
                rexp = 1.0;
            }
            else
            {
                if (jbot == 1)
                {
                    rexp = Math.Pow(10.0, jsgn * jtop);
                }
                else
                {
                    rexp = jsgn * jtop;
                    rexp = rexp / jbot;
                    rexp = Math.Pow(10.0, rexp);
                }
            }

            ret.val = isgn * rexp * rtop / rbot;

            return ret;
        }


        public static r8vec s_to_r8vec(string s, int n)
//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
        {
            int begin = 0;
            int length = s.Length;

            r8vec ret = new r8vec() {rvec = new double[n]};

            for (int i = 0; i < n; i++)
            {
                r8 res = s_to_r8(s.Substring(begin, length));

                ret.rvec[i] = res.val;
                int lchar = res.lchar;

                if (ret.error)
                {
                    return ret;
                }

                begin = begin + lchar;
                length = length - lchar;
            }

            return ret;
        }

        public static double r8vec_sum(int n, double[] a)
//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
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
//    15 October 2004
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
//    Output, double R8VEC_SUM, the sum of the vector.
//
        {
            double value = 0.0;
            for (int i = 0; i < n; i++)
            {
                value = value + a[i];
            }

            return value;
        }


        public static int[] r8vec_sort_heap_index_a_new(int n, double[] a)

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_INDEX_A_NEW does an indexed heap ascending sort of an R8VEC.
//
//  Discussion:
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      A(INDX(I)), I = 1 to N is sorted,
//
//    after which A(I), I = 1 to N is sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], an array to be index-sorted.
//
//    Output, int R8VEC_SORT_HEAP_INDEX_A_NEW[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
        {
            int[] indx = new int[n];

            for (int i = 1; i <= n; i++)
            {
                indx[i - 1] = i;
            }

            int l = n / 2 + 1;
            int ir = n;

            for (;;)
            {
                double aval;
                int indxt;
                if (1 < l)
                {
                    l = l - 1;
                    indxt = indx[l - 1];
                    aval = a[indxt - 1];
                }
                else
                {
                    indxt = indx[ir - 1];
                    aval = a[indxt - 1];
                    indx[ir - 1] = indx[0];
                    ir = ir - 1;

                    if (ir == 1)
                    {
                        indx[0] = indxt;
                        for (int i = 0; i < n; i++)
                        {
                            indx[i] = indx[i] - 1;
                        }

                        break;
                    }
                }

                int i2 = l;
                int j = l + l;

                while (j <= ir)
                {
                    if (j < ir)
                    {
                        if (a[indx[j - 1] - 1] < a[indx[j] - 1])
                        {
                            j = j + 1;
                        }
                    }

                    if (aval < a[indx[j - 1] - 1])
                    {
                        indx[i2 - 1] = indx[j - 1];
                        i2 = j;
                        j = j + j;
                    }
                    else
                    {
                        j = ir + 1;
                    }
                }

                indx[i2 - 1] = indxt;
            }

            return indx;
        }

        public static long r8_nint(double x)
//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Examples:
//
//        X         R8_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
        {
            long value = (long) (Math.Abs(x) + 0.5);

            if (x < 0.0)
            {
                value = -value;
            }

            return value;
        }

        public static void r8pp_print(int n, double[] a, string title)
//****************************************************************************80
//
//  Purpose:
//
//    R8PP_PRINT prints a R8PP matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
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
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the R8PP matrix.
//
//    Input, string TITLE, a title.
//
        {
            r8pp_print_some(n, a, 1, 1, n, n, title);
        }

        public static void r8pp_print_some(int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)
//****************************************************************************80
//
//  Purpose:
//
//    R8PP_PRINT_SOME prints some of a R8PP matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
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
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the R8PP matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
        {
            int INCX = 5;

            Console.WriteLine("");
            Console.WriteLine(title + "");
//
//  Print the columns of the matrix, in strips of 5.
//
            for (int j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX)
            {
                int j2hi = j2lo + INCX - 1;
                j2hi = Math.Min(j2hi, n);
                j2hi = Math.Min(j2hi, jhi);

                Console.WriteLine("");
                string cout = "  Col: ";
                int j;
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
                int i2lo = Math.Max(ilo, 1);
                int i2hi = Math.Min(ihi, n);

                for (int i = i2lo; i <= i2hi; i++)
                {
                    cout = i.ToString().PadLeft(6) + "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
                    for (j = j2lo; j <= j2hi; j++)
                    {
                        double aij;
                        if (i <= j)
                        {
                            aij = a[i - 1 + (j * (j - 1)) / 2];
                        }
                        else
                        {
                            aij = a[j - 1 + (i * (i - 1)) / 2];
                        }

                        cout += aij.ToString().PadLeft(12) + "  ";
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        public static void r8utp_print(int n, double[] a, string title)
//****************************************************************************80
//
//  Purpose:
//
//    R8UTP_PRINT prints a R8UTP matrix.
//
//  Discussion:
//
//    The R8UTP storage format is appropriate for an upper triangular
//    matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the matrix.
//
//    Input, string TITLE, a title.
//
        {
            r8utp_print_some(n, a, 1, 1, n, n, title);

            return;
        }

        public static void r8utp_print_some(int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)
//****************************************************************************80
//
//  Purpose:
//
//    R8UTP_PRINT_SOME prints some of an R8UTP matrix.
//
//  Discussion:
//
//    The R8UTP storage format is appropriate for an upper triangular
//    matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the matrix.
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

            Console.WriteLine("");
            Console.WriteLine(title + "");
//
//  Print the columns of the matrix, in strips of 5.
//
            for (j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX)
            {
                j2hi = j2lo + INCX - 1;
                j2hi = Math.Min(j2hi, n);
                j2hi = Math.Min(j2hi, jhi);

                Console.WriteLine("");
                string cout = "  Col: ";
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
                i2lo = Math.Max(ilo, 1);
                i2hi = Math.Min(ihi, n);

                for (i = i2lo; i <= i2hi; i++)
                {
                    cout = i.ToString().PadLeft(6) + "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
                    for (j = j2lo; j <= j2hi; j++)
                    {
                        if (i <= j)
                        {
                            aij = a[i - 1 + (j * (j - 1)) / 2];
                        }
                        else
                        {
                            aij = 0.0;
                        }

                        cout += aij.ToString().PadLeft(12) + "  ";
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        public static double r8_huge()
//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
        {
            double value;

            value = 1.0E+30;

            return value;
        }

        public static double r8poly_value(int n, double[] a, double x)
//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_VALUE evaluates a double precision polynomial.
//
//  Discussion:
//
//    For sanity's sake, the value of N indicates the NUMBER of 
//    coefficients, or more precisely, the ORDER of the polynomial,
//    rather than the DEGREE of the polynomial.  The two quantities
//    differ by 1, but cause a great deal of confusion.
//
//    Given N and A, the form of the polynomial is:
//
//      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double A[N], the coefficients of the polynomial.
//    A[0] is the constant term.
//
//    Input, double X, the point at which the polynomial is to be evaluated.
//
//    Output, double R8POLY_VALUE, the value of the polynomial at X.
//
        {
            double value = 0.0;

            for (int i = n - 1; 0 <= i; i--)
            {
                value = value * x + a[i];
            }

            return value;
        }

        public static double r8vec_dot_product(int n, double[] a1, double[] a2)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
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
        //    Input, double A1[N], A2[N], the two vectors to be considered.
        //
        //    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
        //
        {
            double value = 0.0;
            for (int i = 0; i < n; i++)
            {
                value = value + a1[i] * a2[i];
            }

            return value;
        }

        public static void r8vec_print(int n, double[] a, string title)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_PRINT prints an R8VEC.
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
        //    16 August 2004
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

            Console.WriteLine(title);
            Console.WriteLine();
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + ": " + a[i].ToString().PadLeft(14) + "");
            }
        }

        public static void r8vec_unit_sum ( int n, ref double[] a )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_UNIT_SUM normalizes an R8VEC to have unit sum.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input/output, double A[N], the vector to be normalized.
        //    On output, the entries of A should have unit sum.  However, if
        //    the input vector has zero sum, the routine halts.
        //
        {
            double a_sum = 0.0;
            for (int i = 0; i < n; i++ )
            {
                a_sum = a_sum + a[i];
            }

            if ( a_sum == 0.0 )
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_UNIT_SUM - Fatal error!");
                Console.WriteLine("  The vector entries sum to 0.");
            }

            for (int i = 0; i < n; i++ )
            {
                a[i] = a[i] / a_sum;
            }
        }

        public static double[] r8mat_mv_new(int m, int n, double[] a, double[] x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MV_NEW multiplies a matrix times a vector.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[M,N], the M by N matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8MAT_MV_NEW[M], the product A*X.
        //
        {
            double[] y = new double[m];

            for (int i = 0; i < m; i++)
            {
                y[i] = 0.0;
                for (int j = 0; j < n; j++)
                {
                    y[i] = y[i] + a[i + j * m] * x[j];
                }
            }

            return y;
        }

        public static double[] r8col_variance(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_VARIANCE returns the variances of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double A[M*N], the array whose variances are desired.
        //
        //    Output, double R8COL_VARIANCE[N], the variances of the rows.
        //
        {
            double[] variance = new double[n];

            for (int j = 0; j < n; j++)
            {
                double mean = 0.0;
                for (int i = 0; i < m; i++)
                {
                    mean = mean + a[i + j * m];
                }

                mean = mean / (double) (m);

                variance[j] = 0.0;
                for (int i = 0; i < m; i++)
                {
                    variance[j] = variance[j] + Math.Pow(a[i + j * m] - mean, 2);
                }

                if (1 < m)
                {
                    variance[j] = variance[j] / (double) (m - 1);
                }
                else
                {
                    variance[j] = 0.0;
                }
            }

            return variance;
        }

        public static double[] r8col_mean(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_MEAN returns the column means of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Example:
        //
        //    A =
        //      1  2  3
        //      2  6  7
        //
        //    R8COL_MEAN =
        //      1.5  4.0  5.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8COL_MEAN[N], the means, or averages, of the columns.
        //
        {
            double[] mean = new double[n];

            for (int j = 0; j < n; j++)
            {
                mean[j] = 0.0;
                for (int i = 0; i < m; i++)
                {
                    mean[j] = mean[j] + a[i + j * m];
                }

                mean[j] = mean[j] / (double) (m);
            }

            return mean;
        }

        public static double r8_zeta(double p)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ZETA estimates the Riemann Zeta function.
        //
        //  Discussion:
        //
        //    For 1 < P, the Riemann Zeta function is defined as:
        //
        //      ZETA ( P ) = Sum ( 1 <= N < oo ) 1 / N^P
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Daniel Zwillinger, editor,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996.
        //
        //  Parameters:
        //
        //    Input, double P, the power to which the integers are raised.
        //    P must be greater than 1.  For integral P up to 20, a
        //    precomputed value is returned; otherwise the infinite
        //    sum is approximated.
        //
        //    Output, double R8_ZETA, an approximation to the Riemann
        //    Zeta function.
        //
        {
            int n;
            const double r8_huge = 1.0E+30;
            const double r8_pi = 3.14159265358979323;
            double value;

            if (p <= 1.0)
            {
                value = r8_huge;
            }
            else if (p == 2.0)
            {
                value = Math.Pow(r8_pi, 2) / 6.0;
            }
            else if (p == 3.0)
            {
                value = 1.2020569032;
            }
            else if (p == 4.0)
            {
                value = Math.Pow(r8_pi, 4) / 90.0;
            }
            else if (p == 5.0)
            {
                value = 1.0369277551;
            }
            else if (p == 6.0)
            {
                value = Math.Pow(r8_pi, 6) / 945.0;
            }
            else if (p == 7.0)
            {
                value = 1.0083492774;
            }
            else if (p == 8.0)
            {
                value = Math.Pow(r8_pi, 8) / 9450.0;
            }
            else if (p == 9.0)
            {
                value = 1.0020083928;
            }
            else if (p == 10.0)
            {
                value = Math.Pow(r8_pi, 10) / 93555.0;
            }
            else if (p == 11.0)
            {
                value = 1.0004941886;
            }
            else if (p == 12.0)
            {
                value = 1.0002460866;
            }
            else if (p == 13.0)
            {
                value = 1.0001227133;
            }
            else if (p == 14.0)
            {
                value = 1.0000612482;
            }
            else if (p == 15.0)
            {
                value = 1.0000305882;
            }
            else if (p == 16.0)
            {
                value = 1.0000152823;
            }
            else if (p == 17.0)
            {
                value = 1.0000076372;
            }
            else if (p == 18.0)
            {
                value = 1.0000038173;
            }
            else if (p == 19.0)
            {
                value = 1.0000019082;
            }
            else if (p == 20.0)
            {
                value = 1.0000009540;
            }
            else
            {
                double zsum = 0.0;
                n = 0;

                for (;;)
                {
                    n = n + 1;
                    double zsum_old = zsum;
                    zsum = zsum + 1.0 / Math.Pow((double) n, p);
                    if (zsum <= zsum_old)
                    {
                        break;
                    }
                }

                value = zsum;
            }

            return value;
        }

        public static double r8_beta(double x, double y)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BETA returns the value of the Beta function.
        //
        //  Discussion:
        //
        //    BETA(X,Y) = ( GAMMA(X) * GAMMA(Y) ) / GAMMA(X+Y)
        //
        //    BETA(X,Y) = BETA(Y,X).
        //    BETA(X,Y) = Integral ( 0 <= T <= 1 ) T^(X-1) (1-T)^(Y-1) dT.
        //    BETA(X,Y) = GAMMA(X) * GAMMA(Y) / GAMMA(X+Y)
        //
        //    Both X and Y must be greater than 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the two parameters that define the Beta function.
        //    X and Y must be greater than 0.
        //
        //    Output, double R8_BETA, the value of the Beta function.
        //
        {
            if (x <= 0.0 || y <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_BETA - Fatal error!");
                Console.WriteLine("  Both X and Y must be greater than 0.");
                return 1.0;
            }

            double value = Math.Exp(
                Helpers.LogGamma(x)
                + Helpers.LogGamma(y)
                - Helpers.LogGamma(x + y));

            return value;
        }

        public static double r8poly_value_horner(int m, double[] c, double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
        //
        //  Discussion:
        //
        //    The polynomial 
        //
        //      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
        //
        //    is to be evaluated at the value X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the degree of the polynomial.
        //
        //    Input, double C[M+1], the coefficients of the polynomial.
        //    A[0] is the constant term.
        //
        //    Input, double X, the point at which the polynomial is to be evaluated.
        //
        //    Output, double R8POLY_VALUE_HORNER, the value of the polynomial at X.
        //
        {
            double value = c[m];

            for (int i = m - 1; 0 <= i; i--)
            {
                value = value * x + c[i];
            }

            return value;
        }

        public static double r8_factorial(int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FACTORIAL computes the factorial of N.
        //
        //  Discussion:
        //
        //    factorial ( N ) = product ( 1 <= I <= N ) I
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 January 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument of the factorial function.
        //    If N is less than 1, the function value is returned as 1.
        //
        //    Output, double R8_FACTORIAL, the factorial of N.
        //
        {
            double value = 1.0;

            for (int i = 1; i <= n; i++)
            {
                value = value * (double) (i);
            }

            return value;
        }
        
        public static double r8_gamma(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA evaluates Gamma(X) for a real argument.
        //
        //  Discussion:
        //
        //    This routine calculates the gamma function for a real argument X.
        //
        //    Computation is based on an algorithm outlined in reference 1.
        //    The program uses rational functions that approximate the gamma
        //    function to at least 20 significant decimal digits.  Coefficients
        //    for the approximation over the interval (1,2) are unpublished.
        //    Those for the approximation for 12 <= X are from reference 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody, Laura Stoltz.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody,
        //    An Overview of Software Development for Special Functions,
        //    in Numerical Analysis Dundee, 1975,
        //    edited by GA Watson,
        //    Lecture Notes in Mathematics 506,
        //    Springer, 1976.
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
        //    Charles Mesztenyi, John Rice, Henry Thatcher,
        //    Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968,
        //    LC: QA297.C64.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the function.
        //
        //    Output, double R8_GAMMA, the value of the function.
        //
        {
            //
            //  Coefficients for minimax approximation over (12, INF).
            //
            double[] c =  {
                -1.910444077728E-03,
                8.4171387781295E-04,
                -5.952379913043012E-04,
                7.93650793500350248E-04,
                -2.777777777777681622553E-03,
                8.333333333333333331554247E-02,
                5.7083835261E-03
            }
            ;
            double eps = 2.22E-16;
            double half = 0.5;
            int i;
            double one = 1.0;
            double[] p =  {
                -1.71618513886549492533811E+00,
                2.47656508055759199108314E+01,
                -3.79804256470945635097577E+02,
                6.29331155312818442661052E+02,
                8.66966202790413211295064E+02,
                -3.14512729688483675254357E+04,
                -3.61444134186911729807069E+04,
                6.64561438202405440627855E+04
            }
            ;
            const double r8_pi = 3.1415926535897932384626434;
            double[] q =  {
                -3.08402300119738975254353E+01,
                3.15350626979604161529144E+02,
                -1.01515636749021914166146E+03,
                -3.10777167157231109440444E+03,
                2.25381184209801510330112E+04,
                4.75584627752788110767815E+03,
                -1.34659959864969306392456E+05,
                -1.15132259675553483497211E+05
            }
            ;
            double res;
            double sqrtpi = 0.9189385332046727417803297;
            double twelve = 12.0;
            double two = 2.0;
            double value;
            double xbig = 171.624;
            double xinf = 1.79E+308;
            double xminin = 2.23E-308;
            double y1;
            double zero = 0.0;
            ;

            bool parity = false;
            double fact = one;
            int n = 0;
            double y = x;
            //
            //  Argument is negative.
            //
            if (y <= zero)
            {
                y = -x;
                y1 = (double) (int) (y);
                res = y - y1;

                if (res != zero)
                {
                    if (y1 != (double) (int) (y1 * half) * two)
                    {
                        parity = true;
                    }

                    fact = -r8_pi / Math.Sin(r8_pi * res);
                    y = y + one;
                }
                else
                {
                    res = xinf;
                    value = res;
                    return value;
                }
            }

            //
            //  Argument is positive.
            //
            if (y < eps)
            {
                //
                //  Argument < EPS.
                //
                if (xminin <= y)
                {
                    res = one / y;
                }
                else
                {
                    res = xinf;
                    value = res;
                    return value;
                }
            }
            else if (y < twelve)
            {
                y1 = y;
                //
                //  0.0 < argument < 1.0.
                //
                double z;
                if (y < one)
                {
                    z = y;
                    y = y + one;
                }
                //
                //  1.0 < argument < 12.0.
                //  Reduce argument if necessary.
                //
                else
                {
                    n = (int) (y) - 1;
                    y = y - (double) (n);
                    z = y - one;
                }

                //
                //  Evaluate approximation for 1.0 < argument < 2.0.
                //
                double xnum = zero;
                double xden = one;
                for (i = 0; i < 8; i++)
                {
                    xnum = (xnum + p[i]) * z;
                    xden = xden * z + q[i];
                }

                res = xnum / xden + one;
                //
                //  Adjust result for case  0.0 < argument < 1.0.
                //
                if (y1 < y)
                {
                    res = res / y1;
                }
                //
                //  Adjust result for case 2.0 < argument < 12.0.
                //
                else if (y < y1)
                {
                    for (i = 1; i <= n; i++)
                    {
                        res = res * y;
                        y = y + one;
                    }
                }
            }
            else
            {
                //
                //  Evaluate for 12.0 <= argument.
                //
                if (y <= xbig)
                {
                    double ysq = y * y;
                    double sum = c[6];
                    for (i = 0; i < 6; i++)
                    {
                        sum = sum / ysq + c[i];
                    }

                    sum = sum / y - y + sqrtpi;
                    sum = sum + (y - half) * Math.Log(y);
                    res = Math.Exp(sum);
                }
                else
                {
                    res = xinf;
                    value = res;
                    return value;
                }
            }

            //
            //  Final adjustments and return.
            //
            if (parity)
            {
                res = -res;
            }

            if (fact != one)
            {
                res = fact / res;
            }

            value = res;

            return value;
        }

        public static double r8_gamma_inc(double p, double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA_INC computes the incomplete Gamma function.
        //
        //  Discussion:
        //
        //    GAMMA_INC(P,       0) = 0,
        //    GAMMA_INC(P,Infinity) = 1.
        //
        //    GAMMA_INC(P,X) = Integral ( 0 <= T <= X ) T^(P-1) EXP(-T) DT / GAMMA(P).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by B L Shea.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    B L Shea,
        //    Chi-squared and Incomplete Gamma Integral,
        //    Algorithm AS239,
        //    Applied Statistics,
        //    Volume 37, Number 3, 1988, pages 466-473.
        //
        //  Parameters:
        //
        //    Input, double P, the exponent parameter.
        //    0.0 < P.
        //
        //    Input, double X, the integral limit parameter.
        //    If X is less than or equal to 0, the value is returned as 0.
        //
        //    Output, double R8_GAMMA_INC, the value of the function.
        //
        {
            double a;
            double arg;
            double c;
            double exp_arg_min = -88.0;
            double overflow = 1.0E+37;
            double plimit = 1000.0;
            double pn1;
            double tol = 1.0E-07;
            double xbig = 1.0E+08;

            double value = 0.0;

            if (p <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("R8_GAMMA_INC - Fatal error!");
                Console.WriteLine("  Parameter P <= 0.");
                return (1);
            }

            if (x <= 0.0)
            {
                value = 0.0;
                return value;
            }

            //
            //  Use a normal approximation if PLIMIT < P.
            //
            if (plimit < p)
            {
                pn1 = 3.0 * Math.Sqrt(p) * (Math.Pow(x / p, 1.0 / 3.0)
                    + 1.0 / (9.0 * p) - 1.0);
                value = Normal.normal_01_cdf(pn1);
                return value;
            }

            //
            //  Is X extremely large compared to P?
            //
            if (xbig < x)
            {
                value = 1.0;
                return value;
            }

            //
            //  Use Pearson's series expansion.
            //  (P is not large enough to force overflow in the log of Gamma.
            //
            if (x <= 1.0 || x < p)
            {
                arg = p * Math.Log(x) - x - Helpers.LogGamma(p + 1.0);
                c = 1.0;
                value = 1.0;
                a = p;

                for (;;)
                {
                    a = a + 1.0;
                    c = c * x / a;
                    value = value + c;

                    if (c <= tol)
                    {
                        break;
                    }
                }

                arg = arg + Math.Log(value);

                if (exp_arg_min <= arg)
                {
                    value = Math.Exp(arg);
                }
                else
                {
                    value = 0.0;
                }
            }
            //
            //  Use a continued fraction expansion.
            //
            else
            {
                arg = p * Math.Log(x) - x - Helpers.LogGamma(p);
                a = 1.0 - p;
                double b = a + x + 1.0;
                c = 0.0;
                pn1 = 1.0;
                double pn2 = x;
                double pn3 = x + 1.0;
                double pn4 = x * b;
                value = pn3 / pn4;

                for (;;)
                {
                    a = a + 1.0;
                    b = b + 2.0;
                    c = c + 1.0;
                    double pn5 = b * pn3 - a * c * pn1;
                    double pn6 = b * pn4 - a * c * pn2;

                    if (0.0 < Math.Abs(pn6))
                    {
                        double rn = pn5 / pn6;

                        if (Math.Abs(value - rn) <= Math.Min(tol, tol * rn))
                        {
                            arg = arg + Math.Log(value);

                            if (exp_arg_min <= arg)
                            {
                                value = 1.0 - Math.Exp(arg);
                            }
                            else
                            {
                                value = 1.0;
                            }

                            return value;
                        }

                        value = rn;
                    }

                    pn1 = pn3;
                    pn2 = pn4;
                    pn3 = pn5;
                    pn4 = pn6;
                    //
                    //  Rescale terms in continued fraction if terms are large.
                    //
                    if (overflow <= Math.Abs(pn5))
                    {
                        pn1 = pn1 / overflow;
                        pn2 = pn2 / overflow;
                        pn3 = pn3 / overflow;
                        pn4 = pn4 / overflow;
                    }
                }
            }

            return value;
        }

        public static double r8_gamma_log(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA_LOG calculates the logarithm of GAMMA ( X ) for positive X.
        //
        //  Discussion:
        //
        //    The program uses rational functions that theoretically approximate
        //    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
        //    approximation for 12 < X is from Hart et al, while approximations
        //    for X < 12.0 are similar to those in Cody and Hillstrom, but are
        //    unpublished.
        //
        //    The accuracy achieved depends on the arithmetic system, the compiler,
        //    intrinsic functions, and proper selection of the machine-dependent
        //    constants.
        //
        //  Modified:
        //
        //    12 May 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody, Laura Stoltz.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody, Kenneth Hillstrom,
        //    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
        //    Mathematics of Computation,
        //    Volume 21, 1967, pages 198-203.
        //
        //    Kenneth Hillstrom,
        //    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
        //    May 1969.
        //
        //    Hart, Ward Cheney, Charles Lawson, Maehly,
        //    Charles Mesztenyi, John Rice, Thacher, Witzgall,
        //    Computer Approximations,
        //    Wiley and sons, New York, 1968.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the Gamma function.  X must be positive.
        //
        //    Output, double R8_GAMMA_LOG, the logarithm of the Gamma function of X.
        //    If X <= 0.0, or if overflow would occur, the program returns the
        //    value XINF, the largest representable floating point number.
        //
        //  Machine-dependent constants:
        //
        //  BETA   - radix for the floating-point representation.
        //
        //  MAXEXP - the smallest positive power of BETA that overflows.
        //
        //  XBIG   - largest argument for which LN(GAMMA(X)) is representable
        //           in the machine, i.e., the solution to the equation
        //             LN(GAMMA(XBIG)) = BETA^MAXEXP.
        //
        //  FRTBIG - Rough estimate of the fourth root of XBIG
        //
        {
            double[] c =
            {
                -1.910444077728E-03,
                8.4171387781295E-04,
                -5.952379913043012E-04,
                7.93650793500350248E-04,
                -2.777777777777681622553E-03,
                8.333333333333333331554247E-02,
                5.7083835261E-03
            };
            double corr;
            double d1 = -5.772156649015328605195174E-01;
            double d2 = 4.227843350984671393993777E-01;
            double d4 = 1.791759469228055000094023;
            double frtbig = 1.42E+09;
            int i;
            double[] p1 =
            {
                4.945235359296727046734888,
                2.018112620856775083915565E+02,
                2.290838373831346393026739E+03,
                1.131967205903380828685045E+04,
                2.855724635671635335736389E+04,
                3.848496228443793359990269E+04,
                2.637748787624195437963534E+04,
                7.225813979700288197698961E+03
            };
            double[] p2 =
            {
                4.974607845568932035012064,
                5.424138599891070494101986E+02,
                1.550693864978364947665077E+04,
                1.847932904445632425417223E+05,
                1.088204769468828767498470E+06,
                3.338152967987029735917223E+06,
                5.106661678927352456275255E+06,
                3.074109054850539556250927E+06
            };
            double[] p4 =
            {
                1.474502166059939948905062E+04,
                2.426813369486704502836312E+06,
                1.214755574045093227939592E+08,
                2.663432449630976949898078E+09,
                2.940378956634553899906876E+010,
                1.702665737765398868392998E+011,
                4.926125793377430887588120E+011,
                5.606251856223951465078242E+011
            };
            double pnt68 = 0.6796875;
            double[] q1 =
            {
                6.748212550303777196073036E+01,
                1.113332393857199323513008E+03,
                7.738757056935398733233834E+03,
                2.763987074403340708898585E+04,
                5.499310206226157329794414E+04,
                6.161122180066002127833352E+04,
                3.635127591501940507276287E+04,
                8.785536302431013170870835E+03
            };
            double[] q2 =
            {
                1.830328399370592604055942E+02,
                7.765049321445005871323047E+03,
                1.331903827966074194402448E+05,
                1.136705821321969608938755E+06,
                5.267964117437946917577538E+06,
                1.346701454311101692290052E+07,
                1.782736530353274213975932E+07,
                9.533095591844353613395747E+06
            };
            double[] q4 =
            {
                2.690530175870899333379843E+03,
                6.393885654300092398984238E+05,
                4.135599930241388052042842E+07,
                1.120872109616147941376570E+09,
                1.488613728678813811542398E+010,
                1.016803586272438228077304E+011,
                3.417476345507377132798597E+011,
                4.463158187419713286462081E+011
            };
            const double r8_huge = 1.0E+30;
            double res;
            double sqrtpi = 0.9189385332046727417803297;
            double xbig = 4.08E+36;
            double xden;
            double xm2;
            double xnum;
            //
            //  Return immediately if the argument is out of range.
            //
            if (x <= 0.0 || xbig < x)
            {
                return r8_huge;
            }

            if (x <= double.Epsilon)
            {
                res = -Math.Log(x);
            }
            else if (x <= 1.5)
            {
                double xm1;
                if (x < pnt68)
                {
                    corr = -Math.Log(x);
                    xm1 = x;
                }
                else
                {
                    corr = 0.0;
                    xm1 = (x - 0.5) - 0.5;
                }

                if (x <= 0.5 || pnt68 <= x)
                {
                    xden = 1.0;
                    xnum = 0.0;

                    for (i = 0; i < 8; i++)
                    {
                        xnum = xnum * xm1 + p1[i];
                        xden = xden * xm1 + q1[i];
                    }

                    res = corr + (xm1 * (d1 + xm1 * (xnum / xden)));
                }
                else
                {
                    xm2 = (x - 0.5) - 0.5;
                    xden = 1.0;
                    xnum = 0.0;
                    for (i = 0; i < 8; i++)
                    {
                        xnum = xnum * xm2 + p2[i];
                        xden = xden * xm2 + q2[i];
                    }

                    res = corr + xm2 * (d2 + xm2 * (xnum / xden));

                }
            }
            else if (x <= 4.0)
            {
                xm2 = x - 2.0;
                xden = 1.0;
                xnum = 0.0;
                for (i = 0; i < 8; i++)
                {
                    xnum = xnum * xm2 + p2[i];
                    xden = xden * xm2 + q2[i];
                }

                res = xm2 * (d2 + xm2 * (xnum / xden));
            }
            else if (x <= 12.0)
            {
                double xm4 = x - 4.0;
                xden = -1.0;
                xnum = 0.0;
                for (i = 0; i < 8; i++)
                {
                    xnum = xnum * xm4 + p4[i];
                    xden = xden * xm4 + q4[i];
                }

                res = d4 + xm4 * (xnum / xden);
            }
            else
            {
                res = 0.0;

                if (x <= frtbig)
                {

                    res = c[6];
                    double xsq = x * x;

                    for (i = 0; i < 6; i++)
                    {
                        res = res / xsq + c[i];
                    }

                }

                res = res / x;
                corr = Math.Log(x);
                res = res + sqrtpi - 0.5 * corr;
                res = res + x * (corr - 1.0);

            }

            return res;
        }

        public static double r8_gamma_log_int(int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA_LOG_INT computes the logarithm of Gamma of an integer N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument of the logarithm of the Gamma function.
        //    0 < N.
        //
        //    Output, double R8_GAMMA_LOG_INT, the logarithm of
        //    the Gamma function of N.
        //
        {
            if (n <= 0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("R8_GAMMA_LOG_INT - Fatal error!");
                Console.WriteLine("  Illegal input value of N = " + n + "");
                Console.WriteLine("  But N must be strictly positive.");
                return (1);
            }

            double value = Helpers.LogGamma((double) (n));

            return value;
        }
        
        public static double r8vec_dot ( int n, double[] a1, double[] a2 )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DOT computes the dot product of a pair of R8VEC's.
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
        //    Input, double A1[N], A2[N], the two vectors to be considered.
        //
        //    Output, double R8VEC_DOT, the dot product of the vectors.
        //
        {
            double value = 0.0;
            for (int i = 0; i < n; i++ )
            {
                value = value + a1[i] * a2[i];
            }

            return value;
        }
        
        public static double r8vec_length ( int dim_num, double[] x )
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
            for (int i = 0; i < dim_num; i++ )
            {
                value = value + Math.Pow ( x[i], 2 );
            }
            value = Math.Sqrt ( value );

            return value;
        }
        
        public static double r8_csc ( double theta )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CSC returns the cosecant of X.
        //
        //  Discussion:
        //
        //    R8_CSC ( THETA ) = 1.0 / SIN ( THETA )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double THETA, the angle, in radians, whose cosecant is desired.
        //    It must be the case that SIN ( THETA ) is not zero.
        //
        //    Output, double R8_CSC, the cosecant of THETA.
        //
        {
            double value = Math.Sin ( theta );

            if ( value == 0.0 )
            {
                Console.WriteLine(" ");
                Console.WriteLine("R8_CSC - Fatal error!");
                Console.WriteLine("  Cosecant undefined for THETA = " + theta + "");
                return ( 1 );
            }

            value = 1.0 / value;

            return value;
        }
        
        public static double r8_sign ( double x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIGN returns the sign of an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose sign is desired.
        //
        //    Output, double R8_SIGN, the sign of X.
        //
        {
            if ( x < 0.0 )
            {
                return ( -1.0 );
            }
            else
            {
                return ( 1.0 );
            }
        }

    }
}