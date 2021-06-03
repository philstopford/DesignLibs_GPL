using System;
using System.ComponentModel.Design;
using System.Linq;
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

        public static double r8_error_f(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ERROR_F evaluates the error function ERF.
        //
        //  Discussion:
        //
        //    Since some compilers already supply a routine named ERF which evaluates
        //    the error function, this routine has been given a distinct, if
        //    somewhat unnatural, name.
        //
        //    The function is defined by:
        //
        //      ERF(X) = ( 2 / sqrt ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( - T^2 ) dT.
        //
        //    Properties of the function include:
        //
        //      Limit ( X -> -Infinity ) ERF(X) =          -1.0;
        //                               ERF(0) =           0.0;
        //                               ERF(0.476936...) = 0.5;
        //      Limit ( X -> +Infinity ) ERF(X) =          +1.0.
        //
        //      0.5 * ( ERF(X/sqrt(2)) + 1 ) = Normal_01_CDF(X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 November 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 versino by William Cody.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody,
        //    "Rational Chebyshev Approximations for the Error Function",
        //    Mathematics of Computation,
        //    1969, pages 631-638.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the error function.
        //
        //    Output, double R8_ERROR_F, the value of the error function.
        //
        {
            double[] a =  {
                3.16112374387056560,
                1.13864154151050156E+02,
                3.77485237685302021E+02,
                3.20937758913846947E+03,
                1.85777706184603153E-01
            }
            ;
            double[] b =  {
                2.36012909523441209E+01,
                2.44024637934444173E+02,
                1.28261652607737228E+03,
                2.84423683343917062E+03
            }
            ;
            double[] c =  {
                5.64188496988670089E-01,
                8.88314979438837594,
                6.61191906371416295E+01,
                2.98635138197400131E+02,
                8.81952221241769090E+02,
                1.71204761263407058E+03,
                2.05107837782607147E+03,
                1.23033935479799725E+03,
                2.15311535474403846E-08
            }
            ;
            double[] d =  {
                1.57449261107098347E+01,
                1.17693950891312499E+02,
                5.37181101862009858E+02,
                1.62138957456669019E+03,
                3.29079923573345963E+03,
                4.36261909014324716E+03,
                3.43936767414372164E+03,
                1.23033935480374942E+03
            }
            ;
            double del;
            double erfxd;
            int i;
            double[] p =  {
                3.05326634961232344E-01,
                3.60344899949804439E-01,
                1.25781726111229246E-01,
                1.60837851487422766E-02,
                6.58749161529837803E-04,
                1.63153871373020978E-02
            }
            ;
            double[] q =  {
                2.56852019228982242,
                1.87295284992346047,
                5.27905102951428412E-01,
                6.05183413124413191E-02,
                2.33520497626869185E-03
            }
            ;
            double sqrpi = 0.56418958354775628695;
            double thresh = 0.46875;
            double xabs;
            double xbig = 26.543;
            double xden;
            double xnum;
            double xsmall = 1.11E-16;
            double xsq;

            xabs = Math.Abs((x));
            //
            //  Evaluate ERF(X) for |X| <= 0.46875.
            //
            if (xabs <= thresh)
            {
                if (xsmall < xabs)
                {
                    xsq = xabs * xabs;
                }
                else
                {
                    xsq = 0.0;
                }

                xnum = a[4] * xsq;
                xden = xsq;

                for (i = 0; i < 3; i++)
                {
                    xnum = (xnum + a[i]) * xsq;
                    xden = (xden + b[i]) * xsq;
                }

                erfxd = x * (xnum + a[3]) / (xden + b[3]);
            }
            //
            //  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
            //
            else if (xabs <= 4.0)
            {
                xnum = c[8] * xabs;
                xden = xabs;
                for (i = 0; i < 7; i++)
                {
                    xnum = (xnum + c[i]) * xabs;
                    xden = (xden + d[i]) * xabs;
                }

                erfxd = (xnum + c[7]) / (xden + d[7]);
                xsq = ((double) ((int) (xabs * 16.0))) / 16.0;
                del = (xabs - xsq) * (xabs + xsq);
                erfxd = Math.Exp(-xsq * xsq) * Math.Exp(-del) * erfxd;

                erfxd = (0.5 - erfxd) + 0.5;

                if (x < 0.0)
                {
                    erfxd = -erfxd;
                }
            }
            //
            //  Evaluate ERFC(X) for 4.0 < |X|.
            //
            else
            {
                if (xbig <= xabs)
                {
                    if (0.0 < x)
                    {
                        erfxd = 1.0;
                    }
                    else
                    {
                        erfxd = -1.0;
                    }
                }
                else
                {
                    xsq = 1.0 / (xabs * xabs);

                    xnum = p[5] * xsq;
                    xden = xsq;

                    for (i = 0; i < 4; i++)
                    {
                        xnum = (xnum + p[i]) * xsq;
                        xden = (xden + q[i]) * xsq;
                    }

                    erfxd = xsq * (xnum + p[4]) / (xden + q[4]);
                    erfxd = (sqrpi - erfxd) / xabs;
                    xsq = ((double) ((int) (xabs * 16.0))) / 16.0;
                    del = (xabs - xsq) * (xabs + xsq);
                    erfxd = Math.Exp(-xsq * xsq) * Math.Exp(-del) * erfxd;

                    erfxd = (0.5 - erfxd) + 0.5;

                    if (x < 0.0)
                    {
                        erfxd = -erfxd;
                    }
                }
            }

            return erfxd;
        }

        public static double r8_error_f_inverse(double y)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ERROR_F_INVERSE inverts the error function ERF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double Y, the value of the error function.
        //
        //    Output, double R8_ERROR_F_INVERSE, the value X such that
        //    ERF(X) = Y.
        //
        {
            double z = (y + 1.0) / 2.0;

            double x = Normal.normal_01_cdf_inv(z);

            double value = x / Math.Sqrt(2.0);

            return value;
        }
        
        public static double r8_modp ( double x, double y )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MODP returns the nonnegative remainder of R8 division.
        //
        //  Discussion:
        //
        //    If
        //      REM = R8_MODP ( X, Y )
        //      RMULT = ( X - REM ) / Y
        //    then
        //      X = Y * RMULT + REM
        //    where REM is always nonnegative.
        //
        //    The MOD function computes a result with the same sign as the
        //    quantity being divided.  Thus, suppose you had an angle A,
        //    and you wanted to ensure that it was between 0 and 360.
        //    Then mod(A,360.0) would do, if A was positive, but if A
        //    was negative, your result would be between -360 and 0.
        //
        //    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
        //
        //  Example:
        //
        //        I         J     MOD R8_MODP  R8_MODP Factorization
        //
        //      107        50       7       7    107 =  2 *  50 + 7
        //      107       -50       7       7    107 = -2 * -50 + 7
        //     -107        50      -7      43   -107 = -3 *  50 + 43
        //     -107       -50      -7      43   -107 =  3 * -50 + 43
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
        //    Input, double X, the number to be divided.
        //
        //    Input, double Y, the number that divides X.
        //
        //    Output, double R8_MODP, the nonnegative remainder when X is divided by Y.
        //
        {
            if ( y == 0.0 )
            {
                Console.WriteLine("");
                Console.WriteLine("R8_MODP - Fatal error!");
                Console.WriteLine("  R8_MODP ( X, Y ) called with Y = " + y + "");
                return ( 1 );
            }

            double value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

            if ( value < 0.0 )
            {
                value = value + Math.Abs ( y );
            }

            return value;
        }

        public static double[] r8row_max(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_MAX returns the maximums of an R8ROW.
        //
        //  Example:
        //
        //    A =
        //      1  2  3
        //      2  6  7
        //
        //    MAX =
        //      3
        //      7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8ROW_MAX[M], the maximums of the rows.
        //
        {
            int i;
            int j;
            double[] amax = new double[m];

            for (i = 0; i < m; i++)
            {
                amax[i] = a[i + 0 * m];

                for (j = 1; j < n; j++)
                {
                    if (amax[i] < a[i + j * m])
                    {
                        amax[i] = a[i + j * m];
                    }
                }
            }

            return amax;
        }

        public static double[] r8row_mean(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_MEAN returns the means of an R8ROW.
        //
        //  Example:
        //
        //    A =
        //      1  2  3
        //      2  6  7
        //
        //    MEAN =
        //      2
        //      5
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8ROW_MEAN[M], the means, or averages, of the rows.
        //
        {
            int i;
            int j;
            double[] mean = new double[m];

            for (i = 0; i < m; i++)
            {
                mean[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    mean[i] = mean[i] + a[i + j * m];
                }

                mean[i] = mean[i] / (double) (n);
            }

            return mean;
        }

        public static double[] r8row_min(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_MIN returns the minimums of an R8ROW.
        //
        //  Example:
        //
        //    A =
        //      1  2  3
        //      2  6  7
        //
        //    MIN =
        //      1
        //      2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8ROW_MIN[M], the minimums of the rows.
        //
        {
            int i;
            int j;
            double[] amin = new double[m];

            for (i = 0; i < m; i++)
            {
                amin[i] = a[i + 0 * m];
                for (j = 1; j < n; j++)
                {
                    if (a[i + j * m] < amin[i])
                    {
                        amin[i] = a[i + j * m];
                    }
                }
            }

            return amin;
        }

        public static double[] r8row_variance(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ROW_VARIANCE returns the variances of an R8ROW.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 October 2004
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
        //    Output, double R8ROW_VARIANCE[M], the variances of the rows.
        //
        {
            int i;
            int j;
            double mean;
            double[] variance = new double[m];

            for (i = 0; i < m; i++)
            {
                mean = 0.0;
                for (j = 0; j < n; j++)
                {
                    mean = mean + a[i + j * m];
                }

                mean = mean / (double) (n);

                variance[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    variance[i] = variance[i] + Math.Pow((a[i + j * m] - mean), 2);
                }

                if (1 < n)
                {
                    variance[i] = variance[i] / (double) (n - 1);
                }
                else
                {
                    variance[i] = 0.0;
                }

            }

            return variance;
        }

        
        public static double r8_fraction ( int i, int j )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FRACTION uses real arithmetic on an integer ratio.
        //
        //  Discussion:
        //
        //    Given integer variables I and J, both FORTRAN and C will evaluate 
        //    an expression such as "I/J" using what is called "integer division",
        //    with the result being an integer.  It is often convenient to express
        //    the parts of a fraction as integers but expect the result to be computed
        //    using real arithmetic.  This function carries out that operation.
        //
        //  Example:
        //
        //       I     J   I/J  R8_FRACTION
        //
        //       1     2     0  0.5
        //       7     4     1  1.75
        //       8     4     2  2.00
        //       9     4     2  2.25
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the arguments.
        //
        //    Output, double R8_FRACTION, the value of the ratio.
        //
        {
            double value = ( double ) ( i ) / ( double ) ( j );

            return value;
        }

        public static double[] r8mat_zero_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_ZERO_NEW returns a new zeroed R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
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
//    Input, int M, N, the number of rows and columns.
//
//    Output, double R8MAT_ZERO[M*N], the new zeroed matrix.
//
        {
            double[] a = new double[m*n];

            for (int j = 0; j < n; j++ )
            {
                for (int i = 0; i < m; i++ )
                {
                    a[i+j*m] = 0.0;
                }
            }
            return a;
        }
        
        public static double r8_error(double x)
//****************************************************************************80
//
//  Purpose:
//
//    R8_ERROR evaluates the error function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ERROR, the error function of X.
//
        {
            double[] erfcs = {
                -0.49046121234691808039984544033376E-01,
                -0.14226120510371364237824741899631,
                +0.10035582187599795575754676712933E-01,
                -0.57687646997674847650827025509167E-03,
                +0.27419931252196061034422160791471E-04,
                -0.11043175507344507604135381295905E-05,
                +0.38488755420345036949961311498174E-07,
                -0.11808582533875466969631751801581E-08,
                +0.32334215826050909646402930953354E-10,
                -0.79910159470045487581607374708595E-12,
                +0.17990725113961455611967245486634E-13,
                -0.37186354878186926382316828209493E-15,
                +0.71035990037142529711689908394666E-17,
                -0.12612455119155225832495424853333E-18,
                +0.20916406941769294369170500266666E-20,
                -0.32539731029314072982364160000000E-22,
                +0.47668672097976748332373333333333E-24,
                -0.65980120782851343155199999999999E-26,
                +0.86550114699637626197333333333333E-28,
                -0.10788925177498064213333333333333E-29,
                +0.12811883993017002666666666666666E-31
            }
            ;
            int nterf = 0;
            double sqeps = 0.0;
            double sqrtpi = 1.77245385090551602729816748334115;
            double value;
            double xbig = 0.0;
            double y;

            if (nterf == 0)
            {
                nterf = inits(erfcs, 21, 0.1 * r8_mach(3));
                xbig = Math.Sqrt(-Math.Log(sqrtpi * r8_mach(3)));
                sqeps = Math.Sqrt(2.0 * r8_mach(3));
            }

            y = Math.Abs(x);

            if (y <= sqeps)
            {
                value = 2.0 * x / sqrtpi;
            }
            else if (y <= 1.0)
            {
                value = x * (1.0 + csevl(2.0 * x * x - 1.0, erfcs, nterf));
            }
            else if (y <= xbig)
            {
                value = 1.0 - r8_errorc(y);
                if (x < 0.0)
                {
                    value = -value;
                }
            }
            else
            {
                value = 1.0;
                if (x < 0.0)
                {
                    value = -value;
                }
            }

            return value;
        }

        public static double r8_errorc(double x)
//****************************************************************************80
//
//  Purpose:
//
//    R8_ERRORC evaluates the co-error function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ERRORC, the co-error function of X.
//
        {
            double[] erc2cs = {
                -0.6960134660230950112739150826197E-01,
                -0.4110133936262089348982212084666E-01,
                +0.3914495866689626881561143705244E-02,
                -0.4906395650548979161280935450774E-03,
                +0.7157479001377036380760894141825E-04,
                -0.1153071634131232833808232847912E-04,
                +0.1994670590201997635052314867709E-05,
                -0.3642666471599222873936118430711E-06,
                +0.6944372610005012589931277214633E-07,
                -0.1371220902104366019534605141210E-07,
                +0.2788389661007137131963860348087E-08,
                -0.5814164724331161551864791050316E-09,
                +0.1238920491752753181180168817950E-09,
                -0.2690639145306743432390424937889E-10,
                +0.5942614350847910982444709683840E-11,
                -0.1332386735758119579287754420570E-11,
                +0.3028046806177132017173697243304E-12,
                -0.6966648814941032588795867588954E-13,
                +0.1620854541053922969812893227628E-13,
                -0.3809934465250491999876913057729E-14,
                +0.9040487815978831149368971012975E-15,
                -0.2164006195089607347809812047003E-15,
                +0.5222102233995854984607980244172E-16,
                -0.1269729602364555336372415527780E-16,
                +0.3109145504276197583836227412951E-17,
                -0.7663762920320385524009566714811E-18,
                +0.1900819251362745202536929733290E-18,
                -0.4742207279069039545225655999965E-19,
                +0.1189649200076528382880683078451E-19,
                -0.3000035590325780256845271313066E-20,
                +0.7602993453043246173019385277098E-21,
                -0.1935909447606872881569811049130E-21,
                +0.4951399124773337881000042386773E-22,
                -0.1271807481336371879608621989888E-22,
                +0.3280049600469513043315841652053E-23,
                -0.8492320176822896568924792422399E-24,
                +0.2206917892807560223519879987199E-24,
                -0.5755617245696528498312819507199E-25,
                +0.1506191533639234250354144051199E-25,
                -0.3954502959018796953104285695999E-26,
                +0.1041529704151500979984645051733E-26,
                -0.2751487795278765079450178901333E-27,
                +0.7290058205497557408997703680000E-28,
                -0.1936939645915947804077501098666E-28,
                +0.5160357112051487298370054826666E-29,
                -0.1378419322193094099389644800000E-29,
                +0.3691326793107069042251093333333E-30,
                -0.9909389590624365420653226666666E-31,
                +0.2666491705195388413323946666666E-31
            }
            ;
            double[] erfccs = {
                +0.715179310202924774503697709496E-01,
                -0.265324343376067157558893386681E-01,
                +0.171115397792085588332699194606E-02,
                -0.163751663458517884163746404749E-03,
                +0.198712935005520364995974806758E-04,
                -0.284371241276655508750175183152E-05,
                +0.460616130896313036969379968464E-06,
                -0.822775302587920842057766536366E-07,
                +0.159214187277090112989358340826E-07,
                -0.329507136225284321486631665072E-08,
                +0.722343976040055546581261153890E-09,
                -0.166485581339872959344695966886E-09,
                +0.401039258823766482077671768814E-10,
                -0.100481621442573113272170176283E-10,
                +0.260827591330033380859341009439E-11,
                -0.699111056040402486557697812476E-12,
                +0.192949233326170708624205749803E-12,
                -0.547013118875433106490125085271E-13,
                +0.158966330976269744839084032762E-13,
                -0.472689398019755483920369584290E-14,
                +0.143587337678498478672873997840E-14,
                -0.444951056181735839417250062829E-15,
                +0.140481088476823343737305537466E-15,
                -0.451381838776421089625963281623E-16,
                +0.147452154104513307787018713262E-16,
                -0.489262140694577615436841552532E-17,
                +0.164761214141064673895301522827E-17,
                -0.562681717632940809299928521323E-18,
                +0.194744338223207851429197867821E-18,
                -0.682630564294842072956664144723E-19,
                +0.242198888729864924018301125438E-19,
                -0.869341413350307042563800861857E-20,
                +0.315518034622808557122363401262E-20,
                -0.115737232404960874261239486742E-20,
                +0.428894716160565394623737097442E-21,
                -0.160503074205761685005737770964E-21,
                +0.606329875745380264495069923027E-22,
                -0.231140425169795849098840801367E-22,
                +0.888877854066188552554702955697E-23,
                -0.344726057665137652230718495566E-23,
                +0.134786546020696506827582774181E-23,
                -0.531179407112502173645873201807E-24,
                +0.210934105861978316828954734537E-24,
                -0.843836558792378911598133256738E-25,
                +0.339998252494520890627359576337E-25,
                -0.137945238807324209002238377110E-25,
                +0.563449031183325261513392634811E-26,
                -0.231649043447706544823427752700E-26,
                +0.958446284460181015263158381226E-27,
                -0.399072288033010972624224850193E-27,
                +0.167212922594447736017228709669E-27,
                -0.704599152276601385638803782587E-28,
                +0.297976840286420635412357989444E-28,
                -0.126252246646061929722422632994E-28,
                +0.539543870454248793985299653154E-29,
                -0.238099288253145918675346190062E-29,
                +0.109905283010276157359726683750E-29,
                -0.486771374164496572732518677435E-30,
                +0.152587726411035756763200828211E-30
            }
            ;
            double[] erfcs = {
                -0.49046121234691808039984544033376E-01,
                -0.14226120510371364237824741899631,
                +0.10035582187599795575754676712933E-01,
                -0.57687646997674847650827025509167E-03,
                +0.27419931252196061034422160791471E-04,
                -0.11043175507344507604135381295905E-05,
                +0.38488755420345036949961311498174E-07,
                -0.11808582533875466969631751801581E-08,
                +0.32334215826050909646402930953354E-10,
                -0.79910159470045487581607374708595E-12,
                +0.17990725113961455611967245486634E-13,
                -0.37186354878186926382316828209493E-15,
                +0.71035990037142529711689908394666E-17,
                -0.12612455119155225832495424853333E-18,
                +0.20916406941769294369170500266666E-20,
                -0.32539731029314072982364160000000E-22,
                +0.47668672097976748332373333333333E-24,
                -0.65980120782851343155199999999999E-26,
                +0.86550114699637626197333333333333E-28,
                -0.10788925177498064213333333333333E-29,
                +0.12811883993017002666666666666666E-31
            }
            ;
            double eta;
            int nterc2 = 0;
            int nterf = 0;
            int nterfc = 0;
            double sqeps = 0.0;
            double sqrtpi = 1.77245385090551602729816748334115;
            double value;
            double xmax = 0.0;
            double xsml = 0.0;
            double y;

            if (nterf == 0)
            {
                eta = 0.1 * r8_mach(3);
                nterf = inits(erfcs, 21, eta);
                nterfc = inits(erfccs, 59, eta);
                nterc2 = inits(erc2cs, 49, eta);

                xsml = -Math.Sqrt(-Math.Log(sqrtpi * r8_mach(3)));
                xmax = Math.Sqrt(-Math.Log(sqrtpi * r8_mach(1)));
                xmax = xmax - 0.5 * Math.Log(xmax) / xmax - 0.01;
                sqeps = Math.Sqrt(2.0 * r8_mach(3));
            }

            if (x <= xsml)
            {
                value = 2.0;
                return value;
            }

            if (xmax < x)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_ERRORC - Warning!");
                Console.WriteLine("  X so big that ERFC underflows.");
                value = 0.0;
                return value;
            }

            y = Math.Abs(x);

            if (y < sqeps)
            {
                value = 1.0 - 2.0 * x / sqrtpi;
                return value;
            }
            else if (y <= 1.0)
            {
                value = 1.0 - x * (1.0
                                   + csevl(2.0 * x * x - 1.0, erfcs, nterf));
                return value;
            }

            y = y * y;

            if (y <= 4.0)
            {
                value = Math.Exp(-y) / Math.Abs(x) * (0.5
                                               + csevl((8.0 / y - 5.0) / 3.0, erc2cs, nterc2));
            }
            else
            {
                value = Math.Exp(-y) / Math.Abs(x) * (0.5
                                               + csevl(8.0 / y - 1.0, erfccs, nterfc));
            }

            if (x < 0.0)
            {
                value = 2.0 - value;
            }

            return value;
        }

        public static double r8_mach(int i)
//****************************************************************************80
//
//  Purpose:
//
//    R8_MACH returns double precision real machine constants.
//
//  Discussion:
//
//    Assuming that the internal representation of a double precision real
//    number is in base B, with T the number of base-B digits in the mantissa,
//    and EMIN the smallest possible exponent and EMAX the largest possible 
//    exponent, then
//
//      R8_MACH(1) = B^(EMIN-1), the smallest positive magnitude.
//      R8_MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
//      R8_MACH(3) = B^(-T), the smallest relative spacing.
//      R8_MACH(4) = B^(1-T), the largest relative spacing.
//      R8_MACH(5) = log10(B).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Phyllis Fox, Andrew Hall, Norman Schryer,
//    Algorithm 528:
//    Framework for a Portable Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, Number 2, June 1978, page 176-188.
//
//  Parameters:
//
//    Input, int I, chooses the parameter to be returned.
//    1 <= I <= 5.
//
//    Output, double R8_MACH, the value of the chosen parameter.
//
        {
            double value;

            if (i == 1)
            {
                value = 4.450147717014403E-308;
            }
            else if (i == 2)
            {
                value = 8.988465674311579E+307;
            }
            else if (i == 3)
            {
                value = 1.110223024625157E-016;
            }
            else if (i == 4)
            {
                value = 2.220446049250313E-016;
            }
            else if (i == 5)
            {
                value = 0.301029995663981E+000;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("R8_MACH - Fatal error!");
                Console.WriteLine("  The input argument I is out of bounds.");
                Console.WriteLine("  Legal values satisfy 1 <= I <= 5.");
                Console.WriteLine("  I = " + i + "");
                value = 0.0;
            }

            return value;
        }

        public static double r8_mop ( int i )
//****************************************************************************80
//
//  Purpose:
//
//    R8_MOP returns the I-th power of -1 as an R8 value.
//
//  Discussion:
//
//    An R8 is an double value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the power of -1.
//
//    Output, double R8_MOP, the I-th power of -1.
//
        {
            double value;

            if ( ( i % 2 ) == 0 )
            {
                value = 1.0;
            }
            else
            {
                value = -1.0;
            }

            return value;
        }

    }
}