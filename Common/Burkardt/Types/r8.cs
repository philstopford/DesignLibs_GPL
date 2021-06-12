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
        public static bool r8_is_integer ( double r )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_IS_INTEGER determines if an R8 represents an integer value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the number to be checked.
            //
            //    Output, bool R8_IS_INTEGER, is TRUE if R is an integer value.
            //
        {
            const int i4_huge = 2147483647;
            bool value;

            if ( ( double ) ( i4_huge ) < r )
            {
                value = false;
            }
            else if ( r < - ( double ) ( i4_huge ) )
            {
                value = false;
            }
            else if ( r == ( double ) ( ( int ) ( r ) ) )
            {
                value = true;
            }
            else
            {
                value = false;
            }
            return value;
        }
        
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


        public static double r8_big ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_BIG returns a "big" R8.
            //
            //  Discussion:
            //
            //    The value returned by this function is NOT required to be the
            //    maximum representable R8.
            //    We simply want a "very large" but non-infinite number.
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
            //    Output, double R8_BIG, a "big" R8 value.
            //
        {
            double value;

            value = 1.0E+30;

            return value;
        }
        

        public static double r8_csc(double theta)
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
            double value = Math.Sin(theta);

            if (value == 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("R8_CSC - Fatal error!");
                Console.WriteLine("  Cosecant undefined for THETA = " + theta + "");
                return (1);
            }

            value = 1.0 / value;

            return value;
        }

        public static double r8_sign(double x)
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
            if (x < 0.0)
            {
                return (-1.0);
            }
            else
            {
                return (1.0);
            }
        }


        public static double r8_modp(double x, double y)
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
            if (y == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_MODP - Fatal error!");
                Console.WriteLine("  R8_MODP ( X, Y ) called with Y = " + y + "");
                return (1);
            }

            double value = x - ((double) ((int) (x / y))) * y;

            if (value < 0.0)
            {
                value = value + Math.Abs(y);
            }

            return value;
        }



        public static double r8_fraction(int i, int j)
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
            double value = (double) (i) / (double) (j);

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

        public static double r8_choose(int n, int k)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
            //
            //  Discussion:
            //
            //    The value is calculated in such a way as to avoid overflow and
            //    roundoff.  The calculation is done in R8 arithmetic.
            //
            //    The formula used is:
            //
            //      C(N,K) = N! / ( K! * (N-K)! )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    ML Wolfson, HV Wright,
            //    Algorithm 160:
            //    Combinatorial of M Things Taken N at a Time,
            //    Communications of the ACM,
            //    Volume 6, Number 4, April 1963, page 161.
            //
            //  Parameters:
            //
            //    Input, int N, K, the values of N and K.
            //
            //    Output, double R8_CHOOSE, the number of combinations of N
            //    things taken K at a time.
            //
        {
            int i;
            int mn;
            int mx;
            double value;

            if (k < n - k)
            {
                mn = k;
                mx = n - k;
            }
            else
            {
                mn = n - k;
                mx = k;
            }

            if (mn < 0)
            {
                value = 0.0;
            }
            else if (mn == 0)
            {
                value = 1.0;
            }
            else
            {
                value = (double) (mx + 1);

                for (i = 2; i <= mn; i++)
                {
                    value = (value * (double) (mx + i)) / (double) i;
                }
            }

            return value;
        }

        public static double r8_mop(int i)
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

            if ((i % 2) == 0)
            {
                value = 1.0;
            }
            else
            {
                value = -1.0;
            }

            return value;
        }

        public static int r8_to_bin_even(int nbin, double a, double b, double c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_TO_BIN_EVEN determines the appropriate "bin" for C in [A,B].
            //
            //  Discussion:
            //
            //    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
            //    An initial bin takes everything less than A, and a final bin takes
            //    everything greater than B.
            //
            //  Example:
            //
            //    NBIN = 7, A = 5, B = 15
            //
            //    C   BIN
            //
            //    1    1
            //    3    1
            //    4.9  1
            //    5    2
            //    6    2
            //    7    3
            //    8    3
            //    9.5  4
            //   13    6
            //   14    6
            //   15    6
            //   15.1  7
            //   99    7
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 January 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NBIN, the number of bins.  NBIN is normally
            //    at least 3.  If NBIN is 1 or 2, then everything is assigned to bin 1.
            //
            //    Input, double A, B, the lower and upper limits of the bin
            //    interval.  While A is expected to be less than B, the code should
            //    return useful results if A is actually greater than B.
            //
            //    Input, double C, a value to be placed in a bin.
            //
            //    Output, inte R8_TO_BIN_EVEN, the index of the bin to which C is
            //    assigned.
            //
        {
            double a2;
            double b2;
            int bin;
            bool swap;
            //
            //  Take care of special cases.
            //
            if (nbin < 1)
            {
                bin = 0;
                return bin;
            }
            else if (nbin == 1 || nbin == 2)
            {
                bin = 1;
                return bin;
            }

            if (b == a)
            {
                bin = 0;
                return bin;
            }

            //
            //  If the limits are descending, then we switch them now, and
            //  unswitch the results at the end.
            //
            if (a < b)
            {
                swap = false;
                a2 = a;
                b2 = b;
            }
            else
            {
                swap = true;
                a2 = b;
                b2 = a;
            }

            //
            //  Compute the bin.
            //
            if (c < a2)
            {
                bin = 1;
            }
            else if (c == a2)
            {
                bin = 2;
            }
            else if (c == b2)
            {
                bin = nbin - 1;
            }
            else if (b2 < c)
            {
                bin = nbin;
            }
            else
            {
                bin = 2 + (int) ((double) (nbin - 2) * (c - a2) / (b2 - a2));
                bin = Math.Max(bin, 2);
                bin = Math.Min(bin, nbin - 1);
            }

            //
            //  Reverse the switching.
            //
            if (swap)
            {
                bin = nbin + 1 - bin;
            }

            return bin;
        }

        public static void r8slmat_print ( int m, int n, double[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SLMAT_PRINT prints a strict lower triangular R8MAT.
        //
        //  Example:
        //
        //    M = 5, N = 5
        //    A = (/ 21, 31, 41, 51, 32, 42, 52, 43, 53, 54 /)
        //
        //    21
        //    31 32
        //    41 42 43
        //    51 52 53 54
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[*], the M by N matrix.  Only the strict
        //    lower triangular elements are stored, in column major order.
        //
        //    Input, string TITLE, a title.
        //
        {
            int i;
            int indx;
            int j;
            int jhi;
            int jlo;
            int jmax;
            int nn;

            Console.WriteLine("");
            Console.WriteLine(title + "");

            jmax = Math.Min(n, m - 1);

            nn = 5;

            for (jlo = 1; jlo <= jmax; jlo = jlo + nn)
            {
                jhi = Math.Min(jlo + nn - 1, Math.Min(m - 1, jmax));
                Console.WriteLine("");
                string cout = "  Col   ";
                for (j = jlo; j <= jhi; j++)
                {
                    cout += j.ToString().PadLeft(7) + "       ";
                }

                Console.WriteLine(cout);
                Console.WriteLine("  Row");
                for (i = jlo + 1; i <= m; i++)
                {
                    cout = i.ToString().PadLeft(5) + ":";
                    jhi = Math.Min(jlo + nn - 1, Math.Min(i - 1, jmax));
                    for (j = jlo; j <= jhi; j++)
                    {
                        indx = (j - 1) * m + i - (j * (j + 1)) / 2;
                        cout += " " + a[indx - 1].ToString().PadLeft(12);
                    }

                    Console.WriteLine(cout);
                }
            }

            return;
        }

        public static void r8_mant ( double x, ref int s, ref double r, ref int l )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MANT computes the "mantissa" or "fraction part" of an R8.
        //
        //  Discussion:
        //
        //    X = S * R * 2^L
        //
        //    S is +1 or -1,
        //    R is a real between 1.0 and 2.0,
        //    L is an integer.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the real number to be decomposed.
        //
        //    Output, int &S, the "sign" of the number.
        //    S will be -1 if X is less than 0, and +1 if X is greater
        //    than or equal to zero.
        //
        //    Output, double &R, the mantissa of X.  R will be greater
        //    than or equal to 1, and strictly less than 2.  The one
        //    exception occurs if X is zero, in which case R will also
        //    be zero.
        //
        //    Output, int &L, the integer part of the logarithm (base 2) of X.
        //
        {
            //
            //  Determine the sign.
            //
            if ( x < 0.0 )
            {
                s = -1;
            }
            else
            {
                s = 1;
            }
            //
            //  Set R to the absolute value of X, and L to zero.
            //  Then force R to lie between 1 and 2.
            //
            if ( x < 0.0 )
            {
                r = -x;
            }
            else
            {
                r = x;
            }

            l = 0;
            //
            //  Time to bail out if X is zero.
            //
            if ( x == 0.0 )
            {
                return;
            }

            while ( 2.0 <= r )
            {
                r = r / 2.0;
                l = l + 1;
            }

            while ( r < 1.0 )
            {
                r = r * 2.0;
                l = l - 1;
            }

            return;
        }

        public static double r8_haversine(double a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_HAVERSINE computes the haversine of an angle.
            //
            //  Discussion:
            //
            //    haversine(A) = ( 1 - cos ( A ) ) / 2
            //
            //    The haversine is useful in spherical trigonometry.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 November 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, the angle.
            //
            //    Output, double R8_HAVERSINE, the haversine of the angle.
            //
        {
            double value;

            value = (1.0 - Math.Cos(a)) / 2.0;

            return value;
        }

        public static double r8_heaviside(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_HEAVISIDE evaluates the Heaviside function.
            //
            //  Discussion:
            //
            //    The Heaviside function is 0 for x < 0, 1 for x > 0, and 1/2 for x = 0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 November 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double R8_HEAVISIDE, the value.
            //
        {
            double value;

            if (x < 0.0)
            {
                value = 0.0;
            }
            else if (x == 0.0)
            {
                value = 0.5;
            }
            else
            {
                value = 1.0;
            }

            return value;
        }

        public static double r8_hypot(double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_HYPOT returns the value of sqrt ( X^2 + Y^2 ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, Y, the arguments.
            //
            //    Output, double R8_HYPOT, the value of sqrt ( X^2 + Y^2 ).
            //
        {
            double a;
            double b;
            double value;

            if (Math.Abs(x) < Math.Abs(y))
            {
                a = Math.Abs(y);
                b = Math.Abs(x);
            }
            else
            {
                a = Math.Abs(x);
                b = Math.Abs(y);
            }

            //
            //  A contains the larger value.
            //
            if (a == 0.0)
            {
                value = 0.0;
            }
            else
            {
                value = a * Math.Sqrt(1.0 + (b / a) * (b / a));
            }

            return value;
        }

        public static bool r8_is_in_01(double a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_IN_01 is TRUE if an R8 is in the range [0,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, the value.
            //
            //    Output, bool R8_IS_IN_01, is TRUE if A is between 0 and 1.
            //
        {
            bool value;

            value = (0.0 <= a && a <= 1.0);

            return value;
        }

        public static bool r8_is_inf(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_IS_INF determines if an R8 represents an infinite value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the number to be checked.
            //
            //    Output, bool R8_IS_INF, is TRUE if R is an infinite value.
            //
        {
            const double r8_huge = 1.79769313486231571E+308;
            bool value;

            if (r < 0.0)
            {
                value = (r < -r8_huge);
            }
            else
            {
                value = (r8_huge < r);
            }

            return value;
        }

        public static bool r8_is_insignificant(double r, double s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_IS_INSIGNIFICANT determines if an R8 is insignificant.
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
            //    Input, double R, the number to be compared against.
            //
            //    Input, double S, the number to be compared.
            //
            //    Output, bool R8_IS_INSIGNIFICANT, is TRUE if S is insignificant
            //    compared to R.
            //
        {
            double t;
            double tol;
            bool value;

            value = true;

            t = r + s;
            tol = double.Epsilon * Math.Abs(r);

            if (tol < Math.Abs(r - t))
            {
                value = false;
            }

            return value;
        }
        
        public static bool r8_is_nan(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_IS_NAN determines if an R8 represents a NaN value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the number to be checked.
            //
            //    Output, bool R8_IS_NAN, is TRUE if R is a NaN
            //
        {
            bool value;

            value = (r != r);

            return value;
        }
        
        public static double r8_wrap ( double r, double rlo, double rhi )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_WRAP forces an R8 to lie between given limits by wrapping.
            //
            //  Discussion:
            //
            //    An R8 is a double value.
            //
            //  Example:
            //
            //    RLO = 4.0, RHI = 8.0
            //
            //     R  Value
            //
            //    -2     8
            //    -1     4
            //     0     5
            //     1     6
            //     2     7
            //     3     8
            //     4     4
            //     5     5
            //     6     6
            //     7     7
            //     8     8
            //     9     4
            //    10     5
            //    11     6
            //    12     7
            //    13     8
            //    14     4
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 December 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, a value.
            //
            //    Input, double RLO, RHI, the desired bounds.
            //
            //    Output, double R8_WRAP, a "wrapped" version of the value.
            //
        {
            int n;
            double rhi2;
            double rlo2;
            double rwide;
            double value;
            //
            //  Guarantee RLO2 < RHI2.
            //
            if ( rlo <= rhi )
            {
                rlo2 = rlo;
                rhi2 = rhi;
            }
            else
            {
                rlo2 = rhi;
                rhi2 = rlo;
            }
            //
            //  Find the width.
            //
            rwide = rhi2 - rlo2;
            //
            //  Add enough copies of (RHI2-RLO2) to R so that the
            //  result ends up in the interval RLO2 - RHI2.
            //
            if ( rwide == 0.0 )
            {
                value = rlo;
            }
            else if ( r < rlo2 )
            {
                n = ( int ) ( ( rlo2 - r ) / rwide ) + 1;
                value = r + n * rwide;
                if ( value == rhi )
                {
                    value = rlo;
                }
            }
            else
            {
                n = ( int ) ( ( r - rlo2 ) / rwide );
                value = r - n * rwide;
                if ( value == rlo )
                {
                    value = rhi;
                }
            }
            return value;
        }

    }
}