using System;

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
			int i;
			double value;

			value = 0.0;

			for (i = n - 1; 0 <= i; i--)
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
			int i;
			double value;

			value = 0.0;
			for (i = 0; i < n; i++)
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
			int i;
			int j;
			double[] y = new double[m];

			for (i = 0; i < m; i++)
			{
				y[i] = 0.0;
				for (j = 0; j < n; j++)
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
			int i;
			int j;
			double mean;
			double[] variance = new double[n];

			for (j = 0; j < n; j++)
			{
				mean = 0.0;
				for (i = 0; i < m; i++)
				{
					mean = mean + a[i + j * m];
				}

				mean = mean / (double) (m);

				variance[j] = 0.0;
				for (i = 0; i < m; i++)
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
			int i;
			int j;
			double[] mean = new double[n];

			for (j = 0; j < n; j++)
			{
				mean[j] = 0.0;
				for (i = 0; i < m; i++)
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
			double zsum;
			double zsum_old;

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
				zsum = 0.0;
				n = 0;

				for (;;)
				{
					n = n + 1;
					zsum_old = zsum;
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
	}
}