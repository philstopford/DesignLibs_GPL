﻿using System;

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
        public static r8 s_to_r8 ( string s )
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

            int nchar = s_len_trim ( s );
            int isgn = 1;
            double rtop = 0.0;
            double rbot = 1.0;
            int jsgn = 1;
            int jtop = 0;
            int jbot = 1;
            int ihave = 1;
            int iterm = 0;

            for ( ; ; )
            {
                char c = s[ret.lchar+1];
                ret.lchar = ret.lchar + 1;
                //
                //  Blank or TAB character.
                //
                if ( c == ' ' || c == TAB )
                {
                    if ( ihave == 2 )
                    {
                    }
                    else if ( ihave == 6 || ihave == 7 )
                    {
                        iterm = 1;
                    }
                    else if ( 1 < ihave )
                    {
                        ihave = 11;
                    }
                }
                //
                //  Comma.
                //
                else if ( c == ',' || c == ';' )
                {
                    if ( ihave != 1 )
                    {
                        iterm = 1;
                        ihave = 12;
                        ret.lchar = ret.lchar + 1;
                    }
                }
                //
                //  Minus sign.
                //
                else if ( c == '-' )
                {
                    if ( ihave == 1 )
                    {
                        ihave = 2;
                        isgn = -1;
                    }
                    else if ( ihave == 6 )
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
                else if ( c == '+' )
                {
                    if ( ihave == 1 )
                    {
                        ihave = 2;
                    }
                    else if ( ihave == 6 )
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
                else if ( c == '.' )
                {
                    if ( ihave < 4 )
                    {
                        ihave = 4;
                    }
                    else if ( 6 <= ihave && ihave <= 8 )
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
                else if ( ( Char.ToUpper(c) == 'E' ) || ( Char.ToUpper(c) == 'D' ) )
                {
                    if ( ihave < 6 )
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
                else if ( ihave < 11 && '0' <= c && c <= '9' )
                {
                    if ( ihave <= 2 )
                    {
                        ihave = 3;
                    }
                    else if ( ihave == 4 )
                    {
                        ihave = 5;
                    }
                    else if ( ihave == 6 || ihave == 7 )
                    {
                        ihave = 8;
                    }
                    else if ( ihave == 9 )
                    {
                        ihave = 10;
                    }

                    int ndig = ch_to_digit ( c );

                    if ( ihave == 3 )
                    {
                        rtop = 10.0 * rtop + ( double ) ndig;
                    }
                    else if ( ihave == 5 )
                    {
                        rtop = 10.0 * rtop + ( double ) ndig;
                        rbot = 10.0 * rbot;
                    }
                    else if ( ihave == 8 )
                    {
                        jtop = 10 * jtop + ndig;
                    }
                    else if ( ihave == 10 )
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
                if ( iterm == 1 || nchar <= ret.lchar + 1 )
                {
                    break;
                }

            }
            //
            //  If we haven't seen a terminator, and we have examined the
            //  entire string, then we're done, and LCHAR is equal to NCHAR.
            //
            if ( iterm != 1 && (ret.lchar) + 1 == nchar )
            {
                ret.lchar = nchar;
            }
            //
            //  Number seems to have terminated.  Have we got a legal number?
            //  Not if we terminated in states 1, 2, 6 or 7!
            //
            if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
            {
                ret.error = true;
                return ret;
            }
            //
            //  Number seems OK.  Form it.
            //
            if ( jtop == 0 )
            {
                rexp = 1.0;
            }
            else
            {
                if ( jbot == 1 )
                {
                    rexp = Math.Pow ( 10.0, jsgn * jtop );
                }
                else
                {
                    rexp = jsgn * jtop;
                    rexp = rexp / jbot;
                    rexp = Math.Pow ( 10.0, rexp );
                }
            }

            ret.val = isgn * rexp * rtop / rbot;

            return ret;
        }

        
        public static r8vec s_to_r8vec ( string s, int n )
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

            r8vec ret = new r8vec() {rvec = new double[n]} ;

            for (int i = 0; i < n; i++ )
            {
                r8 res = s_to_r8 ( s.Substring(begin,length));

                ret.rvec[i] = res.val;
                int lchar = res.lchar;
                
                if ( ret.error )
                {
                    return ret;
                }
                begin = begin + lchar;
                length = length - lchar;
            }

            return ret;
        }
        
        public static double r8vec_sum ( int n, double[] a )
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
	        for (int i = 0; i < n; i++ )
	        {
		        value = value + a[i];
	        }
	        return value;
        }
        
        
		public static int[] r8vec_sort_heap_index_a_new ( int n, double[] a )

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

			for (int i = 1; i <= n; i++ )
			{
				indx[i-1] = i;
			}

			int l = n / 2 + 1;
			int ir = n;

			for ( ; ; )
			{
				double aval;
				int indxt;
				if ( 1 < l )
				{
					l = l - 1;
					indxt = indx[l-1];
					aval = a[indxt-1];
				}
				else
				{
					indxt = indx[ir-1];
					aval = a[indxt-1];
					indx[ir-1] = indx[0];
					ir = ir - 1;

					if ( ir == 1 )
					{
						indx[0] = indxt;
						for (int  i = 0; i < n; i++ )
						{
							indx[i] = indx[i] - 1;
						}
						break;
					}
				}

				int i2 = l;
				int j = l + l;

				while ( j <= ir )
				{
					if ( j < ir )
					{
						if ( a[indx[j-1]-1] < a[indx[j]-1] )
						{
							j = j + 1;
						}
					}

					if ( aval < a[indx[j-1]-1] )
					{
						indx[i2-1] = indx[j-1];
						i2 = j;
						j = j + j;
					}
					else
					{
						j = ir + 1;
					}
				}
				indx[i2-1] = indxt;
			}
			return indx;
		}
		
		public static long r8_nint ( double x )
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
			long value =  (long) ( Math.Abs ( x ) + 0.5 );

			if ( x < 0.0 )
			{
				value = -value;
			}

			return value;
		}
    }
}