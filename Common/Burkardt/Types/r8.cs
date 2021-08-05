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
        public static double r8_abs(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_ABS returns the absolute value of an R8.
            //
            //  Discussion:
            //
            //    The C++ math library provides the function fabs() which is preferred.
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
            //    Input, double X, the quantity whose absolute value is desired.
            //
            //    Output, double R8_ABS, the absolute value of X.
            //
        {
            double value;

            if (0.0 <= x)
            {
                value = +x;
            }
            else
            {
                value = -x;
            }

            return value;
        }

        public static double r8_acos(double c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_ACOS computes the arc cosine function, with argument truncation.
            //
            //  Discussion:
            //
            //    If you call your system ACOS routine with an input argument that is
            //    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
            //    This routine truncates arguments outside the range.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2002
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double C, the argument, the cosine of an angle.
            //
            //    Output, double R8_ACOS, an angle whose cosine is C.
            //
        {
            const double r8_pi = 3.141592653589793;
            double value;

            if (c <= -1.0)
            {
                value = r8_pi;
            }
            else if (1.0 <= c)
            {
                value = 0.0;
            }
            else
            {
                value = Math.Acos(c);
            }

            return value;
        }

        public static double r8_acosh(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_ACOSH returns the inverse hyperbolic cosine of a number.
            //
            //  Discussion:
            //
            //    Applying the inverse function
            //
            //      Y = R8_ACOSH(X)
            //
            //    implies that
            //
            //      X = COSH(Y) = 0.5 * ( EXP(Y) + EXP(-Y) ).
            //
            //    For every X greater than or equal to 1, there are two possible
            //    choices Y such that X = COSH(Y), differing only in sign.  It
            //    is usual to resolve this choice by taking the value of ACOSH(X)
            //    to be nonnegative.
            //
            //  Method:
            //
            //    One formula is:
            //
            //      R8_ACOSH = LOG ( X + SQRT ( X^2 - 1.0 ) )
            //
            //    but this formula suffers from roundoff and overflow problems.
            //    The formula used here was recommended by W Kahan, as discussed
            //    by Moler.
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
            //  Reference:
            //
            //    Cleve Moler,
            //    Trigonometry is a Complex Subject,
            //    MATLAB News and Notes,
            //    Summer 1998.
            //
            //  Parameters:
            //
            //    Input, double X, the number whose inverse hyperbolic cosine is desired.
            //    X should be greater than or equal to 1.
            //
            //    Output, double R8_ACOSH, the inverse hyperbolic cosine of X.  The
            //    principal value (that is, the positive value of the two ) is returned.
            //
        {
            double value;

            if (x < 1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_ACOSH - Fatal error!");
                Console.WriteLine("  Argument X must satisfy 1 <= X.");
                Console.WriteLine("  The input X = " + x + "");
                return (1);
            }

            value = 2.0 * Math.Log(
                Math.Sqrt(0.5 * (x + 1.0)) + Math.Sqrt(0.5 * (x - 1.0)));

            return value;
        }

        public static double r8_add(double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_ADD adds two R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 August 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, Y, the numbers to be added.
            //
            //    Output, double R8_ADD, the sum of X and Y.
            //
        {
            double value;

            value = x + y;

            return value;
        }

        public static double r8_atanh(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_ATANH returns the inverse hyperbolic tangent of a number.
            //
            //  Discussion:
            //
            //    Y = R8_ATANH ( X )
            //
            //    implies that
            //
            //    X = TANH(Y) = ( EXP(Y) - EXP(-Y) ) / ( EXP(Y) + EXP(-Y) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 November 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the number whose inverse hyperbolic 
            //    tangent is desired.  The absolute value of X should be less than 
            //    or equal to 1.
            //
            //    Output, double R8_ATANH, the inverse hyperbolic tangent of X.
            //
        {
            const double r8_huge = 1.79769313486231571E+308;
            double value;

            if (x <= -1.0)
            {
                value = -r8_huge;
            }
            else if (1.0 <= x)
            {
                value = +r8_huge;
            }
            else
            {
                value = 0.5 * Math.Log((1.0 + x) / (1.0 - x));
            }

            return value;
        }


        public static double r8_agm(double a, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_AGM computes the arithmetic-geometric mean of A and B.
            //
            //  Discussion:
            //
            //    The AGM is defined for nonnegative A and B.
            //
            //    The AGM of numbers A and B is defined by setting
            //
            //      A(0) = A,
            //      B(0) = B
            //
            //      A(N+1) = ( A(N) + B(N) ) / 2
            //      B(N+1) = sqrt ( A(N) * B(N) )
            //
            //    The two sequences both converge to AGM(A,B).
            //
            //    In Mathematica, the AGM can be evaluated by
            //
            //      ArithmeticGeometricMean [ a, b ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //  Parameters:
            //
            //    Input, double A, B, the arguments whose AGM is to be computed.
            //    0 <= A, 0 <= B.
            //
            //    Output, double R8_AGM, the arithmetic-geometric mean of A and B.
            //
        {
            double a1;
            double a2;
            double b1;
            double b2;
            int it;
            int it_max = 1000;
            double tol;
            double value;

            if (a < 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_AGM - Fatal error!");
                Console.WriteLine("  A < 0.");
                return 1;
            }

            if (b < 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_AGM - Fatal error!");
                Console.WriteLine("  B < 0.");
                return 1;
            }

            if (a == 0.0 || b == 0.0)
            {
                value = 0.0;
                return value;
            }

            if (a == b)
            {
                value = a;
                return value;
            }

            a1 = a;
            b1 = b;

            it = 0;
            tol = 100.0 * double.Epsilon;

            for (;;)
            {
                it = it + 1;

                a2 = (a1 + b1) / 2.0;
                b2 = Math.Sqrt(a1 * b1);

                if (Math.Abs(a2 - b2) <= tol * (a2 + b2))
                {
                    break;
                }

                if (it_max < it)
                {
                    break;
                }

                a1 = a2;
                b1 = b2;
            }

            value = a2;

            return value;
        }

        public static double r8_aint(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_AINT truncates an R8 argument to an integer.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    1 September 2011
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double R8_AINT, the truncated version of X.
            //
        {
            double value;

            if (x < 0.0)
            {
                value = -(double) ((int) (Math.Abs(x)));
            }
            else
            {
                value = (double) ((int) (Math.Abs(x)));
            }

            return value;
        }

        public static double r8_asin(double s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_ASIN computes the arc sine function, with argument truncation.
            //
            //  Discussion:
            //
            //    If you call your system ASIN routine with an input argument that is
            //    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
            //    This routine truncates arguments outside the range.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2002
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double S, the argument, the sine of an angle.
            //
            //    Output, double R8_ASIN, an angle whose sine is S.
            //
        {
            double angle;
            const double r8_pi = 3.141592653589793;

            if (s <= -1.0)
            {
                angle = -r8_pi / 2.0;
            }
            else if (1.0 <= s)
            {
                angle = r8_pi / 2.0;
            }
            else
            {
                angle = Math.Asin(s);
            }

            return angle;
        }

        public static double r8_asinh(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_ASINH returns the inverse hyperbolic sine of a number.
            //
            //  Discussion:
            //
            //    The assertion that:
            //
            //      Y = R8_ASINH ( X )
            //
            //    implies that
            //
            //      X = SINH(Y) = 0.5 * ( EXP(Y) - EXP(-Y) ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 November 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the number whose inverse hyperbolic 
            //    sine is desired.
            //
            //    Output, double R8_ASINH, the inverse hyperbolic sine of X.
            //
        {
            double value;

            value = Math.Log(x + Math.Sqrt(x * x + 1.0));

            return value;
        }

        public static double r8_radians(double degrees)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_RADIANS converts an angle from degree to radian measure.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 May 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double DEGREES, the angle measurement in degrees.
            //
            //    Output, double R8_RADIANS, the angle measurement in radians.
            //
        {
            const double r8_pi = 3.1415926535897932384626434;
            double value;

            value = degrees * r8_pi / 180.0;

            return value;
        }

        public static double r8_relu(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_RELU evaluates the ReLU function of an R8.
            //
            //  Discussion:
            //
            //    An R8 is a double precision real value.
            //
            //    The ReLU function is max(x,0.0).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 January 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double VALUE, the function value.
            //
        {
            double value;

            if (x <= 0.0)
            {
                value = 0.0;
            }
            else
            {
                value = x;
            }

            return value;
        }

        public static double r8_reverse_bytes ( double x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_REVERSE_BYTES reverses the bytes in an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, a value whose bytes are to be reversed.
        //
        //    Output, R8_REVERSE_BYTES, a value with bytes in reverse order;
        //
        {
            byte[] t = BitConverter.GetBytes(x);
            Array.Reverse(t);
            return BitConverter.ToDouble(t);
        }

        public static double r8_secd(double degrees)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_SECD returns the secant of an angle given in degrees.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double DEGREES, the angle in degrees.
            //
            //    Output, double R8_SECD, the secant of the angle.
            //
        {
            const double r8_pi = 3.141592653589793;
            double radians;
            double value;

            radians = r8_pi * (degrees / 180.0);

            value = 1.0 / Math.Cos(radians);

            return value;
        }

        public static double r8_sech(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_SECH evaluates the hyperbolic secant, while avoiding COSH overflow.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument of the function.
            //
            //    Output, double R8_SECH, the value of the function.
            //
        {
            const double log_huge = 80.0;
            double value;

            if (log_huge < Math.Abs(x))
            {
                value = 0.0;
            }
            else
            {
                value = 1.0 / Math.Cosh(x);
            }

            return value;
        }

        public static double r8_sigmoid(double l, double b, double m, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8_sigmoid evaluates the sigmoid or logistic function.
            //
            //  Discussion:
            //
            //    An R8 is a double value.
            //
            //    The sigmoid function is useful for classification problems in
            //    machine learning.  Its value is always between 0 and 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 October 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double l, the maximum value of the function.  This is often 1.
            //
            //    double b, the cutoff value, where the function equals l/2.
            //    This is often 0.
            //
            //    double m, the slope, which determines the steepness of the curve
            //    and the width of the uncertainty interval.  This is often 1.
            //
            //    double x, the argument.
            //
            //  Output:
            //
            //    double r8_sigmoid, the value.
            //
        {
            double value;

            value = l / (1.0 + Math.Exp(-m * (x - b)));

            return value;
        }


        public static void r8_sincos_sum(double a, double b, ref double d, ref double e, ref double f)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_SINCOS_SUM simplifies a*sin(cx)+b*cos(cx).
            //
            //  Discussion:
            //
            //    The expression
            //      a * sin ( c * x ) + b * cos ( c * x )
            //    can be rewritten as
            //      d * sin ( c * x + e )
            //    or
            //      d * cos ( c * x + f ) 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, B, the coefficients in the linear combination.
            //
            //    Output, double &D, &E, &F, the new coefficient, and the shift for
            //    sine or for cosine.
            //
        {
            const double r8_pi = 3.141592653589793E+00;

            d = Math.Sqrt(a * a + b * b);
            e = Math.Atan2(b, a);
            f = Math.Atan2(b, a) - r8_pi / 2.0E+00;
            if (f < -r8_pi)
            {
                f = f + 2.0E+00 * r8_pi;
            }
        }

        public static double r8_sind(double degrees)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_SIND returns the sine of an angle given in degrees.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double DEGREES, the angle in degrees.
            //
            //    Output, double R8_SIND, the sine of the angle.
            //
        {
            const double r8_pi = 3.141592653589793;
            double radians;
            double value;

            radians = r8_pi * (degrees / 180.0);

            value = Math.Sin(radians);

            return value;
        }

        public static double r8_softplus(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_SOFTPLUS evaluates the softplus function of an R8.
            //
            //  Discussion:
            //
            //    An R8 is a double precision real value.
            //
            //    The softplus function is a smoothed (differentiable) version of max(x,0.0).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 September 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double VALUE, the function value.
            //
        {
            double value;

            if (x <= -36.841)
            {
                value = 0.0;
            }
            else if (+36.841 <= x)
            {
                value = x;
            }
            else
            {
                value = Math.Log(1.0 + Math.Exp(x));
            }

            return value;
        }

        public static double r8_sqrt_i4(int i)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_SQRT_I4 returns the square root of an I4 as an R8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int I, the number whose square root is desired.
            //
            //    Output, double R8_SQRT_I4, the value of sqrt(I).
            //
        {
            double value;

            value = Math.Sqrt((double) (i));

            return value;
        }


        public static double r8_square(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8_square returns the square of an R8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double x: the argument.
            //
            //  Output:
            //
            //    double r8_square: the square of x.
            //
        {
            double value = x * x;

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
            double value;

            if (p <= 1.0)
            {
                value = r8_huge;
            }
            else if (p == 2.0)
            {
                value = Math.Pow(Math.PI, 2) / 6.0;
            }
            else if (p == 3.0)
            {
                value = 1.2020569032;
            }
            else if (p == 4.0)
            {
                value = Math.Pow(Math.PI, 4) / 90.0;
            }
            else if (p == 5.0)
            {
                value = 1.0369277551;
            }
            else if (p == 6.0)
            {
                value = Math.Pow(Math.PI, 6) / 945.0;
            }
            else if (p == 7.0)
            {
                value = 1.0083492774;
            }
            else if (p == 8.0)
            {
                value = Math.Pow(Math.PI, 8) / 9450.0;
            }
            else if (p == 9.0)
            {
                value = 1.0020083928;
            }
            else if (p == 10.0)
            {
                value = Math.Pow(Math.PI, 10) / 93555.0;
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



        public static double r8_big()

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

        public static double r8_chop(int place, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CHOP chops an R8 to a given number of binary places.
            //
            //  Example:
            //
            //    3.875 = 2 + 1 + 1/2 + 1/4 + 1/8.
            //
            //    The following values would be returned for the 'chopped' value of
            //    3.875:
            //
            //    PLACE  Value
            //
            //       1      2
            //       2      3     = 2 + 1
            //       3      3.5   = 2 + 1 + 1/2
            //       4      3.75  = 2 + 1 + 1/2 + 1/4
            //       5+     3.875 = 2 + 1 + 1/2 + 1/4 + 1/8
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int PLACE, the number of binary places to preserve.
            //    PLACE = 0 means return the integer part of X.
            //    PLACE = 1 means return the value of X, correct to 1/2.
            //    PLACE = 2 means return the value of X, correct to 1/4.
            //    PLACE = -1 means return the value of X, correct to 2.
            //
            //    Input, double X, the number to be chopped.
            //
            //    Output, double R8_CHOP, the chopped number.
            //
        {
            double fac;
            int temp;
            double value;

            temp = (int) (r8_log_2(x));
            fac = Math.Pow(2.0, (temp - place + 1));
            value = (double) ((int) (x / fac)) * fac;

            return value;
        }

        public static double r8_cosd(double degrees)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_COSD returns the cosine of an angle given in degrees.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double DEGREES, the angle in degrees.
            //
            //    Output, double R8_COSD, the cosine of the angle.
            //
        {
            const double r8_pi = 3.141592653589793;
            double radians;
            double value;

            radians = r8_pi * (degrees / 180.0);

            value = Math.Cos(radians);

            return value;
        }

        public static double r8_cot(double angle)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_COT returns the cotangent of an angle.
            //
            //  Discussion:
            //
            //    R8_COT ( THETA ) = COS ( THETA ) / SIN ( THETA )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double ANGLE, the angle, in radians.
            //
            //    Output, double R8_COT, the cotangent of the angle.
            //
        {
            double value;

            value = Math.Cos(angle) / Math.Sin(angle);

            return value;
        }

        public static double r8_cotd(double degrees)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_COTD returns the cotangent of an angle given in degrees.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double DEGREES, the angle in degrees.
            //
            //    Output, double R8_COTD, the cotangent of the angle.
            //
        {
            const double r8_pi = 3.141592653589793;
            double radians;
            double value;

            radians = r8_pi * (degrees / 180.0);

            value = Math.Cos(radians) / Math.Sin(radians);

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

        public static double r8_cscd(double degrees)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CSCD returns the cosecant of an angle given in degrees.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double DEGREES, the angle in degrees.
            //
            //    Output, double R8_CSCD, the cosecant of the angle.
            //
        {
            const double r8_pi = 3.141592653589793;
            double radians;
            double value;

            radians = r8_pi * (degrees / 180.0);

            value = 1.0 / Math.Sin(radians);

            return value;
        }

        public static double r8_cube_root(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CUBE_ROOT returns the cube root of an R8.
            //
            //  Discussion:
            //
            //    This routine is designed to avoid the possible problems that can occur
            //    when formulas like 0.0^(1/3) or (-1.0)^(1/3) are to be evaluated.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input double X, the number whose cube root is desired.
            //
            //    Output, double R8_CUBE_ROOT, the cube root of X.
            //
        {
            double value;

            if (0.0 < x)
            {
                value = Math.Pow((double) x, (1.0 / 3.0));
            }
            else if (x == 0.0)
            {
                value = 0.0;
            }
            else
            {
                value = -Math.Pow((double) (Math.Abs(x)), (1.0 / 3.0));
            }

            return value;
        }

        public static double r8_degrees(double radians)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_DEGREES converts an angle from radian to degree measure.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 May 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double RADIANS, the angle measurement in radians.
            //
            //    Output, double R8_DEGREES, the angle measurement in degrees.
            //
        {
            const double r8_pi = 3.1415926535897932384626434;
            double value;

            value = radians * 180.0 / r8_pi;

            return value;
        }

        public static double r8_diff(double x, double y, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_DIFF computes (X-Y) to a specified accuracy.
            //
            //  Discussion:
            //
            //    The user controls how many binary digits of accuracy
            //    are to be used.
            //
            //    N determines the accuracy of the value.  If N = 10,
            //    for example, only 11 binary places will be used in the arithmetic.
            //    In general, only N+1 binary places will be used.
            //
            //    N may be zero.  However, a negative value of N should
            //    not be used, since this will cause both X and Y to look like 0.
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
            //    Input, double X, Y, the two values whose difference is desired.
            //
            //    Input, int N, the number of binary digits to use.
            //
            //    Output, double R8_DIFF, the value of X-Y.
            //
        {
            double cx;
            double cy;
            double pow2;
            double size;
            double value;

            if (x == y)
            {
                value = 0.0;
                return value;
            }

            pow2 = Math.Pow(2.0, n);
            //
            //  Compute the magnitude of X and Y, and take the larger of the
            //  two.  At least one of the two values is not zero//
            //
            size = Math.Max(Math.Abs(x), Math.Abs(y));
            //
            //  Make normalized copies of X and Y.  One of the two values will
            //  actually be equal to 1.
            //
            cx = x / size;
            cy = y / size;
            //
            //  Here's where rounding comes in.  We know that the larger of the
            //  the two values equals 1.  We multiply both values by 2^N,
            //  where N+1 is the number of binary digits of accuracy we want
            //  to use, truncate the values, and divide back by 2^N.
            //
            cx = (double) ((int) (cx * pow2 + 0.5 * r8_sign(cx))) / pow2;
            cy = (double) ((int) (cy * pow2 + 0.5 * r8_sign(cy))) / pow2;
            //
            //  Take the difference now.
            //
            value = cx - cy;
            //
            //  Undo the scaling.
            //
            value = value * size;

            return value;
        }

        public static int r8_digit(double x, int idigit)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_DIGIT returns a particular decimal digit of an R8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the number whose IDIGIT-th decimal digit is desired.
            //    Note that if X is zero, all digits will be returned as 0.
            //
            //    Input, int IDIGIT, the position of the desired decimal digit.
            //    A value of 1 means the leading digit, a value of 2 the second digit
            //    and so on.
            //
            //    Output, int R8_DIGIT, the value of the IDIGIT-th decimal digit of X.
            //
        {
            int digit;
            int i;
            int ival = 0;

            if (x == 0.0)
            {
                digit = 0;
                return digit;
            }

            if (idigit <= 0)
            {
                digit = 0;
                return digit;
            }

            //
            //  Force X to lie between 1 and 10.
            //
            x = Math.Abs(x);

            while (x < 1.0)
            {
                x = x * 10.0;
            }

            while (10.0 <= x)
            {
                x = x / 10.0;
            }

            for (i = 1; i <= idigit; i++)
            {
                ival = (int) (x);
                x = (x - (double) ival) * 10.0;
            }

            digit = ival;

            return digit;
        }

        public static double r8_divide_i4(int i, int j)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_DIVIDE_I4 returns an I4 fraction as an R8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int I, J, the numerator and denominator.
            //
            //    Output, double R8_DIVIDE_I4, the value of (I/J).
            //
        {
            double value;

            value = (double) (i) / (double) (j);

            return value;
        }

        public static double r8_e()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_E returns the value of the base of the natural logarithm system.
            //
            //  Definition:
            //
            //    E = Limit ( N -> +oo ) ( 1 + 1 / N )^N
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double R8_E, the base of the natural logarithm system.
            //
        {
            const double r8_e_save = 2.718281828459045235360287;
            double value;

            value = r8_e_save;

            return value;
        }

        public static double r8_epsilon()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_EPSILON returns the R8 roundoff unit.
            //
            //  Discussion:
            //
            //    The roundoff unit is a number R which is a power of 2 with the
            //    property that, to the precision of the computer's arithmetic,
            //      1 < 1 + R
            //    but
            //      1 = ( 1 + R / 2 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double R8_EPSILON, the R8 round-off unit.
            //
        {
            return double.Epsilon;
        }

        public class r8EpsilonData
        {
            public double value = 0.0;
            
        }
        
        public static double r8_epsilon_compute(ref r8EpsilonData data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_EPSILON_COMPUTE computes the R8 roundoff unit.
            //
            //  Discussion:
            //
            //    The roundoff unit is a number R which is a power of 2 with the
            //    property that, to the precision of the computer's arithmetic,
            //      1 < 1 + R
            //    but
            //      1 = ( 1 + R / 2 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double R8_EPSILON_COMPUTE, the R8 round-off unit.
            //
        {
            double one;
            double temp;
            double test;

            if (data.value == 0.0)
            {
                one = (double) (1);

                data.value = one;
                temp = data.value / 2.0;
                test = r8_add(one, temp);

                while (one < test)
                {
                    data.value = temp;
                    temp = data.value / 2.0;
                    test = r8_add(one, temp);
                }
            }

            return data.value;
        }

        public static double r8_exp(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_EXP computes the exponential function, avoiding overflow and underflow.
            //
            //  Discussion:
            //
            //    For arguments of very large magnitude, the evaluation of the
            //    exponential function can cause computational problems.  Some languages
            //    and compilers may return an infinite value or a "Not-a-Number".  
            //    An alternative, when dealing with a wide range of inputs, is simply
            //    to truncate the calculation for arguments whose magnitude is too large.
            //    Whether this is the right or convenient approach depends on the problem
            //    you are dealing with, and whether or not you really need accurate
            //    results for large magnitude inputs, or you just want your code to
            //    stop crashing.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument of the exponential function.
            //
            //    Output, double R8_EXP, the value of exp ( X ).
            //
        {
            const double r8_big = 1.0E+30;
            const double r8_log_max = +69.0776;
            const double r8_log_min = -69.0776;
            double value;

            if (x <= r8_log_min)
            {
                value = 0.0;
            }
            else if (x < r8_log_max)
            {
                value = Math.Exp(x);
            }
            else
            {
                value = r8_big;
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


        public static double r8_cas(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CAS returns the "casine" of an R8.
            //
            //  Discussion:
            //
            //    The "casine", used in the discrete Hartley transform, is abbreviated
            //    CAS(X), and defined by:
            //
            //      CAS(X) = cos ( X ) + sin( X )
            //             = sqrt ( 2 ) * sin ( X + pi/4 )
            //             = sqrt ( 2 ) * cos ( X - pi/4 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the number whose casine is desired.
            //
            //    Output, double R8_CAS, the casine of X, which will be between
            //    plus or minus the square root of 2.
            //
        {
            double value;

            value = Math.Cos(x) + Math.Sin(x);

            return value;
        }

        public static double r8_ceiling(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CEILING rounds an R8 up to the nearest integral R8.
            //
            //  Example:
            //
            //    X        R8_CEILING(X)
            //
            //   -1.1      -1.0
            //   -1.0      -1.0
            //   -0.9       0.0
            //   -0.1       0.0
            //    0.0       0.0
            //    0.1       1.0
            //    0.9       1.0
            //    1.0       1.0
            //    1.1       2.0
            //    2.9       3.0
            //    3.0       3.0
            //    3.14159   4.0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the number whose ceiling is desired.
            //
            //    Output, double R8_CEILING, the ceiling of X.
            //
        {
            double value;

            value = (double) ((int) x);

            if (value < x)
            {
                value = value + 1.0;
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

        public static void r8_to_dhms(double r, ref int d, ref int h, ref int m, ref int s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_TO_DHMS converts an R8 day value into days, hours, minutes, seconds.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, a real number representing a time period measured in days.
            //
            //    Output, int &D, &H, &M, &S, the equivalent number of days, hours,
            //    minutes and seconds.
            //
        {
            int sign;

            if (0.0 <= r)
            {
                sign = 1;
            }
            else
            {
                sign = -1;
                r = -r;
            }

            d = (int) r;

            r = r - (double) d;
            r = 24.0 * r;
            h = (int) r;

            r = r - (double) h;
            r = 60.0 * r;
            m = (int) r;

            r = r - (double) m;
            r = 60.0 * r;
            s = (int) r;

            if (sign == -1)
            {
                d = -d;
                h = -h;
                m = -m;
                s = -s;
            }
        }

        public static int r8_to_i4(double xmin, double xmax, double x, int ixmin, int ixmax)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_TO_I4 maps real X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
            //
            //  Discussion:
            //
            //    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
            //    IX := min ( IX, max ( IXMIN, IXMAX ) )
            //    IX := max ( IX, min ( IXMIN, IXMAX ) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double XMIN, XMAX, the real range.  XMAX and XMIN must not be
            //    equal.  It is not necessary that XMIN be less than XMAX.
            //
            //    Input, double X, the real number to be converted.
            //
            //    Input, int IXMIN, IXMAX, the allowed range of the output
            //    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
            //    It is not necessary that IXMIN be less than IXMAX.
            //
            //    Output, int R8_TO_I4, the value in the range [IXMIN,IXMAX] that
            //    corresponds to X.
            //
        {
            int ix;
            double temp;

            if (xmax == xmin)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_TO_I4 - Fatal error!");
                Console.WriteLine("  XMAX = XMIN, making a zero divisor.");
                Console.WriteLine("  XMAX = " + xmax + "");
                Console.WriteLine("  XMIN = " + xmin + "");
                return (1);
            }

            temp =
                ((xmax - x) * (double) ixmin
                 + (x - xmin) * (double) ixmax)
                / (xmax - xmin);

            if (0.0 <= temp)
            {
                temp = temp + 0.5;
            }
            else
            {
                temp = temp - 0.5;
            }

            ix = (int) temp;

            return ix;
        }

        public static double r8_sum(double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_SUM returns the sum of two R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, Y, the quantities to add.
            //
            //    Output, double R8_SUM, the sum of X and Y.
            //
        {
            double value;

            value = x + y;

            return value;
        }

        public static void r8_swap(ref double x, ref double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_SWAP switches two R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, double &X, &Y.  On output, the values of X and
            //    Y have been interchanged.
            //
        {
            double z;

            z = x;
            x = y;
            y = z;
        }

        public static void r8_swap3(ref double x, ref double y, ref double z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_SWAP3 swaps three R8's.
            //
            //  Example:
            //
            //    Input:
            //
            //      X = 1, Y = 2, Z = 3
            //
            //    Output:
            //
            //      X = 2, Y = 3, Z = 1
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
            //    Input/output, double &X, &Y, &Z, three values to be swapped.
            //
        {
            double w;

            w = x;
            x = y;
            y = z;
            z = w;
        }

        public static double r8_tand(double degrees)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_TAND returns the tangent of an angle given in degrees.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double DEGREES, the angle in degrees.
            //
            //    Output, double R8_TAND, the tangent of the angle.
            //
        {
            const double r8_pi = 3.141592653589793;
            double radians;
            double value;

            radians = r8_pi * (degrees / 180.0);

            value = Math.Sin(radians) / Math.Cos(radians);

            return value;
        }

        public static double r8_tiny()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_TINY returns a "tiny" R8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double R8_TINY, a "tiny" R8 value.
            //
        {
            const double value = 0.4450147717014E-307;

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

        public static double r8_to_r8_discrete(double r, double rmin, double rmax, int nr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_TO_R8_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.
            //
            //  Discussion:
            //
            //    if ( R < RMIN ) then
            //      RD = RMIN
            //    else if ( RMAX < R ) then
            //      RD = RMAX
            //    else
            //      T = nint ( ( NR - 1 ) * ( R - RMIN ) / ( RMAX - RMIN ) )
            //      RD = RMIN + T * ( RMAX - RMIN ) / real ( NR - 1 )
            //
            //    In the special case where NR = 1, when
            //
            //      XD = 0.5 * ( RMAX + RMIN )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the number to be converted.
            //
            //    Input, double RMAX, RMIN, the maximum and minimum
            //    values for RD.
            //
            //    Input, int NR, the number of allowed values for XD.
            //    NR should be at least 1.
            //
            //    Output, double RD, the corresponding discrete value.
            //
        {
            int f;
            double rd;
            //
            //  Check for errors.
            //
            if (nr < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_TO_R8_DISCRETE - Fatal error!");
                Console.WriteLine("  NR = " + nr + "");
                Console.WriteLine("  but NR must be at least 1.");
                return (1);
            }

            if (nr == 1)
            {
                rd = 0.5 * (rmin + rmax);
                return rd;
            }

            if (rmax == rmin)
            {
                rd = rmax;
                return rd;
            }

            f = (int)((double) (nr) * (rmax - r) / (rmax - rmin));
            f = Math.Max(f, 0);
            f = Math.Min(f, nr);

            rd = ((double) (f) * rmin
                  + (double) (nr - f) * rmax)
                 / (double) (nr);

            return rd;
        }

        public static void r8_unswap3(ref double x, ref double y, ref double z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_UNSWAP3 unswaps three R8's.
            //
            //  Example:
            //
            //    Input:
            //
            //      X = 2, Y = 3, Z = 1
            //
            //    Output:
            //
            //      X = 1, Y = 2, Z = 3
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, double &X, &Y, &Z, three values to be swapped.
            //
        {
            double w;

            w = z;
            z = y;
            y = x;
            x = w;
        }


        public static void r8slmat_print(int m, int n, double[] a, string title)

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

        public static void r8_mant(double x, ref int s, ref double r, ref int l)

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
            if (x < 0.0)
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
            if (x < 0.0)
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
            if (x == 0.0)
            {
                return;
            }

            while (2.0 <= r)
            {
                r = r / 2.0;
                l = l + 1;
            }

            while (r < 1.0)
            {
                r = r * 2.0;
                l = l - 1;
            }

            return;
        }

        public static double r8_fractional(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FRACTIONAL returns the fractional part of an R8.
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
            //    Input, double X, the argument.
            //
            //    Output, double R8_FRACTIONAL, the fractional part of X.
            //
        {
            double value;

            value = Math.Abs(x) - (double) ((int) Math.Abs(x));

            return value;
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


        public static double r8_walsh_1d ( double x, int digit )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_WALSH_1D evaluates the Walsh function of a real scalar argument.
            //
            //  Discussion:
            //
            //    Consider the binary representation of X, and number the digits
            //    in descending order, from leading to lowest, with the units digit
            //    being numbered 0.
            //
            //    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument of the Walsh function.
            //
            //    Input, int DIGIT, the index of the Walsh function.
            //
            //    Output, double R8_WALSH_1D, the value of the Walsh function.
            //
        {
            int n;
            double value;
            //
            //  Hide the effect of the sign of X.
            //
            x = Math.Abs ( x );
            //
            //  If DIGIT is positive, divide by 2 DIGIT times.
            //  If DIGIT is negative, multiply by 2 (-DIGIT) times.
            //
            x = x / Math.Pow ( 2.0, digit );
            //
            //  Make it an integer.
            //  Because it's positive, and we're using INT, we don't change the
            //  units digit.
            //
            n = ( int ) x;
            //
            //  Is the units digit odd or even?
            //
            if ( ( n % 2 ) == 0 )
            {
                value = 0.0;
            }
            else
            {
                value = 1.0;
            }

            return value;
        }

        public static double r8_wrap(double r, double rlo, double rhi)

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
            if (rlo <= rhi)
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
            if (rwide == 0.0)
            {
                value = rlo;
            }
            else if (r < rlo2)
            {
                n = (int) ((rlo2 - r) / rwide) + 1;
                value = r + n * rwide;
                if (value == rhi)
                {
                    value = rlo;
                }
            }
            else
            {
                n = (int) ((r - rlo2) / rwide);
                value = r - n * rwide;
                if (value == rlo)
                {
                    value = rhi;
                }
            }

            return value;
        }
        
        public static double r8_atan(double y, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_ATAN computes the inverse tangent of the ratio Y / X.
            //
            //  Discussion:
            //
            //    R8_ATAN returns an angle whose tangent is ( Y / X ), a job which
            //    the built in functions ATAN and ATAN2 already do.
            //
            //    However:
            //
            //    * R8_ATAN always returns a positive angle, between 0 and 2 PI,
            //      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
            //      and [-PI,+PI] respectively;
            //
            //    * R8_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
            //     function by contrast always returns an angle in the first or fourth
            //     quadrants.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Y, X, two quantities which represent the tangent of
            //    an angle.  If Y is not zero, then the tangent is (Y/X).
            //
            //    Output, double R8_ATAN, an angle between 0 and 2 * PI, whose tangent is
            //    (Y/X), and which lies in the appropriate quadrant so that the signs
            //    of its cosine and sine match those of X and Y.
            //
        {
            double abs_x;
            double abs_y;
            double theta = 0;
            double theta_0;
            //
            //  Special cases:
            //
            if (x == 0.0)
            {
                if (0.0 < y)
                {
                    theta = Math.PI / 2.0;
                }
                else if (y < 0.0)
                {
                    theta = 3.0 * Math.PI / 2.0;
                }
                else if (y == 0.0)
                {
                    theta = 0.0;
                }
            }
            else if (y == 0.0)
            {
                if (0.0 < x)
                {
                    theta = 0.0;
                }
                else if (x < 0.0)
                {
                    theta = Math.PI;
                }
            }
            //
            //  We assume that ATAN2 is correct when both arguments are positive.
            //
            else
            {
                abs_y = Math.Abs(y);
                abs_x = Math.Abs(x);

                theta_0 = Math.Atan2(abs_y, abs_x);

                if (0.0 < x && 0.0 < y)
                {
                    theta = theta_0;
                }
                else if (x < 0.0 && 0.0 < y)
                {
                    theta = Math.PI - theta_0;
                }
                else if (x < 0.0 && y < 0.0)
                {
                    theta = Math.PI + theta_0;
                }
                else if (0.0 < x && y < 0.0)
                {
                    theta = 2.0 * Math.PI - theta_0;
                }
            }

            return theta;
        }


        
        public static double r8_nth_root ( double x, int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_NTH_ROOT returns the nth-root of an R8.
            //
            //  Discussion:
            //
            //    The nth root of X is x^(1/n)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 August 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, real X, the number whose nth root is desired.
            //
            //    Input, integer N, the index of the root.
            //
            //    Output, real VALUE, the Nth root of X.
            //
        {
            double e;
            double value;
            //
            //  Potential Error 1: 0^0
            //  But we will use it as 1.
            //
            if ( x == 0.0 && n == 0 )
            {
                value = 1.0;
                return value;
            }
            //
            //  Error 2: 0^(negative power)
            //
            if ( x == 0.0 && n < 0 )
            {
                value = double.NaN;
                return value;
            }
            //
            //  Error 3: (negative)^(even strictly positive root)
            //
            if ( x < 0.0 && ( n % 2 ) == 0 && 0 < n )
            {
                value = double.NaN;
                return value;
            }
            //
            //  X^0 = 1
            //
            if ( n == 0 )
            {
                value = 1.0;
            }
            //
            //  X^1 = X
            //
            else if ( n == 1 )
            {
                value = x;
            }
            //
            //  X^(-1) = 1/X
            //
            else if ( n == -1 )
            {
                value = 1.0 / x;
            }
            else
            {
                e = 1.0 / ( double ) ( Math.Abs ( n ) );

                if ( 0.0 < x )
                {
                    value = Math.Pow ( x, e );
                }
                else if ( x == 0.0 )
                {
                    value = 0.0;
                }
                else
                {
                    value = - Math.Pow ( - x, e );
                }

                if ( n < 0 )
                {
                    value = 1.0 / value;
                }
            }

            return value;
        }

        public static double r8_pi ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_PI returns the value of PI as an R8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double R8_PI, the value of PI.
            //
        {
            return Math.PI;
        }

        public static double r8_pi_sqrt ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_PI_SQRT returns the square root of PI as an R8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double R8_PI_SQRT, the square root of PI.
            //
        {
            return Math.Sqrt(Math.PI);
        }


        public static void r8_print(double r, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_PRINT prints an R8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the value to print.
            //
            //    Input, string TITLE, a title.
            //
        {
            Console.WriteLine(title + "  " + r + "");
        }

        public static double r8_random(int[] iseed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_RANDOM returns a uniformly distributed random number between 0 and 1.
            //
            //  Discussion:
            //
            //    This routine uses a multiplicative congruential method with modulus
            //    2**48 and multiplier 33952834046453 (see G.S.Fishman,
            //    'Multiplicative congruential random number generators with modulus
            //    2**b: an exhaustive analysis for b = 32 and a partial analysis for
            //    b = 48', Math. Comp. 189, pp 331-344, 1990).
            //
            //    48-bit integers are stored in 4 integer array elements with 12 bits
            //    per element. Hence the routine is portable across machines with
            //    integers of 32 bits or more.
            //
            //  Parameters:
            //
            //    Input/output, integer ISEED(4).
            //    On entry, the seed of the random number generator; the array
            //    elements must be between 0 and 4095, and ISEED(4) must be odd.
            //    On exit, the seed is updated.
            //
            //    Output, double R8_RANDOM, the next pseudorandom number.
            //
        {
            int ipw2 = 4096;
            int it1;
            int it2;
            int it3;
            int it4;
            int m1 = 494;
            int m2 = 322;
            int m3 = 2508;
            int m4 = 2549;
            double r = 1.0 / 4096.0;
            double value;
            //
            //  Multiply the seed by the multiplier modulo 2^48.
            //
            it4 = iseed[3] * m4;
            it3 = it4 / ipw2;
            it4 = it4 - ipw2 * it3;
            it3 = it3 + iseed[2] * m4 + iseed[3] * m3;
            it2 = it3 / ipw2;
            it3 = it3 - ipw2 * it2;
            it2 = it2 + iseed[1] * m4 + iseed[2] * m3 + iseed[3] * m2;
            it1 = it2 / ipw2;
            it2 = it2 - ipw2 * it1;
            it1 = it1 + iseed[0] * m4 + iseed[1] * m3 + iseed[2] * m2 + iseed[3] * m1;
            it1 = (it1 % ipw2);
            //
            //  Return updated seed
            //
            iseed[0] = it1;
            iseed[1] = it2;
            iseed[2] = it3;
            iseed[3] = it4;
            //
            //  Convert 48-bit integer to a real number in the interval (0,1)
            //
            value =
                r * ((double) (it1)
                     + r * ((double) (it2)
                            + r * ((double) (it3)
                                   + r * ((double) (it4)))));

            return value;
        }
    }
}