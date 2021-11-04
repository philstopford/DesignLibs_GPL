using System;
using System.Numerics;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double c8_abs(Complex x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ABS returns the absolute value of a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex X, the value whose norm is desired.
            //
            //    Output, double C8_ABS, the magnitude of X.
            //
        {
            double value;

            value = Math.Sqrt(Math.Pow((x.Real), 2) + Math.Pow((x.Imaginary), 2));

            return value;
        }

        public static Complex c8_acos(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ACOS evaluates the inverse cosine of a C8.
            //
            //  Discussion:
            //
            //    Here we use the relationship:
            //
            //      C8_ACOS ( Z ) = pi/2 - C8_ASIN ( Z ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_ACOS, the function value.
            //
        {
            Complex c2;
            double c2_imag;
            double c2_real;

            c2 = c8_asin(c1);

            c2_real = (Math.PI * 0.5) - (c2.Real);
            c2_imag = -(c2.Imaginary);

            c2 = new Complex(c2_real, c2_imag);

            return c2;
        }

        public static Complex c8_acosh(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ACOSH evaluates the inverse hyperbolic cosine of a C8.
            //
            //  Discussion:
            //
            //    Here we use the relationship:
            //
            //      C8_ACOSH ( Z ) = i * C8_ACOS ( Z ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_ACOSH, the function value.
            //
        {
            Complex c2;

            c2 = c8_i() * c8_acos(c1);

            return c2;
        }

        public static Complex c8_add(Complex c1, Complex c2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ADD adds two C8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 March 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex C1, C2, the arguments.
            //
            //    Output, Complex C8_ADD, the sum of C1 and C2.
            //
        {
            Complex c3;

            c3 = c1 + c2;

            return c3;
        }

        public static double c8_arg(Complex x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ARG returns the argument of a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
            //
            //    The value returned by this function is always between 0 and 2*PI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex X, the value whose argument is desired.
            //
            //    Output, double C8_ARG, the argument of X.
            //
        {
            double value;

            if ((x.Imaginary) == 0.0 && (x.Real) == 0.0)
            {
                value = 0.0;
            }
            else
            {
                value = r8_atan((x.Imaginary), (x.Real));
            }

            return value;
        }

        public static Complex c8_asin(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ASIN evaluates the inverse sine of a C8.
            //
            //  Discussion:
            //
            //    Here we use the relationship:
            //
            //      C8_ASIN ( Z ) = - i * log ( i * z + Math.Sqrt ( 1 - z * z ) )
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_ASIN, the function value.
            //
        {
            Complex c2;
            Complex c3;
            Complex c4;
            Complex ce;

            c2 = c8_i();
            c3 = Complex.Sqrt(1.0 - c1 * c1);
            c4 = c8_log(c3 + c2 * c1);
            ce = -c2 * c4;

            return ce;
        }

        public static Complex c8_asinh(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ASINH evaluates the inverse hyperbolic sine of a C8.
            //
            //  Discussion:
            //
            //    Here we use the relationship:
            //
            //      C8_ASINH ( Z ) = - i * C8_ASIN ( i * Z ).
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_ASINH, the function value.
            //
        {
            Complex c2;
            Complex c3;
            Complex c4;
            Complex c5;
            Complex c6;

            c2 = c8_i();
            c3 = c2 * c1;
            c4 = c8_asin(c3);
            c5 = c2 * c4;
            c6 = -c5;

            return c6;
        }

        public static Complex c8_atan(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ATAN evaluates the inverse tangent of a C8.
            //
            //  Discussion:
            //
            //    Here we use the relationship:
            //
            //      C8_ATAN ( Z ) = ( i / 2 ) * log ( ( 1 - i * z ) / ( 1 + i * z ) )
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_ATAN, the function value.
            //
        {
            Complex c2;
            Complex c3;
            Complex c4;
            Complex c5;
            Complex c6;
            Complex c7;
            Complex c8;
            Complex c9;
            Complex cx;

            c2 = c8_i();
            c3 = c8_one();
            c4 = c8_mul(c2, c1);
            c5 = c8_sub(c3, c4);
            c6 = c8_add(c3, c4);
            c7 = c8_div(c5, c6);

            c8 = c8_log(c7);
            c9 = c8_mul(c2, c8);
            cx = c9 / 2.0;

            return cx;
        }

        public static Complex c8_atanh(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ATANH evaluates the inverse hyperbolic tangent of a C8.
            //
            //  Discussion:
            //
            //    Here we use the relationship:
            //
            //      C8_ATANH ( Z ) = - i * C8_ATAN ( i * Z ).
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_ATANH, the function value.
            //
        {
            Complex c2;
            Complex c3;
            Complex c4;
            Complex c5;
            Complex c6;

            c2 = c8_i();

            c3 = c8_mul(c2, c1);
            c4 = c8_atan(c3);
            c5 = c8_mul(c2, c4);
            c6 = c8_neg(c5);

            return c6;
        }

        public static Complex c8_conj(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_CONJ conjugates a C8.
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_CONJ, the function value.
            //
        {
            Complex c2;

            c2 = Complex.Conjugate(c1);

            return c2;
        }

        public static void c8_copy(ref Complex c1, Complex c2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_COPY copies a C8.
            //
            //  Discussion:
            //
            //    The order of the arguments may seem unnatural, but it is arranged so
            //    that the call
            //
            //      c8_copy ( c1, c2 )
            //
            //    mimics the assignment
            //
            //      c1 = c2.
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
            //    Output, Complex C1, the copy of C2.
            //
            //    Input, Complex C2, the value to be copied.
            //
        {
            c1 = c2;

        }

        public static Complex c8_cos(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_COS evaluates the cosine of a C8.
            //
            //  Discussion:
            //
            //    We use the relationship:
            //
            //      C8_COS ( C ) = ( C8_EXP ( i * C ) + C8_EXP ( - i * C ) ) / 2
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_COS, the function value.
            //
        {
            Complex c2;

            c2 = (Complex.Exp(c1 * c8_i()) + Complex.Exp(-c1 * c8_i())) / 2.0;

            return c2;
        }

        public static Complex c8_cosh(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_COSH evaluates the hyperbolic cosine of a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_COSH, the function value.
            //
        {
            Complex c2;
            Complex c3;
            Complex c4;
            Complex c5;
            Complex c6;

            c2 = c8_exp(c1);

            c3 = c8_neg(c1);
            c4 = c8_exp(c3);

            c5 = c8_add(c2, c4);
            c6 = c8_div_r8(c5, 2.0);

            return c6;
        }

        public static Complex c8_cube_root(Complex x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_CUBE_ROOT returns the principal cube root of a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
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
            //    Input, Complex X, the number whose cube root is desired.
            //
            //    Output, Complex C8_CUBE_ROOT, the cube root of X.
            //
        {
            double argument;
            double magnitude;
            Complex value;

            argument = c8_arg(x);
            magnitude = c8_mag(x);

            if (magnitude == 0.0)
            {
                value = new Complex(0.0, 0.0);
            }
            else
            {
                value = Math.Pow(magnitude, (double) (1.0 / 3.0))
                        * new Complex(Math.Cos(argument / 3.0), Math.Sin(argument / 3.0));
            }

            return value;
        }

        public static Complex c8_div(Complex c1, Complex c2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_DIV divides two C8's.
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
            //    Input, Complex C1, C2, the arguments.
            //
            //    Output, Complex C8_DIV, the function value.
            //
        {
            double c2_norm;
            Complex c3;
            double c3_imag;
            double c3_real;

            c2_norm = c8_abs(c2);

            c3_real = ((c1.Real) * (c2.Real)
                       + (c1.Imaginary) * (c2.Imaginary)) / c2_norm / c2_norm;

            c3_imag = ((c1.Imaginary) * (c2.Real)
                       - (c1.Real) * (c2.Imaginary)) / c2_norm / c2_norm;

            c3 = new Complex(c3_real, c3_imag);

            return c3;
        }

        public static Complex c8_div_r8(Complex c1, double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_DIV_R8 divides a C8 by an R8.
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
            //    Input, Complex C1, the value to be divided.
            //
            //    Input, double R, the divisor.
            //
            //    Output, Complex C8_DIV_R8, the function value.
            //
        {
            Complex c2;

            c2 = c1 / r;

            return c2;
        }

        public static Complex c8_exp(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_EXP exponentiates a C8.
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_EXP, the function value.
            //
        {
            Complex c2;

            c2 = Complex.Exp(c1);

            return c2;
        }

        public static Complex c8_i()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_I returns the value of the imaginary unit, i as a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
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
            //    Output, Complex C8_I, the value of complex i.
            //
        {
            Complex value;

            value = new Complex(0.0, 1.0);

            return value;
        }

        public static double c8_imag(Complex c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_IMAG returns the imaginary part of a C8.
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
            //    Input, Complex C, the argument.
            //
            //    Output, double C8_IMAG, the function value.
            //
        {
            double value;

            value = (c.Imaginary);

            return value;
        }

        public static Complex c8_inv(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_INV inverts a C8.
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_INV, the function value;
            //
        {
            Complex c2;

            c2 = 1.0 / c1;

            return c2;
        }

        public static bool c8_le_l1(Complex x, Complex y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_LE_L1 := X <= Y for C8 values, and the L1 norm.
            //
            //  Discussion:
            //
            //    A C8 is a complex double precision value.
            //
            //    The L1 norm can be defined here as:
            //
            //      C8_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex X, Y, the values to be compared.
            //
            //    Output, bool C8_LE_L1, is TRUE if X <= Y.
            //
        {
            bool value;

            if (Math.Abs((x.Real)) + Math.Abs((x.Imaginary)) <=
                Math.Abs((y.Real)) + Math.Abs((y.Imaginary)))
            {
                value = true;
            }
            else
            {
                value = false;
            }

            return value;
        }

        public static bool c8_le_l2(Complex x, Complex y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_LE_L2 := X <= Y for C8 values, and the L2 norm.
            //
            //  Discussion:
            //
            //    A C8 is a complex double precision value.
            //
            //    The L2 norm can be defined here as:
            //
            //      C8_NORM_L2(X) = Math.Sqrt ( ( real (X) )^2 + ( imag (X) )^2 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex X, Y, the values to be compared.
            //
            //    Output, bool C8_LE_L2, is TRUE if X <= Y.
            //
        {
            bool value;

            if (Math.Pow((x.Real), 2) + Math.Pow((x.Imaginary), 2) <=
                Math.Pow((y.Real), 2) + Math.Pow((y.Imaginary), 2))
            {
                value = true;
            }
            else
            {
                value = false;
            }

            return value;
        }

        public static bool c8_le_li(Complex x, Complex y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_LE_LI := X <= Y for C8 values, and the L-oo norm.
            //
            //  Discussion:
            //
            //    A C8 is a complex double precision value.
            //
            //    The L-oo norm can be defined here as:
            //
            //      C8_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex X, Y, the values to be compared.
            //
            //    Output, bool C8_LE_LI, is TRUE if X <= Y.
            //
        {
            bool value;

            if (Math.Max(Math.Abs((x.Real)), Math.Abs((x.Imaginary))) <=
                Math.Max(Math.Abs((y.Real)), Math.Abs((y.Imaginary))))
            {
                value = true;
            }
            else
            {
                value = false;
            }

            return value;
        }

        public static Complex c8_log(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_LOG evaluates the logarithm of a C8.
            //
            //  Discussion:
            //
            //    Here we use the relationship:
            //
            //      C8_LOG ( Z ) = LOG ( MAG ( Z ) ) + i * ARG ( Z )
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_LOG, the function value.
            //
        {
            double arg;
            Complex c2;
            double mag;

            arg = c8_arg(c1);
            mag = c8_mag(c1);

            c2 = new Complex(Math.Log(mag), arg);

            return c2;
        }

        public static double c8_mag(Complex x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_MAG returns the magnitude of a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex X, the value whose norm is desired.
            //
            //    Output, double C8_MAG, the magnitude of X.
            //
        {
            double magnitude;

            magnitude = Math.Sqrt(Math.Pow((x.Real), 2) + Math.Pow((x.Imaginary), 2));

            return magnitude;
        }

        public static Complex c8_mul(Complex c1, Complex c2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_MUL multiplies two C8's.
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
            //    Input, Complex C1, C2, the arguments.
            //
            //    Output, Complex C8_MUL, the function value.
            //
        {
            Complex c3;

            c3 = c1 * c2;

            return c3;
        }

        public static Complex c8_neg(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_NEG negates a C8.
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_NEG, the function value.
            //
        {
            Complex c2;

            c2 = -c1;

            return c2;
        }

        public static Complex c8_nint(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_NINT returns the nearest complex integer of a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex C1, the value to be NINT'ed.
            //
            //    Output, Complex C8_NINT, the NINT'ed value.
            //
        {
            double r;
            double r_min;
            double x;
            double x_min;
            double xc;
            double y;
            double y_min;
            double yc;
            Complex value;

            xc = (c1.Real);
            yc = (c1.Imaginary);
            //
            //  Lower left.
            //
            x = Math.Floor((c1.Real));
            y = Math.Floor((c1.Imaginary));
            r = Math.Pow(x - xc, 2) + Math.Pow(y - yc, 2);
            r_min = r;
            x_min = x;
            y_min = y;
            //
            //  Lower right.
            //
            x = Math.Floor((c1.Real)) + 1.0;
            y = Math.Floor((c1.Imaginary));
            r = Math.Pow(x - xc, 2) + Math.Pow(y - yc, 2);
            if (r < r_min)
            {
                r_min = r;
                x_min = x;
                y_min = y;
            }

            //
            //  Upper right.
            //
            x = Math.Floor((c1.Real)) + 1.0;
            y = Math.Floor((c1.Imaginary)) + 1.0;
            r = Math.Pow(x - xc, 2) + Math.Pow(y - yc, 2);
            if (r < r_min)
            {
                r_min = r;
                x_min = x;
                y_min = y;
            }

            //
            //  Upper left.
            //
            x = Math.Floor((c1.Real));
            y = Math.Floor((c1.Imaginary)) + 1.0;
            r = Math.Pow(x - xc, 2) + Math.Pow(y - yc, 2);
            if (r < r_min)
            {
                r_min = r;
                x_min = x;
                y_min = y;
            }

            value = new Complex(x_min, y_min);

            return value;
        }

        public static double c8_norm_l1(Complex x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_NORM_L1 evaluates the L1 norm of a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
            //
            //    Numbers of equal norm lie along diamonds centered at (0,0).
            //
            //    The L1 norm can be defined here as:
            //
            //      C8_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex X, the value whose norm is desired.
            //
            //    Output, double C8_NORM_L1, the norm of X.
            //
        {
            double value;

            value = Math.Abs((x.Real)) + Math.Abs((x.Imaginary));

            return value;
        }

        public static double c8_norm_l2(Complex x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_NORM_L2 evaluates the L2 norm of a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
            //
            //    Numbers of equal norm lie on circles centered at (0,0).
            //
            //    The L2 norm can be defined here as:
            //
            //      C8_NORM_L2(X) = Math.Sqrt ( ( real (X) )^2 + ( imag ( X ) )^2 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex X, the value whose norm is desired.
            //
            //    Output, double C8_NORM_L2, the 2-norm of X.
            //
        {
            double value;

            value = Math.Sqrt(Math.Pow((x.Real), 2)
                              + Math.Pow((x.Imaginary), 2));

            return value;
        }

        public static double c8_norm_li(Complex x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_NORM_LI evaluates the L-oo norm of a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
            //
            //    Numbers of equal norm lie along squares whose centers are at (0,0).
            //
            //    The L-oo norm can be defined here as:
            //
            //      C8_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, Complex X, the value whose norm is desired.
            //
            //    Output, double C8_NORM_LI, the L-oo norm of X.
            //
        {
            double value;

            value = Math.Max(Math.Abs((x.Real)), Math.Abs((x.Imaginary)));

            return value;
        }

        public static Complex c8_normal_01(ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_NORMAL_01 returns a unit pseudonormal C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
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
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, Complex C8_NORMAL_01, a unit pseudornormal value.
            //
        {
            
            double v1;
            double v2;
            double x_c;
            double x_r;
            Complex value;

            v1 = UniformRNG.r8_uniform_01(ref seed);
            v2 = UniformRNG.r8_uniform_01(ref seed);

            x_r = Math.Sqrt(-2.0 * Math.Log(v1)) * Math.Cos(2.0 * Math.PI * v2);
            x_c = Math.Sqrt(-2.0 * Math.Log(v1)) * Math.Sin(2.0 * Math.PI * v2);

            value = new Complex(x_r, x_c);

            return value;
        }

        public static Complex c8_one()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ONE returns the value of complex 1.
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
            //    Output, Complex C8_ONE, the value of complex 1.
            //
        {
            Complex value;

            value = new Complex(1.0, 0.0);

            return value;
        }

        public static void c8_print(Complex a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_PRINT prints a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
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
            //    Input, Complex A, the value to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            Console.WriteLine(title
                              + "  ( " + (a.Real.ToString().PadLeft(14))
                              + ", " + (a.Imaginary.ToString().PadLeft(14)) + " )");

            return;
        }

        public static double c8_real(Complex c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_REAL returns the real part of a C8.
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
            //    Input, Complex C, the complex number.
            //
            //    Output, double C8_REAL, the function value.
            //
        {
            double value;

            value = (c.Real);

            return value;
        }

        public static Complex c8_sin(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_SIN evaluates the sine of a C8.
            //
            //  Discussion:
            //
            //    We use the relationship:
            //
            //      C8_SIN ( C ) = - i * ( C8_EXP ( i * C ) - C8_EXP ( - i * C ) ) / 2
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_SIN, the function value.
            //
        {
            Complex c2;
            Complex c3;
            Complex c4;
            Complex c5;
            Complex c6;
            Complex c7;
            Complex c8;
            Complex c9;
            Complex cx;
            double r;

            c2 = c8_i();

            c3 = c8_mul(c2, c1);
            c4 = c8_exp(c3);

            c5 = c8_neg(c3);
            c6 = c8_exp(c5);

            c7 = c8_sub(c4, c6);

            r = 2.0;
            c8 = c8_div_r8(c7, r);
            c9 = c8_mul(c8, c2);
            cx = c8_neg(c9);

            return cx;
        }

        public static Complex c8_sinh(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_SINH evaluates the hyperbolic sine of a C8.
            //
            //  Discussion:
            //
            //    We use the relationship:
            //
            //      C8_SINH ( C ) = ( C8_EXP ( C ) - C8_EXP ( - i * C ) ) / 2
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_SINH, the function value.
            //
        {
            Complex c2;
            Complex c3;
            Complex c4;
            Complex c5;
            Complex c6;
            double r;

            c2 = c8_exp(c1);

            c3 = c8_neg(c1);
            c4 = c8_exp(c3);

            c5 = c8_sub(c2, c4);

            r = 2.0;
            c6 = c8_div_r8(c5, r);

            return c6;
        }

        public static Complex c8_sqrt(Complex x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_SQRT returns the principal square root of a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
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
            //    Input, Complex X, the number whose square root is desired.
            //
            //    Output, Complex C8_SQRT, the square root of X.
            //
        {
            double argument;
            double magnitude;
            Complex value;

            argument = c8_arg(x);
            magnitude = c8_mag(x);

            if (magnitude == 0.0)
            {
                value = new Complex(0.0, 0.0);
            }
            else
            {
                value = Math.Sqrt(magnitude)
                        * new Complex(Math.Cos(argument / 2.0), Math.Sin(argument / 2.0));
            }

            return value;
        }

        public static Complex c8_sub(Complex c1, Complex c2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_SUB subtracts two C8's.
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
            //    Input, Complex C1, C2, the arguments.
            //
            //    Output, Complex C8_SUB, the function value.
            //
        {
            Complex c3;

            c3 = c1 - c2;

            return c3;
        }

        public static void c8_swap(ref Complex x, ref Complex y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_SWAP swaps two C8's.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, Complex &X, &Y.  On output, the values of X and
            //    Y have been interchanged.
            //
        {
            Complex z;

            z = x;
            x = y;
            y = z;

            return;
        }

        public static Complex c8_tan(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_TAN evaluates the tangent of a C8.
            //
            //  Discussion:
            //
            //    We use the relationship:
            //
            //      C8_TAN ( C ) = - i * ( C8_EXP ( i * C ) - C8_EXP ( - i * C ) ) 
            //                         / ( C8_EXP ( I * C ) + C8_EXP ( - i * C ) )
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_TAN, the function value.
            //
        {
            Complex c2;
            Complex c3;
            Complex c4;
            Complex c5;
            Complex c6;
            Complex c7;
            Complex c8;
            Complex c9;
            Complex cx;
            Complex ce;

            c2 = c8_i();
            c3 = c8_mul(c2, c1);
            c4 = c8_neg(c3);

            c5 = c8_exp(c3);
            c6 = c8_exp(c4);

            c7 = c8_sub(c5, c6);
            c8 = c8_add(c5, c6);

            c9 = c8_div(c7, c8);
            cx = c8_mul(c2, c9);
            ce = c8_neg(cx);

            return ce;
        }

        public static Complex c8_tanh(Complex c1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_TANH evaluates the hyperbolic tangent of a C8.
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
            //    Input, Complex C1, the argument.
            //
            //    Output, Complex C8_TANH, the function value.
            //
        {
            Complex c2;
            Complex c3;
            Complex c4;
            Complex c5;
            Complex c6;
            Complex c7;

            c2 = c8_exp(c1);

            c3 = c8_neg(c1);
            c4 = c8_exp(c3);

            c5 = c8_sub(c2, c4);
            c6 = c8_add(c2, c4);

            c7 = c8_div(c5, c6);

            return c7;
        }

        public static void c8_to_cartesian(Complex c, ref double x, ref double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_TO_CARTESIAN converts a C8 to Cartesian form.
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
            //    Input, Complex C, the argument.
            //
            //    Output, double &X, &Y, the Cartesian form.
            //
        {
            x = (c.Real);
            y = (c.Imaginary);
        }

        public static void c8_to_polar(Complex c, ref double r, ref double theta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_TO_POLAR converts a C8 to polar form.
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
            //    Input, Complex C, the argument.
            //
            //    Output, double &R, &THETA, the polar form.
            //
        {
            r = c8_abs(c);
            theta = c8_arg(c);
        }

        public static Complex c8_uniform_01(ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_UNIFORM_01 returns a unit pseudorandom C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
            //
            //    The angle should be uniformly distributed between 0 and 2 * PI,
            //    the square root of the radius uniformly distributed between 0 and 1.
            //
            //    This results in a uniform distribution of values in the unit circle.
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
            //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
            //    On output, SEED has been updated.
            //
            //    Output, Complex C8_UNIFORM_01, a pseudorandom complex value.
            //
        {
            double r;
            int k;
            
            double theta;
            Complex value;

            k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            if (seed < 0)
            {
                seed = seed + 2147483647;
            }

            r = Math.Sqrt((double) (seed) * 4.656612875E-10);

            k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            if (seed < 0)
            {
                seed = seed + 2147483647;
            }

            theta = 2.0 * Math.PI * ((double) (seed) * 4.656612875E-10);

            value = r * new Complex(Math.Cos(theta), Math.Sin(theta));

            return value;
        }

        public static Complex c8_zero()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_ZERO returns the value of 0 as a C8.
            //
            //  Discussion:
            //
            //    A C8 is a Complex value.
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
            //    Output, Complex C8_ZERO, the value of complex 0.
            //
        {
            Complex value;

            value = new Complex(0.0, 0.0);

            return value;
        }

    }
}