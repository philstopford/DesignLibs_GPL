using System;
using System.Numerics;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static Complex cartesian_to_c8(double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CARTESIAN_TO_C8 converts a Cartesian form to a C8.
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
            //    Input, double X, Y, the Cartesian form.
            //
            //    Output, Complex CARTESIAN_TO_C8, the complex number.
            //
        {
            Complex c;

            c = new Complex(x, y);

            return c;
        }

        public static Complex polar_to_c8(double r, double theta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLAR_TO_C8 converts a polar form to a C8.
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
            //    Input, double R, THETA, the polar form.
            //
            //    Output, Complex POLAR_TO_C8, the complex number.
            //
        {
            Complex c;

            c = new Complex(r * Math.Cos(theta), r * Math.Sin(theta));

            return c;
        }

        public static Complex r8_csqrt(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CSQRT returns the complex square root of an R8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the number whose square root is desired.
            //
            //    Output, Complex R8_CSQRT, the square root of X:
            //
        {
            double argument = 0;
            double magnitude = 1;
            const double r8_pi = 3.141592653589793;
            Complex value;

            if (0.0 < x)
            {
                magnitude = x;
                argument = 0.0;
            }
            else if (0.0 == x)
            {
                magnitude = 0.0;
                argument = 0.0;
            }
            else if (x < 0.0)
            {
                magnitude = -x;
                argument = r8_pi;
            }

            magnitude = Math.Sqrt(magnitude);
            argument = argument / 2.0;

            value = magnitude * new Complex(Math.Cos(argument), Math.Sin(argument));

            return value;
        }

        public static void r8poly2_root(double a, double b, double c, ref Complex r1,
        ref Complex r2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY2_ROOT returns the two roots of a quadratic polynomial.
        //
        //  Discussion:
        //
        //    The polynomial has the form:
        //
        //      A * X^2 + B * X + C = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 October 2005
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the coefficients of the polynomial.
        //    A must not be zero.
        //
        //    Output, Complex &R1, &R2, the roots of the polynomial, which
        //    might be real and distinct, real and equal, or complex conjugates.
        //
        {
            double disc;
            Complex q;

            if (a == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8POLY2_ROOT - Fatal error!");
                Console.WriteLine("  The coefficient A is zero.");
                return;
            }

            disc = b * b - 4.0 * a * c;
            q = -0.5 * (b + r8_sign(b) * r8_csqrt(disc));
            r1 = q / a;
            r2 = c / q;
        }

        public static void r8poly3_root(double a, double b, double c, double d,
            ref Complex r1, ref Complex r2, ref Complex r3 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY3_ROOT returns the three roots of a cubic polynomial.
        //
        //  Discussion:
        //
        //    The polynomial has the form
        //
        //      A * X^3 + B * X^2 + C * X + D = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 October 2005
        //
        //  Parameters:
        //
        //    Input, double A, B, C, D, the coefficients of the polynomial.
        //    A must not be zero.
        //
        //    Output, Complex &R1, &R2, &R3, the roots of the polynomial, which
        //    will include at least one real root.
        //
        {
            Complex i;
            const double r8_pi = 3.141592653589793;
            double q;
            double r;
            double s1;
            double s2;
            double temp;
            double theta;

            if (a == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8POLY3_ROOT - Fatal error!");
                Console.WriteLine("  A must not be zero.");
                return;
            }

            i = new Complex(0.0, 1.0);

            q = (Math.Pow(b / a, 2) - 3.0 * (c / a)) / 9.0;

            r = (2.0 * Math.Pow(b / a, 3) - 9.0 * (b / a) * (c / a)
                 + 27.0 * (d / a)) / 54.0;

            if (r * r < q * q * q)
            {
                theta = Math.Acos(r / Math.Sqrt(Math.Pow(q, 3)));
                r1 = -2.0 * Math.Sqrt(q) * Math.Cos(theta / 3.0);
                r2 = -2.0 * Math.Sqrt(q) * Math.Cos((theta + 2.0 * r8_pi) / 3.0);
                r3 = -2.0 * Math.Sqrt(q) * Math.Cos((theta + 4.0 * r8_pi) / 3.0);
            }
            else if (q * q * q <= r * r)
            {
                temp = -r + Math.Sqrt(r * r - q * q * q);
                s1 = r8_sign(temp) * Math.Pow(Math.Abs(temp), 1.0 / 3.0);

                temp = -r - Math.Sqrt(r * r - q * q * q);
                s2 = r8_sign(temp) * Math.Pow(Math.Abs(temp), 1.0 / 3.0);

                r1 = s1 + s2;
                r2 = -0.5 * (s1 + s2) + i * 0.5 * Math.Sqrt(3.0) * (s1 - s2);
                r3 = -0.5 * (s1 + s2) - i * 0.5 * Math.Sqrt(3.0) * (s1 - s2);
            }

            r1 = r1 - b / (3.0 * a);
            r2 = r2 - b / (3.0 * a);
            r3 = r3 - b / (3.0 * a);
        }

        public static void r8poly4_root(double a, double b, double c, double d, double e,
            ref Complex r1, ref Complex r2, ref Complex r3,
        ref Complex r4 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY4_ROOT returns the four roots of a quartic polynomial.
        //
        //  Discussion:
        //
        //    The polynomial has the form:
        //
        //      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2005
        //
        //  Parameters:
        //
        //    Input, double A, B, C, D, the coefficients of the polynomial.
        //    A must not be zero.
        //
        //    Output, Complex &R1, &R2, &R3, &R4, the roots of the polynomial.
        //
        {
            double a3;
            double a4;
            double b3;
            double b4;
            double c3;
            double c4;
            double d3;
            double d4;
            Complex p;
            Complex q;
            Complex r;
            Complex zero;

            zero = 0.0;

            if (a == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8POLY4_ROOT - Fatal error!");
                Console.WriteLine("  A must not be zero.");
                return;
            }

            a4 = b / a;
            b4 = c / a;
            c4 = d / a;
            d4 = e / a;
            //
            //  Set the coefficients of the resolvent cubic equation.
            //
            a3 = 1.0;
            b3 = -b4;
            c3 = a4 * c4 - 4.0 * d4;
            d3 = -a4 * a4 * d4 + 4.0 * b4 * d4 - c4 * c4;
            //
            //  Find the roots of the resolvent cubic.
            //
            r8poly3_root(a3, b3, c3, d3, ref r1, ref r2, ref r3);
            //
            //  Choose one root of the cubic, here R1.
            //
            //  Set R = Math.Sqrt ( 0.25 * A4^2 - B4 + R1 )
            //
            r = c8_sqrt(0.25 * a4 * a4 - b4 + r1);

            if ((r.Real) != 0.0 || (r.Imaginary) != 0.0)
            {
                p = c8_sqrt(0.75 * a4 * a4 - r * r - 2.0 * b4
                            + 0.25 * (4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4) / r);

                q = c8_sqrt(0.75 * a4 * a4 - r * r - 2.0 * b4
                            - 0.25 * (4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4) / r);
            }
            else
            {
                p = c8_sqrt(0.75 * a4 * a4 - 2.0 * b4
                            + 2.0 * c8_sqrt(r1 * r1 - 4.0 * d4));

                q = c8_sqrt(0.75 * a4 * a4 - 2.0 * b4
                                           - 2.0 * c8_sqrt(r1 * r1 - 4.0 * d4));
            }

            //
            //  Set the roots.
            //
            r1 = -0.25 * a4 + 0.5 * r + 0.5 * p;
            r2 = -0.25 * a4 + 0.5 * r - 0.5 * p;
            r3 = -0.25 * a4 - 0.5 * r + 0.5 * q;
            r4 = -0.25 * a4 - 0.5 * r - 0.5 * q;
        }
    }
}