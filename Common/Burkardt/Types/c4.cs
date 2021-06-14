using System;
using System.Numerics;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static float cabs1(Complex z)
//****************************************************************************80
//
//  Purpose:
//
//    CABS1 returns the L1 norm of a number.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L1 norm of a complex number is the sum of the absolute values
//    of the real and imaginary components.
//
//    CABS1 ( Z ) = abs ( real ( Z ) ) + abs ( imaginary ( Z ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2007
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, complex <float> Z, the number whose norm is desired.
//
//    Output, float CABS1, the L1 norm of Z.
//
        {
            float value = (float) (Math.Abs(z.Real) + Math.Abs(z.Imaginary));

            return value;
        }

        public static float cabs2(Complex z)

//****************************************************************************80
//
//  Purpose:
//
//    CABS2 returns the L2 norm of a number.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L2 norm of a complex number is the square root of the sum
//    of the squares of the real and imaginary components.
//
//    CABS2 ( Z ) = sqrt ( ( real ( Z ) )**2 + ( imaginary ( Z ) )**2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, complex <float> Z, the number whose norm is desired.
//
//    Output, float CABS2, the L2 norm of Z.
//
        {
            float value = (float) Math.Sqrt(Math.Pow(z.Real, 2)
                                            + Math.Pow(z.Imaginary, 2));

            return value;
        }

        public static float cmach(int job)

//****************************************************************************80
//
//  Purpose:
//
//    CMACH computes machine parameters for complex arithmetic.
//
//  Discussion:
//
//    Assume the computer has
//
//      B = base of arithmetic;
//      T = number of base B digits;
//      L = smallest possible exponent;
//      U = largest possible exponent;
//
//    then
//
//      EPS = B^(1-T)
//      TINY = 100.0 * B^(-L+T)
//      HUGE = 0.01 * B^(U-T)
//
//    If complex division is done by
//
//      1 / (X+i*Y) = (X-i*Y) / (X^2+Y^2)
//
//    then
//
//      TINY = sqrt ( TINY )
//      HUGE = sqrt ( HUGE )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, int JOB:
//    1, EPS is desired;
//    2, TINY is desired;
//    3, HUGE is desired.
//
//    Output, float CMACH, the requested value.
//
        {
            float eps;
            float huge;
            float s;
            Complex temp1;
            Complex temp2;
            Complex temp3;
            float tiny;
            float value;

            eps = 1.0f;

            for (;;)
            {
                eps = eps / 2.0f;
                s = 1.0f + eps;
                if (s <= 1.0)
                {
                    break;
                }
            }

            eps = 2.0f * eps;

            s = 1.0f;

            for (;;)
            {
                tiny = s;
                s = s / 16.0f;

                if (s * 1.0 == 0.0)
                {
                    break;
                }
            }

            tiny = (tiny / eps) * 100.0f;
//
//  Had to insert this manually!
//
            tiny = (float) Math.Sqrt(tiny);

            if (false)
            {
                temp1 = new Complex(1.0, 0.0);
                temp2 = new Complex(tiny, 0.0);
                temp3 = temp1 / temp2;

                s = (float) temp3.Real;

                if (s != 1.0 / tiny)
                {
                    tiny = (float) Math.Sqrt(tiny);
                }
            }

            huge = 1.0f / tiny;

            if (job == 1)
            {
                value = eps;
            }
            else if (job == 2)
            {
                value = tiny;
            }
            else if (job == 3)
            {
                value = huge;
            }
            else
            {
                value = 0.0f;
            }

            return value;
        }

        public static Complex csign1(Complex z1, Complex z2)

//****************************************************************************80
//
//  Purpose:
//
//    CSIGN1 is a transfer-of-sign function.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L1 norm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> Z1, Z2, the arguments.
//
//    Output, complex <float> CSIGN1,  a complex value, with the magnitude of
//    Z1, and the argument of Z2.
//
        {
            Complex value;

            if (cabs1(z2) == 0.0)
            {
                value = new Complex(0.0, 0.0);
            }
            else
            {
                value = cabs1(z1) * (z2 / cabs1(z2));
            }

            return value;
        }

        public static Complex csign2(Complex z1, Complex z2)

//****************************************************************************80
//
//  Purpose:
//
//    CSIGN2 is a transfer-of-sign function.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L2 norm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> Z1, Z2, the arguments.
//
//    Output, complex <float> CSIGN2,  a complex value, with the magnitude of
//    Z1, and the argument of Z2.
//
        {
            Complex value;

            if (cabs2(z2) == 0.0)
            {
                value = new Complex(0.0, 0.0);
            }
            else
            {
                value = cabs2(z1) * (z2 / cabs2(z2));
            }

            return value;
        }


        public static int triangle_upper_to_i4(int i, int j)

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UPPER_TO_I4 converts an upper triangular coordinate to an integer.
//
//  Discussion:
//
//    Triangular coordinates are handy when storing a naturally triangular
//    array (such as the upper half of a matrix) in a linear array.
//
//    Thus, for example, we might consider storing
//
//    (0,0) (0,1) (0,2) (0,3)
//          (1,1) (1,2) (1,3)
//                (2,2) (2,3)
//                      (3,3)
//
//    as the linear array
//
//    (0,0) (0,1) (1,1) (0,2) (1,2) (2,2) (0,3) (1,3) (2,3) (3,3)
//
//    Here, the quantities in parenthesis represent the natural row and
//    column indices of a single number when stored in a rectangular array.
//
//    Thus, our goal is, given the row I and column J of the data,
//    to produce the value K which indicates its position in the linear
//    array.
//
//    The triangular numbers are the indices associated with the
//    diagonal elements of the original array, T(0,0), T(1,1), T(2,2), 
//    T(3,3) and so on.
//
//    The formula is:
//
//      K = I + ( ( J * ( J + 1 ) ) / 2
//
//  Example:
//
//    I  J  K
//
//    0  0  0
//    0  1  1
//    1  1  2
//    0  2  3
//    1  2  4
//    2  2  5
//    0  3  6
//    1  3  7
//    2  3  8
//    3  3  9
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2017
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the row and column indices.  I and J must
//    be nonnegative, and I must not be greater than J.
//
//    Output, int TRIANGLE_UPPER_TO_I4, the linear index of the (I,J) element.
//
        {
            int value;

            if (i < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_UPPER_TO_I4 - Fatal error!");
                Console.WriteLine("  I < 0.");
                Console.WriteLine("  I = " + i + "");
                return (1);
            }
            else if (j < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_UPPER_TO_I4 - Fatal error!");
                Console.WriteLine("  J < 0.");
                Console.WriteLine("  J = " + j + "");
                return (1);
            }
            else if (j < i)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_UPPER_TO_I4 - Fatal error!");
                Console.WriteLine("  J < I.");
                Console.WriteLine("  I = " + i + "");
                Console.WriteLine("  J = " + j + "");
                return (1);
            }

            value = i + (j * (j + 1)) / 2;

            return value;
        }

        public static void xerbla(string srname, int info)

//****************************************************************************80
//
//  Purpose:
//
//    XERBLA is an error handler for the LAPACK routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, string SRNAME, the name of the routine
//    which called XERBLA.
//
//    Input, int INFO, the position of the invalid parameter in
//    the parameter list of the calling routine.
//
        {
            Console.WriteLine("");
            Console.WriteLine("XERBLA - Fatal error!");
            Console.WriteLine("  On entry to routine '" + srname + "',");
            Console.WriteLine("  input parameter number " + info + " had an illegal value.");
        }

        public static double zabs1(Complex z)

//****************************************************************************80
//
//  Purpose:
//
//    ZABS1 returns the L1 norm of a complex <double>.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    The L1 norm of a complex number is the sum of the absolute values
//    of the real and imaginary components.
//
//    ZABS1 ( Z ) = abs ( real ( Z ) ) + abs ( imaginary ( Z ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, complex <double> Z, the number whose norm is desired.
//
//    Output, double ZABS1, the L1 norm of Z.
//
        {
            double value = Math.Abs(z.Real) + Math.Abs(z.Imaginary);

            return value;
        }

        public static double zabs2(Complex z)

//****************************************************************************80
//
//  Purpose:
//
//    ZABS2 returns the L2 norm of a complex <double>.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    The L2 norm of a complex number is the square root of the sum
//    of the squares of the real and imaginary components.
//
//    ZABS2 ( Z ) = sqrt ( ( real ( Z ) )**2 + ( imaginary ( Z ) )**2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, complex <double> Z, the number whose norm is desired.
//
//    Output, float ZABS2, the L2 norm of Z.
//
        {
            double zr = z.Real;
            double zi = z.Imaginary;
            double value = Math.Sqrt(zr * zr + zi * zi);

            return value;
        }

        public static double zmach(int job)

//****************************************************************************80
//
//  Purpose:
//
//    ZMACH computes machine parameters for complex <double> arithmetic.
//
//  Discussion:
//
//    Assume the computer has
//
//      B = base of arithmetic;
//      T = number of base B digits;
//      L = smallest possible exponent;
//      U = largest possible exponent;
//
//    then
//
//      EPS = B^(1-T)
//      TINY = 100.0 * B^(-L+T)
//      HUGE = 0.01 * B^(U-T)
//
//    If complex division is done by
//
//      1 / (X+i*Y) = (X-i*Y) / (X^2+Y^2)
//
//    then
//
//      TINY = sqrt ( TINY )
//      HUGE = sqrt ( HUGE )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, int JOB:
//    1, EPS is desired;
//    2, TINY is desired;
//    3, HUGE is desired.
//
//    Output, double ZMACH, the requested value.
//
        {
            double s;
            double tiny;
            double value;

            double eps = 1.0;

            for (;;)
            {
                eps = eps / 2.0;
                s = 1.0 + eps;
                if (s <= 1.0)
                {
                    break;
                }
            }

            eps = 2.0 * eps;

            s = 1.0;

            for (;;)
            {
                tiny = s;
                s = s / 16.0;

                if (s * 1.0 == 0.0)
                {
                    break;
                }
            }

            tiny = (tiny / eps) * 100.0;
//
//  Had to insert this manually!
//
            tiny = Math.Sqrt(tiny);

            if (false)
            {
                Complex temp1 = new Complex(1.0, 0.0);
                Complex temp2 = new Complex(tiny, 0.0);
                Complex temp3 = temp1 / temp2;

                s = temp3.Real;

                if (s != 1.0 / tiny)
                {
                    tiny = Math.Sqrt(tiny);
                }
            }

            double huge = 1.0 / tiny;

            if (job == 1)
            {
                value = eps;
            }
            else if (job == 2)
            {
                value = tiny;
            }
            else if (job == 3)
            {
                value = huge;
            }
            else
            {
                value = 0.0;
            }

            return value;
        }

        public static Complex zsign1(Complex z1, Complex z2)

//****************************************************************************80
//
//  Purpose:
//
//    ZSIGN1 is a complex <double> transfer-of-sign function.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    The L1 norm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <double> Z1, Z2, the arguments.
//
//    Output, complex <double> ZSIGN1,  a complex value, with the magnitude of
//    Z1, and the argument of Z2.
//
        {
            Complex value;

            if (zabs1(z2) == 0.0)
            {
                value = new Complex(0.0, 0.0);
            }
            else
            {
                value = zabs1(z1) * (z2 / zabs1(z2));
            }

            return value;
        }

        public static Complex zsign2(Complex z1, Complex z2)

//****************************************************************************80
//
//  Purpose:
//
//    ZSIGN2 is a complex <double> transfer-of-sign function.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    The L2 norm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <double> Z1, Z2, the arguments.
//
//    Output, complex <double> ZSIGN2,  a complex value, with the magnitude of
//    Z1, and the argument of Z2.
//
        {
            Complex value;

            if (zabs2(z2) == 0.0)
            {
                value = new Complex(0.0, 0.0);
            }
            else
            {
                value = zabs2(z1) * (z2 / zabs2(z2));
            }

            return value;
        }


    }
}
