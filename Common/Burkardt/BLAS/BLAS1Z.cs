using System;
using System.Numerics;

namespace Burkardt.BLAS
{
    public static class BLAS1Z
    {
        public static double dzasum(int n, Complex[] x, int incx, int index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DZASUM takes the sum of the absolute values of a complex <double> vector.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
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
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, complex <double> X[], the vector.
            //
            //    Input, int INCX, the increment between successive entries of X.
            //
            //    Output, double DZASUM, the sum of the absolute values.
            //
        {
            int i;
            int ix;
            double value;

            value = 0.0;

            if (n <= 0 || incx <= 0)
            {
                return value;
            }

            if (incx == 1)
            {
                for (i = 0; i < n; i++)
                {
                    value = value + Math.Abs(x[index + (i)].Real)
                                  + Math.Abs(x[index + (i)].Imaginary);
                }
            }
            else
            {
                ix = 0;
                for (i = 0; i < n; i++)
                {
                    value = value + Math.Abs(x[index = (ix)].Real)
                                  + Math.Abs(x[index + (ix)].Imaginary);
                    ix = ix + incx;
                }
            }

            return value;
        }

        public static double dznrm2(int n, Complex[] x, int incx, int index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DZNRM2 returns the euclidean norm of a complex <double> vector.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
            //
            //    DZNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
            //            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
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
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, complex <double> X[], the vector.
            //
            //    Input, int INCX, the increment between successive entries of X.
            //
            //    Output, double DZNRM2, the norm of the vector.
            //
        {
            int i;
            int ix;
            double scale;
            double ssq;
            double temp;
            double value;

            if (n < 1 || incx < 1)
            {
                value = 0.0;
            }
            else
            {
                scale = 0.0;
                ssq = 1.0;
                ix = 0;

                for (i = 0; i < n; i++)
                {
                    if (x[index + ix].Real != 0.0)
                    {
                        temp = Math.Abs(x[index + ix].Real);
                        if (scale < temp)
                        {
                            ssq = 1.0 + ssq * Math.Pow(scale / temp, 2);
                            scale = temp;
                        }
                        else
                        {
                            ssq = ssq + Math.Pow(temp / scale, 2);
                        }
                    }

                    if (x[index + ix].Imaginary != 0.0)
                    {
                        temp = Math.Abs(x[index + ix].Imaginary);
                        if (scale < temp)
                        {
                            ssq = 1.0 + ssq * Math.Pow(scale / temp, 2);
                            scale = temp;
                        }
                        else
                        {
                            ssq = ssq + Math.Pow(temp / scale, 2);
                        }
                    }

                    ix = ix + incx;
                }

                value = scale * Math.Sqrt(ssq);
            }

            return value;
        }

        public static int izamax(int n, Complex[] x, int incx, int index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    IZAMAX indexes the complex <double> vector element of maximum absolute value.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
            //
            //    WARNING: This index is a 1-based index, not a 0-based index!
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
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, complex <double> X[], the vector.
            //
            //    Input, int INCX, the increment between successive entries of X.
            //
            //    Output, int IZAMAX, the index of the element of maximum
            //    absolute value.
            //
        {
            int i;
            int ix;
            double smax;
            int value;

            value = 0;

            if (n < 1 || incx <= 0)
            {
                return value;
            }

            value = 1;

            if (n == 1)
            {
                return value;
            }

            if (incx != 1)
            {
                ix = 0;
                smax = BLAS0.zabs1(x[index + 0]);
                ix = ix + incx;

                for (i = 1; i < n; i++)
                {
                    if (smax < BLAS0.zabs1(x[index + ix]))
                    {
                        value = i + 1;
                        smax = BLAS0.zabs1(x[index + ix]);
                    }

                    ix = ix + incx;
                }
            }
            else
            {
                smax = BLAS0.zabs1(x[index + 0]);
                for (i = 1; i < n; i++)
                {
                    if (smax < BLAS0.zabs1(x[index + i]))
                    {
                        value = i + 1;
                        smax = BLAS0.zabs1(x[index + i]);
                    }
                }
            }

            return value;
        }

        public static void zaxpy(int n, Complex ca, Complex[] cx,
                int incx, ref Complex[] cy, int incy, int xIndex = 0, int yIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZAXPY adds a multiple of one complex <double> vector to another.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
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
            //    Basic Linear Algebra Subprograms for Fortran Usage,
            //    Algorithm 539,
            //    ACM Transactions on Mathematical Software,
            //    Volume 5, Number 3, September 1979, pages 308-323.
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements in CX and CY.
            //
            //    Input, complex <double> CA, the multiplier of CX.
            //
            //    Input, complex <double> CX[], the first vector.
            //
            //    Input, int INCX, the increment between successive entries of CX.
            //
            //    Input/output, complex <double> CY[], the second vector.
            //    On output, CY(*) has been replaced by CY(*) + CA * CX(*).
            //
            //    Input, int INCY, the increment between successive entries of CY.
            //
        {
            int i;
            int ix;
            int iy;

            if (n <= 0)
            {
                return;
            }

            if (BLAS0.zabs1(ca) == 0.0)
            {
                return;
            }

            if (incx != 1 || incy != 1)
            {
                if (0 <= incx)
                {
                    ix = 0;
                }
                else
                {
                    ix = (-n + 1) * incx;
                }

                if (0 <= incy)
                {
                    iy = 0;
                }
                else
                {
                    iy = (-n + 1) * incy;
                }

                for (i = 0; i < n; i++)
                {
                    cy[yIndex + iy] = cy[yIndex + iy] + ca * cx[xIndex + ix];
                    ix = ix + incx;
                    iy = iy + incy;
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    cy[yIndex + i] = cy[yIndex + i] + ca * cx[xIndex + i];
                }

            }
        }

        public static void zcopy(int n, Complex[] cx, int incx, ref Complex[] cy,
                int incy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZCOPY copies a complex <double> vector.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
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
            //    Basic Linear Algebra Subprograms for Fortran Usage,
            //    Algorithm 539,
            //    ACM Transactions on Mathematical Software,
            //    Volume 5, Number 3, September 1979, pages 308-323.
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements in CX and CY.
            //
            //    Input, complex <double> CX[], the first vector.
            //
            //    Input, int INCX, the increment between successive entries of CX.
            //
            //    Output, complex <double> CY[], the second vector.
            //
            //    Input, int INCY, the increment between successive entries of CY.
            //
        {
            int i;
            int ix;
            int iy;

            if (n <= 0)
            {
                return;
            }

            if (incx != 1 || incy != 1)
            {
                if (0 <= incx)
                {
                    ix = 0;
                }
                else
                {
                    ix = (-n + 1) * incx;
                }

                if (0 <= incy)
                {
                    iy = 0;
                }
                else
                {
                    iy = (-n + 1) * incy;
                }

                for (i = 0; i < n; i++)
                {
                    cy[iy] = cx[ix];
                    ix = ix + incx;
                    iy = iy + incy;
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    cy[i] = cx[i];
                }
            }

            return;
        }

        public static Complex zdotc(int n, Complex[] cx, int incx,
                Complex[] cy, int incy, int xIndex = 0, int yIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZDOTC forms the conjugated dot product of two complex <double> vectors.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
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
            //    Basic Linear Algebra Subprograms for Fortran Usage,
            //    Algorithm 539,
            //    ACM Transactions on Mathematical Software,
            //    Volume 5, Number 3, September 1979, pages 308-323.
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, complex <double> CX[], the first vector.
            //
            //    Input, int INCX, the increment between successive entries in CX.
            //
            //    Input, complex <double> CY[], the second vector.
            //
            //    Input, int INCY, the increment between successive entries in CY.
            //
            //    Output, complex <double> ZDOTC, the conjugated dot product of
            //    the corresponding entries of CX and CY.
            //
        {
            int i;
            int ix;
            int iy;
            Complex value;

            value = new Complex(0.0, 0.0);

            if (n <= 0)
            {
                return value;
            }

            if (incx == 1 && incy == 1)
            {
                for (i = 0; i < n; i++)
                {
                    value = value + Complex.Conjugate(cx[xIndex + i]) * cy[yIndex + i];
                }
            }
            else
            {
                if (0 <= incx)
                {
                    ix = 0;
                }
                else
                {
                    ix = (-n + 1) * incx;
                }

                if (0 <= incy)
                {
                    iy = 0;
                }
                else
                {
                    iy = (-n + 1) * incy;
                }

                for (i = 0; i < n; i++)
                {
                    value = value + Complex.Conjugate(cx[xIndex + ix]) * cy[yIndex + iy];
                    ix = ix + incx;
                    iy = iy + incy;
                }
            }

            return value;
        }

        public static Complex zdotu(int n, Complex[] cx, int incx,
                Complex[] cy, int incy, int xIndex = 0, int yIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZDOTU forms the unconjugated dot product of two complex <double> vectors.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
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
            //    Basic Linear Algebra Subprograms for Fortran Usage,
            //    Algorithm 539,
            //    ACM Transactions on Mathematical Software,
            //    Volume 5, Number 3, September 1979, pages 308-323.
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, complex <double> CX[], the first vector.
            //
            //    Input, int INCX, the increment between successive entries in CX.
            //
            //    Input, complex <double> CY[], the second vector.
            //
            //    Input, int INCY, the increment between successive entries in CY.
            //
            //    Output, complex <double> ZDOTU, the unconjugated dot product of
            //    the corresponding entries of CX and CY.
            //
        {
            int i;
            int ix;
            int iy;
            Complex value;

            value = new Complex(0.0, 0.0);

            if (n <= 0)
            {
                return value;
            }

            if (incx == 1 && incy == 1)
            {
                for (i = 0; i < n; i++)
                {
                    value = value + cx[xIndex + i] * cy[yIndex + i];
                }
            }
            else
            {
                if (0 <= incx)
                {
                    ix = 0;
                }
                else
                {
                    ix = (-n + 1) * incx;
                }

                if (0 <= incy)
                {
                    iy = 0;
                }
                else
                {
                    iy = (-n + 1) * incy;
                }

                for (i = 0; i < n; i++)
                {
                    value = value + cx[xIndex + ix] * cy[yIndex + iy];
                    ix = ix + incx;
                    iy = iy + incy;
                }
            }

            return value;
        }

        public static void zdrot(int n, ref Complex[] cx, int incx, ref Complex[] cy,
                int incy, double c, double s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZDROT applies a plane rotation to complex <double> vectors.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
            //
            //    The cosine C and sine S are real and the vectors CX and CY are complex.
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
            //    Basic Linear Algebra Subprograms for Fortran Usage,
            //    Algorithm 539,
            //    ACM Transactions on Mathematical Software,
            //    Volume 5, Number 3, September 1979, pages 308-323.
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input/output, complex <double> CX[], one of the vectors to be rotated.
            //
            //    Input, int INCX, the increment between successive entries of CX.
            //
            //    Input/output, complex <double> CY[], one of the vectors to be rotated.
            //
            //    Input, int INCY, the increment between successive elements of CY.
            //
            //    Input, double C, S, parameters (presumably the cosine and sine of
            //    some angle) that define a plane rotation.
            //
        {
            Complex ctemp;
            int i;
            int ix;
            int iy;

            if (n <= 0)
            {
                return;
            }

            if (incx == 1 && incy == 1)
            {
                for (i = 0; i < n; i++)
                {
                    ctemp = c * cx[i] + s * cy[i];
                    cy[i] = c * cy[i] - s * cx[i];
                    cx[i] = ctemp;
                }
            }
            else
            {
                if (0 <= incx)
                {
                    ix = 0;
                }
                else
                {
                    ix = (-n + 1) * incx;
                }

                if (0 <= incy)
                {
                    iy = 0;
                }
                else
                {
                    iy = (-n + 1) * incy;
                }

                for (i = 0; i < n; i++)
                {
                    ctemp = c * cx[ix] + s * cy[iy];
                    cy[iy] = c * cy[iy] - s * cx[ix];
                    cx[ix] = ctemp;
                    ix = ix + incx;
                    iy = iy + incy;
                }
            }
        }

        public static void zdscal(int n, double sa, ref Complex[] cx, int incx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZDSCAL scales a complex <double> vector by a double.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
            //
            //    The scaling constant is double precision real.
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
            //    Basic Linear Algebra Subprograms for Fortran Usage,
            //    Algorithm 539,
            //    ACM Transactions on Mathematical Software,
            //    Volume 5, Number 3, September 1979, pages 308-323.
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double SA, the multiplier.
            //
            //    Input/output, complex <double> CX[], the vector to be scaled.
            //
            //    Input, int INCX, the increment between successive entries of
            //    the vector CX.
            //
        {
            int i;

            if (n <= 0 || incx <= 0)
            {
                return;
            }

            if (incx == 1)
            {
                for (i = 0; i < n; i++)
                {
                    cx[i] = sa * cx[i];
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    cx[i * incx] = sa * cx[i * incx];
                }
            }

            return;
        }

        public static void zrotg(ref Complex ca, Complex cb, ref double c,
                ref Complex s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZROTG determines a complex <double> Givens rotation.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
            //
            //    Given values A and B, this routine computes:
            //
            //    If A = 0:
            //
            //      R = B
            //      C = 0
            //      S = (1,0).
            //
            //    If A /= 0:
            //
            //      ALPHA = A / abs ( A )
            //      NORM  = sqrt ( ( abs ( A ) )^2 + ( abs ( B ) )^2 )
            //      R     = ALPHA * NORM
            //      C     = abs ( A ) / NORM
            //      S     = ALPHA * conj ( B ) / NORM
            //
            //    In either case, the computed numbers satisfy the equation:
            //
            //    (         C    S ) * ( A ) = ( R )
            //    ( -conj ( S )  C )   ( B ) = ( 0 )
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
            //    Input/output, complex <double> *CA, on input, the value A.  On output,
            //    the value R.
            //
            //    Input, complex <double> CB, the value B.
            //
            //    Output, double *C, the cosine of the Givens rotation.
            //
            //    Output, complex <double> *S, the sine of the Givens rotation.
            //
        {
            Complex alpha;
            double norm;
            double scale;

            if (BLAS0.zabs2(ca) == 0.0)
            {
                c = 0.0;
                s = new Complex(1.0, 0.0);
                ca = cb;
            }
            else
            {
                scale = BLAS0.zabs2(ca) + BLAS0.zabs2(cb);
                norm = scale * Math.Sqrt(Math.Pow(BLAS0.zabs2(ca / scale), 2)
                                         + Math.Pow(BLAS0.zabs2(cb / scale), 2));
                alpha = ca / BLAS0.zabs2(ca);
                c = BLAS0.zabs2(ca) / norm;
                s = alpha * Complex.Conjugate(cb) / norm;
                ca = alpha * norm;
            }
        }

        public static void zscal(int n, Complex ca, ref Complex[] cx, int incx, int index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZSCAL scales a complex <double> vector by a complex <double>.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
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
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
            //    LINPACK User's Guide,
            //    SIAM, 1979.
            //
            //    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
            //    Basic Linear Algebra Subprograms for Fortran Usage,
            //    Algorithm 539,
            //    ACM Transactions on Mathematical Software,
            //    Volume 5, Number 3, September 1979, pages 308-323.
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, complex <double> CA, the multiplier.
            //
            //    Input/output, complex <double> CX[], the vector to be scaled.
            //
            //    Input, int INCX, the increment between successive entries of CX.
            //
        {
            int i;

            if (n <= 0 || incx <= 0)
            {
                return;
            }

            if (incx == 1)
            {
                for (i = 0; i < n; i++)
                {
                    cx[index + i] = ca * cx[index + i];
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    cx[index + (i * incx)] = ca * cx[index + (i * incx)];
                }
            }
        }

        public static void zswap(int n, ref Complex[] cx, int incx, ref Complex[] cy,
                int incy, int xIndex = 0, int yIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZSWAP interchanges two complex <double> vectors.
            //
            //  Discussion:
            //
            //    This routine uses double precision complex arithmetic.
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
            //    Basic Linear Algebra Subprograms for Fortran Usage,
            //    Algorithm 539,
            //    ACM Transactions on Mathematical Software,
            //    Volume 5, Number 3, September 1979, pages 308-323.
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input/output, complex <double> CX[], one of the vectors to swap.
            //
            //    Input, int INCX, the increment between successive entries of CX.
            //
            //    Input/output, complex <double> CY[], one of the vectors to swap.
            //
            //    Input, int INCY, the increment between successive elements of CY.
            //
        {
            Complex ctemp;
            int i;
            int ix;
            int iy;

            if (n <= 0)
            {
                return;
            }

            if (incx == 1 && incy == 1)
            {
                for (i = 0; i < n; i++)
                {
                    ctemp = cx[xIndex + i];
                    cx[xIndex + i] = cy[yIndex + i];
                    cy[yIndex + i] = ctemp;
                }
            }
            else
            {
                if (0 <= incx)
                {
                    ix = 0;
                }
                else
                {
                    ix = (-n + 1) * incx;
                }

                if (0 <= incy)
                {
                    iy = 0;
                }
                else
                {
                    iy = (-n + 1) * incy;
                }

                for (i = 0; i < n; i++)
                {
                    ctemp = cx[xIndex + ix];
                    cx[xIndex + ix] = cy[yIndex + iy];
                    cy[yIndex + iy] = ctemp;
                    ix = ix + incx;
                    iy = iy + incy;
                }
            }
        }
    }
}