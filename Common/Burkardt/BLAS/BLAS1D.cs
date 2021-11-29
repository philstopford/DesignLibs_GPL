using System;
using Burkardt.Types;

namespace Burkardt.BLAS;

public static class BLAS1D
{
    public static double dasum(int n, double[] x, int incx, int index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DASUM takes the sum of the absolute values of a vector.
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
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
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
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double X[*], the vector to be examined.
        //
        //    Input, int INCX, the increment between successive entries of X.
        //    INCX must not be negative.
        //
        //    Output, double DASUM, the sum of the absolute values of X.
        //
    {
        int i;

        double value = 0.0;
        int j = 0;

        for (i = 0; i < n; i++)
        {
            value += Math.Abs(x[index + j]);
            j += incx;
        }

        return value;
    }

    public static void daxpy(int n, double da, double[] dx, int incx, ref double[] dy,
            int incy, int xIndex = 0, int yIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAXPY computes constant times a vector plus a vector.
        //
        //  Discussion:
        //
        //    This routine uses unrolled loops for increments equal to one.
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
        //    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
        //    David Kincaid, Fred Krogh.
        //    C++ version by John Burkardt.
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
        //    Basic Linear Algebra Subprograms for Fortran Usage,
        //    Algorithm 539,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, September 1979, pages 308-323.
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements in DX and DY.
        //
        //    Input, double DA, the multiplier of DX.
        //
        //    Input, double DX[*], the first vector.
        //
        //    Input, int INCX, the increment between successive entries of DX.
        //
        //    Input/output, double DY[*], the second vector.
        //    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
        //
        //    Input, int INCY, the increment between successive entries of DY.
        //
    {
        switch (n)
        {
            case <= 0:
                return;
        }

        switch (da)
        {
            case 0.0:
                return;
        }

        //
        //  Code for unequal increments or equal increments
        //  not equal to 1.
        //
        if (incx != 1 || incy != 1)
        {
            int ix = incx switch
            {
                >= 0 => 0,
                _ => (-n + 1) * incx
            };

            int iy = incy switch
            {
                >= 0 => 0,
                _ => (-n + 1) * incy
            };

            for (int i = 0; i < n; i++)
            {
                dy[iy + yIndex] += da * dx[ix + xIndex];
                ix += incx;
                iy += incy;
            }
        }
        //
        //  Code for both increments equal to 1.
        //
        else
        {
            int m = n % 4;

            for (int i = 0; i < m; i++)
            {
                dy[i + yIndex] += da * dx[i + xIndex];
            }

            for (int i = m; i < n; i += 4)
            {
                dy[i + yIndex] += da * dx[i + xIndex];
                dy[i + 1 + yIndex] += da * dx[i + 1 + xIndex];
                dy[i + 2 + yIndex] += da * dx[i + 2 + xIndex];
                dy[i + 3 + yIndex] += da * dx[i + 3 + xIndex];
            }

        }
    }

    public static void dcopy(int n, double[] dx, int incx, ref double[] dy, int incy, int xIndex = 0, int yIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DCOPY copies a vector X to a vector Y.
        //
        //  Discussion:
        //
        //    The routine uses unrolled loops for increments equal to one.
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
        //    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
        //    David Kincaid, Fred Krogh.
        //    C++ version by John Burkardt.
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
        //    Basic Linear Algebra Subprograms for Fortran Usage,
        //    Algorithm 539,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, September 1979, pages 308-323.
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements in DX and DY.
        //
        //    Input, double DX[*], the first vector.
        //
        //    Input, int INCX, the increment between successive entries of DX.
        //
        //    Output, double DY[*], the second vector.
        //
        //    Input, int INCY, the increment between successive entries of DY.
        //
    {
        int i;

        switch (n)
        {
            case <= 0:
                return;
        }

        switch (incx)
        {
            case 1 when incy == 1:
            {
                int 
                    m = n % 7;

                if (m != 0)
                {
                    for (i = 0; i < m; i++)
                    {
                        dy[yIndex + i] = dx[xIndex + i];
                    }
                }

                for (i = m; i < n; i += 7)
                {
                    dy[yIndex + i] = dx[xIndex + i];
                    dy[yIndex + i + 1] = dx[xIndex + i + 1];
                    dy[yIndex + i + 2] = dx[xIndex + i + 2];
                    dy[yIndex + i + 3] = dx[xIndex + i + 3];
                    dy[yIndex + i + 4] = dx[xIndex + i + 4];
                    dy[yIndex + i + 5] = dx[xIndex + i + 5];
                    dy[yIndex + i + 6] = dx[xIndex + i + 6];
                }

                break;
            }
            default:
            {
                int ix = incx switch
                {
                    >= 0 => 0,
                    _ => (-n + 1) * incx
                };

                int iy = incy switch
                {
                    >= 0 => 0,
                    _ => (-n + 1) * incy
                };

                for (i = 0; i < n; i++)
                {
                    dy[yIndex + iy] = dx[xIndex + ix];
                    ix += incx;
                    iy += incy;
                }

                break;
            }
        }
    }

    public static double ddot(int n, double[] dx, int incx, double[] dy,
            int incy, int xIndex = 0, int yIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DDOT forms the dot product of two vectors.
        //
        //  Discussion:
        //
        //    This routine uses unrolled loops for increments equal to one.
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
        //    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
        //    David Kincaid, Fred Krogh.
        //    C++ version by John Burkardt.
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
        //    Basic Linear Algebra Subprograms for Fortran Usage,
        //    Algorithm 539,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, September 1979, pages 308-323.
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input, double DX[*], the first vector.
        //
        //    Input, int INCX, the increment between successive entries in DX.
        //
        //    Input, double DY[*], the second vector.
        //
        //    Input, int INCY, the increment between successive entries in DY.
        //
        //    Output, double DDOT, the sum of the product of the corresponding
        //    entries of DX and DY.
        //
    {
        int i;

        double dtemp = 0.0;

        switch (n)
        {
            case <= 0:
                return dtemp;
        }

        //
        //  Code for unequal increments or equal increments
        //  not equal to 1.
        //
        if (incx != 1 || incy != 1)
        {
            int ix = incx switch
            {
                >= 0 => 0,
                _ => (-n + 1) * incx
            };

            int iy = incy switch
            {
                >= 0 => 0,
                _ => (-n + 1) * incy
            };

            for (i = 0; i < n; i++)
            {
                dtemp += dx[xIndex + ix] * dy[yIndex + iy];
                ix += incx;
                iy += incy;
            }
        }
        //
        //  Code for both increments equal to 1.
        //
        else
        {
            int m = n % 5;

            for (i = 0; i < m; i++)
            {
                dtemp += dx[xIndex + i] * dy[yIndex + i];
            }

            for (i = m; i < n; i += 5)
            {
                dtemp = dtemp + dx[xIndex + i] * dy[yIndex + i]
                              + dx[xIndex + i + 1] * dy[yIndex + i + 1]
                              + dx[xIndex + i + 2] * dy[yIndex + i + 2]
                              + dx[xIndex + i + 3] * dy[yIndex + i + 3]
                              + dx[xIndex + i + 4] * dy[yIndex + i + 4];
            }

        }

        return dtemp;
    }

    public static double dnrm2(int n, double[] x, int incx, int index = 0)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DNRM2 returns the euclidean norm of a vector.
        //
        //  Discussion:
        //
        //     DNRM2 ( X ) = sqrt ( X' * X )
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
        //    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
        //    David Kincaid, Fred Krogh.
        //    C++ version by John Burkardt.
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
        //    Basic Linear Algebra Subprograms for Fortran Usage,
        //    Algorithm 539,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, September 1979, pages 308-323.
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double X[*], the vector whose norm is to be computed.
        //
        //    Input, int INCX, the increment between successive entries of X.
        //
        //    Output, double DNRM2, the Euclidean norm of X.
        //
    {
        double norm;

        if (n < 1 || incx < 1)
        {
            norm = 0.0;
        }
        else
        {
            switch (n)
            {
                case 1:
                    norm = Math.Abs(x[index + 0]);
                    break;
                default:
                {
                    double scale = 0.0;
                    double ssq = 1.0;
                    int ix = 0;

                    int i;
                    for (i = 0; i < n; i++)
                    {
                        if (x[index + ix] != 0.0)
                        {
                            double absxi = Math.Abs(x[index + ix]);
                            if (scale < absxi)
                            {
                                ssq = 1.0 + ssq * (scale / absxi) * (scale / absxi);
                                scale = absxi;
                            }
                            else
                            {
                                ssq += absxi / scale * (absxi / scale);
                            }
                        }

                        ix += incx;
                    }

                    norm = scale * Math.Sqrt(ssq);
                    break;
                }
            }
        }

        return norm;
    }

    public static void drot(int n, ref double[] x, int incx, ref double[] y, int incy, double c,
            double s, int xindex = 0, int yindex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DROT applies a plane rotation.
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
        //    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
        //    David Kincaid, Fred Krogh.
        //    C++ version by John Burkardt.
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
        //    Basic Linear Algebra Subprograms for Fortran Usage,
        //    Algorithm 539,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, September 1979, pages 308-323.
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input/output, double X[*], one of the vectors to be rotated.
        //
        //    Input, int INCX, the increment between successive entries of X.
        //
        //    Input/output, double Y[*], one of the vectors to be rotated.
        //
        //    Input, int INCY, the increment between successive elements of Y.
        //
        //    Input, double C, S, parameters (presumably the cosine and
        //    sine of some angle) that define a plane rotation.
        //
    {
        switch (n)
        {
            case <= 0:
                break;
            default:
            {
                double stemp;
                int i;
                switch (incx)
                {
                    case 1 when incy == 1:
                    {
                        for (i = 0; i < n; i++)
                        {
                            stemp = c * x[i + xindex] + s * y[i + yindex];
                            y[i + yindex] = c * y[i + yindex] - s * x[i + xindex];
                            x[i + xindex] = stemp;
                        }

                        break;
                    }
                    default:
                    {
                        int ix = incx switch
                        {
                            >= 0 => 0,
                            _ => (-n + 1) * incx
                        };

                        int iy = incy switch
                        {
                            >= 0 => 0,
                            _ => (-n + 1) * incy
                        };

                        for (i = 0; i < n; i++)
                        {
                            stemp = c * x[ix + xindex] + s * y[iy + yindex];
                            y[iy + yindex] = c * y[iy + yindex] - s * x[ix + xindex];
                            x[ix + xindex] = stemp;
                            ix += incx;
                            iy += incy;
                        }

                        break;
                    }
                }

                break;
            }
        }
    }

    public static void drotg(ref double sa, ref double sb, ref double c, ref double s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DROTG constructs a Givens plane rotation.
        //
        //  Discussion:
        //
        //    Given values A and B, this routine computes
        //
        //    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
        //          = sign ( B ) if abs ( A ) <= abs ( B );
        //
        //    R     = SIGMA * ( A * A + B * B );
        //
        //    C = A / R if R is not 0
        //      = 1     if R is 0;
        //
        //    S = B / R if R is not 0,
        //        0     if R is 0.
        //
        //    The computed numbers then satisfy the equation
        //
        //    (  C  S ) ( A ) = ( R )
        //    ( -S  C ) ( B ) = ( 0 )
        //
        //    The routine also computes
        //
        //    Z = S     if abs ( A ) > abs ( B ),
        //      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
        //      = 1     if C is 0.
        //
        //    The single value Z encodes C and S, and hence the rotation:
        //
        //    If Z = 1, set C = 0 and S = 1;
        //    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
        //    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
        //    David Kincaid, Fred Krogh.
        //    C++ version by John Burkardt.
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
        //    Basic Linear Algebra Subprograms for Fortran Usage,
        //    Algorithm 539,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, September 1979, pages 308-323.
        //
        //  Parameters:
        //
        //    Input/output, double *SA, *SB,  On input, SA and SB are the values
        //    A and B.  On output, SA is overwritten with R, and SB is
        //    overwritten with Z.
        //
        //    Output, double *C, *S, the cosine and sine of the Givens rotation.
        //
    {
        double r;

        double roe = Math.Abs(sb) < Math.Abs(sa) ? sa : sb;

        double scale = Math.Abs(sa) + Math.Abs(sb);

        switch (scale)
        {
            case 0.0:
                c = 1.0;
                s = 0.0;
                r = 0.0;
                break;
            default:
                r = scale * Math.Sqrt(sa / scale * (sa / scale)
                                      + sb / scale * (sb / scale));
                r = typeMethods.r8_sign(roe) * r;
                c = sa / r;
                s = sb / r;
                break;
        }

        double z = Math.Abs(c) switch
        {
            > 0.0 when Math.Abs(c) <= s => 1.0 / c,
            _ => s
        };

        sa = r;
        sb = z;
    }

    public static void dscal(int n, double sa, ref double[] x, int incx, int index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSCAL scales a vector by a constant.
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
        //    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
        //    David Kincaid, Fred Krogh.
        //    C++ version by John Burkardt.
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
        //    Input/output, double X[*], the vector to be scaled.
        //
        //    Input, int INCX, the increment between successive entries of X.
        //
    {
        switch (n)
        {
            case <= 0:
                break;
            default:
            {
                switch (incx)
                {
                    case 1:
                    {
                        int m = n % 5;

                        for (int i = 0; i < m; i++)
                        {
                            x[i + index] = sa * x[i + index];
                        }

                        for (int i = m; i < n; i += 5)
                        {
                            x[i + index] = sa * x[i + index];
                            x[i + 1 + index] = sa * x[i + 1 + index];
                            x[i + 2 + index] = sa * x[i + 2 + index];
                            x[i + 3 + index] = sa * x[i + 3 + index];
                            x[i + 4 + index] = sa * x[i + 4 + index];
                        }

                        break;
                    }
                    default:
                    {
                        int ix = incx switch
                        {
                            >= 0 => 0,
                            _ => (-n + 1) * incx
                        };

                        for (int i = 0; i < n; i++)
                        {
                            x[ix + index] = sa * x[ix + index];
                            ix += incx;
                        }

                        break;
                    }
                }

                break;
            }
        }
    }

    public static int dsvdc(ref double[] a, int lda, int m, int n, ref double[] s, ref double[] e,
            ref double[] u, int ldu, ref double[] v, int ldv, double[] work, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSVDC computes the singular value decomposition of a real rectangular matrix.
        //
        //  Discussion:
        //
        //    This routine reduces an M by N matrix A to diagonal form by orthogonal
        //    transformations U and V.  The diagonal elements S(I) are the singular
        //    values of A.  The columns of U are the corresponding left singular
        //    vectors, and the columns of V the right singular vectors.
        //
        //    The form of the singular value decomposition is then
        //
        //      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 May 2007
        //
        //  Author:
        //
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input/output, double A[LDA*N].  On input, the M by N matrix whose
        //    singular value decomposition is to be computed.  On output, the matrix
        //    has been destroyed.  Depending on the user's requests, the matrix may 
        //    contain other useful information.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //    LDA must be at least M.
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix A.
        //
        //    Output, double S[MM], where MM = min(M+1,N).  The first
        //    min(M,N) entries of S contain the singular values of A arranged in
        //    descending order of magnitude.
        //
        //    Output, double E[MM], where MM = min(M+1,N), ordinarily contains zeros.
        //    However see the discussion of INFO for exceptions.
        //
        //    Output, double U[LDU*K].  If JOBA = 1 then K = M;
        //    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of left singular
        //    vectors.  U is not referenced if JOBA = 0.  If M <= N or if JOBA = 2, then
        //    U may be identified with A in the subroutine call.
        //
        //    Input, int LDU, the leading dimension of the array U.
        //    LDU must be at least M.
        //
        //    Output, double V[LDV*N], the N by N matrix of right singular vectors.
        //    V is not referenced if JOB is 0.  If N <= M, then V may be identified
        //    with A in the subroutine call.
        //
        //    Input, int LDV, the leading dimension of the array V.
        //    LDV must be at least N.
        //
        //    Workspace, double WORK[M].
        //
        //    Input, int JOB, controls the computation of the singular
        //    vectors.  It has the decimal expansion AB with the following meaning:
        //      A =  0, do not compute the left singular vectors.
        //      A =  1, return the M left singular vectors in U.
        //      A >= 2, return the first min(M,N) singular vectors in U.
        //      B =  0, do not compute the right singular vectors.
        //      B =  1, return the right singular vectors in V.
        //
        //    Output, int *DSVDC, status indicator INFO.
        //    The singular values (and their corresponding singular vectors)
        //    S(*INFO+1), S(*INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
        //    Thus if *INFO is 0, all the singular values and their vectors are
        //    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
        //    matrix with the elements of S on its diagonal and the elements of E on
        //    its superdiagonal.  Thus the singular values of A and B are the same.
        //
    {
        double cs = 0;
        int i;
        int j;
        int l;
        int ll;
        int ls = 0;
        const int maxit = 30;
        double sn = 0;
        double t;
        //
        //  Determine what is to be computed.
        //
        int info = 0;
        bool wantu = false;
        bool wantv = false;
        int jobu = job % 100 / 10;

        int ncu = jobu switch
        {
            > 1 => Math.Min(m, n),
            _ => m
        };

        if (jobu != 0)
        {
            wantu = true;
        }

        if (job % 10 != 0)
        {
            wantv = true;
        }

        //
        //  Reduce A to bidiagonal form, storing the diagonal elements
        //  in S and the super-diagonal elements in E.
        //
        int nct = Math.Min(m - 1, n);
        int nrt = Math.Max(0, Math.Min(m, n - 2));
        int lu = Math.Max(nct, nrt);

        for (l = 1; l <= lu; l++)
        {
            //
            //  Compute the transformation for the L-th column and
            //  place the L-th diagonal in S(L).
            //
            if (l <= nct)
            {
                s[l - 1] = dnrm2(m - l + 1, a, 1, + l - 1 + (l - 1) * lda);

                if (s[l - 1] != 0.0)
                {
                    if (a[l - 1 + (l - 1) * lda] != 0.0)
                    {
                        s[l - 1] = typeMethods.r8_sign(a[l - 1 + (l - 1) * lda]) * Math.Abs(s[l - 1]);
                    }

                    dscal(m - l + 1, 1.0 / s[l - 1], ref a, 1, + l - 1 + (l - 1) * lda);
                    a[l - 1 + (l - 1) * lda] = 1.0 + a[l - 1 + (l - 1) * lda];
                }

                s[l - 1] = -s[l - 1];
            }

            for (j = l + 1; j <= n; j++)
            {
                //
                //  Apply the transformation.
                //
                if (l <= nct && s[l - 1] != 0.0)
                {
                    t = -ddot(m - l + 1, a, 1, a, 1, + l - 1 + (l - 1) * lda, + l - 1 + (j - 1) * lda)
                        / a[l - 1 + (l - 1) * lda];
                    daxpy(m - l + 1, t, a, 1, ref a, 1, + l - 1 + (l - 1) * lda, + l - 1 + (j - 1) * lda);
                }

                //
                //  Place the L-th row of A into E for the
                //  subsequent calculation of the row transformation.
                //
                e[j - 1] = a[l - 1 + (j - 1) * lda];
            }

            switch (wantu)
            {
                //
                //  Place the transformation in U for subsequent back multiplication.
                //
                case true when l <= nct:
                {
                    for (i = l; i <= m; i++)
                    {
                        u[i - 1 + (l - 1) * ldu] = a[i - 1 + (l - 1) * lda];
                    }

                    break;
                }
            }

            if (l <= nrt)
            {
                //
                //  Compute the L-th row transformation and place the
                //  L-th superdiagonal in E(L).
                //
                e[l - 1] = dnrm2(n - l, e, 1, + l);

                if (e[l - 1] != 0.0)
                {
                    if (e[l] != 0.0)
                    {
                        e[l - 1] = typeMethods.r8_sign(e[l]) * Math.Abs(e[l - 1]);
                    }

                    dscal(n - l, 1.0 / e[l - 1], ref e, 1, + l);
                    e[l] = 1.0 + e[l];
                }

                e[l - 1] = -e[l - 1];
                //
                //  Apply the transformation.
                //
                if (l + 1 <= m && e[l - 1] != 0.0)
                {
                    for (j = l + 1; j <= m; j++)
                    {
                        work[j - 1] = 0.0;
                    }

                    for (j = l + 1; j <= n; j++)
                    {
                        daxpy(m - l, e[j - 1], a, 1, ref work, 1, + l + (j - 1) * lda, + l);
                    }

                    for (j = l + 1; j <= n; j++)
                    {
                        daxpy(m - l, -e[j - 1] / e[l], work, 1, ref a, 1, + l, + l + (j - 1) * lda);
                    }
                }

                switch (wantv)
                {
                    //
                    //  Place the transformation in V for subsequent back multiplication.
                    //
                    case true:
                    {
                        for (j = l + 1; j <= n; j++)
                        {
                            v[j - 1 + (l - 1) * ldv] = e[j - 1];
                        }

                        break;
                    }
                }
            }
        }

        //
        //  Set up the final bidiagonal matrix of order MN.
        //
        int mn = Math.Min(m + 1, n);
        int nctp1 = nct + 1;
        int nrtp1 = nrt + 1;

        if (nct < n)
        {
            s[nctp1 - 1] = a[nctp1 - 1 + (nctp1 - 1) * lda];
        }

        if (m < mn)
        {
            s[mn - 1] = 0.0;
        }

        if (nrtp1 < mn)
        {
            e[nrtp1 - 1] = a[nrtp1 - 1 + (mn - 1) * lda];
        }

        e[mn - 1] = 0.0;
        switch (wantu)
        {
            //
            //  If required, generate U.
            //
            case true:
            {
                for (i = 1; i <= m; i++)
                {
                    for (j = nctp1; j <= ncu; j++)
                    {
                        u[i - 1 + (j - 1) * ldu] = 0.0;
                    }
                }

                for (j = nctp1; j <= ncu; j++)
                {
                    u[j - 1 + (j - 1) * ldu] = 1.0;
                }

                for (ll = 1; ll <= nct; ll++)
                {
                    l = nct - ll + 1;

                    if (s[l - 1] != 0.0)
                    {
                        for (j = l + 1; j <= ncu; j++)
                        {
                            t = -ddot(m - l + 1, u, 1, u, 1, + (l - 1) + (l - 1) * ldu, + (l - 1) + (j - 1) * ldu)
                                / u[l - 1 + (l - 1) * ldu];
                            daxpy(m - l + 1, t, u, 1, ref u, 1, + (l - 1) + (l - 1) * ldu, + (l - 1) + (j - 1) * ldu);
                        }

                        dscal(m - l + 1, -1.0, ref u, 1, + (l - 1) + (l - 1) * ldu);
                        u[l - 1 + (l - 1) * ldu] = 1.0 + u[l - 1 + (l - 1) * ldu];
                        for (i = 1; i <= l - 1; i++)
                        {
                            u[i - 1 + (l - 1) * ldu] = 0.0;
                        }
                    }
                    else
                    {
                        for (i = 1; i <= m; i++)
                        {
                            u[i - 1 + (l - 1) * ldu] = 0.0;
                        }

                        u[l - 1 + (l - 1) * ldu] = 1.0;
                    }
                }

                break;
            }
        }

        switch (wantv)
        {
            //
            //  If it is required, generate V.
            //
            case true:
            {
                for (ll = 1; ll <= n; ll++)
                {
                    l = n - ll + 1;

                    if (l <= nrt && e[l - 1] != 0.0)
                    {
                        for (j = l + 1; j <= n; j++)
                        {
                            t = -ddot(n - l, v, 1, v, 1, + l + (l - 1) * ldv, + l + (j - 1) * ldv)
                                / v[l + (l - 1) * ldv];
                            daxpy(n - l, t, v, 1, ref v, 1, + l + (l - 1) * ldv, + l + (j - 1) * ldv);
                        }

                    }

                    for (i = 1; i <= n; i++)
                    {
                        v[i - 1 + (l - 1) * ldv] = 0.0;
                    }

                    v[l - 1 + (l - 1) * ldv] = 1.0;
                }

                break;
            }
        }

        //
        //  Main iteration loop for the singular values.
        //
        int mm = mn;
        int iter = 0;

        while (0 < mn)
        {
            //
            //  If too many iterations have been performed, set flag and return.
            //
            if (maxit <= iter)
            {
                info = mn;
                return info;
            }

            //
            //  This section of the program inspects for
            //  negligible elements in the S and E arrays.
            //
            //  On completion the variables KASE and L are set as follows:
            //
            //  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
            //  KASE = 2     if S(L) is negligible and L < MN
            //  KASE = 3     if E(L-1) is negligible, L < MN, and
            //               S(L), ..., S(MN) are not negligible (QR step).
            //  KASE = 4     if E(MN-1) is negligible (convergence).
            //
            double ztest = 0;
            double test = 0;
            for (ll = 1; ll <= mn; ll++)
            {
                l = mn - ll;

                if (l == 0)
                {
                    break;
                }

                test = Math.Abs(s[l - 1]) + Math.Abs(s[l]);
                ztest = test + Math.Abs(e[l - 1]);

                if (!(Math.Abs(ztest - test) <= double.Epsilon))
                {
                    continue;
                }

                e[l - 1] = 0.0;
                break;
            }

            int kase = 0;
            if (l == mn - 1)
            {
                kase = 4;
            }
            else
            {
                int lls = 0;
                for (lls = l + 1; lls <= mn + 1; lls++)
                {
                    ls = mn - lls + l + 1;

                    if (ls == l)
                    {
                        break;
                    }

                    test = 0.0;
                    if (ls != mn)
                    {
                        test += Math.Abs(e[ls - 1]);
                    }

                    if (ls != l + 1)
                    {
                        test += Math.Abs(e[ls - 2]);
                    }

                    ztest = test + Math.Abs(s[ls - 1]);

                    if (!(Math.Abs(ztest - test) <= double.Epsilon))
                    {
                        continue;
                    }

                    s[ls - 1] = 0.0;
                    break;

                }

                if (ls == l)
                {
                    kase = 3;
                }
                else if (ls == mn)
                {
                    kase = 1;
                }
                else
                {
                    kase = 2;
                    l = ls;
                }
            }

            l += 1;
            //
            //  Deflate negligible S(MN).
            //
            double f = 0;
            int mm1 = 0;
            int k = 0;
            double t1 = 0;
            switch (kase)
            {
                case 1:
                {
                    mm1 = mn - 1;
                    f = e[mn - 2];
                    e[mn - 2] = 0.0;

                    int kk = 0;
                    for (kk = 1; kk <= mm1; kk++)
                    {
                        k = mm1 - kk + l;
                        t1 = s[k - 1];
                        drotg(ref t1, ref f, ref cs, ref sn);
                        s[k - 1] = t1;

                        if (k != l)
                        {
                            f = -sn * e[k - 2];
                            e[k - 2] = cs * e[k - 2];
                        }

                        switch (wantv)
                        {
                            case true:
                                drot(n, ref v, 1, ref v, 1, cs, sn, + 0 + (k - 1) * ldv, + 0 + (mn - 1) * ldv);
                                break;
                        }
                    }

                    break;
                }
                //
                //  Split at negligible S(L).
                //
                case 2:
                {
                    f = e[l - 2];
                    e[l - 2] = 0.0;

                    for (k = l; k <= mn; k++)
                    {
                        t1 = s[k - 1];
                        drotg(ref t1, ref f, ref cs, ref sn);
                        s[k - 1] = t1;
                        f = -sn * e[k - 1];
                        e[k - 1] = cs * e[k - 1];
                        switch (wantu)
                        {
                            case true:
                                drot(m, ref u, 1, ref u, 1, cs, sn, + 0 + (k - 1) * ldu, + 0 + (l - 2) * ldu);
                                break;
                        }
                    }

                    break;
                }
                //
                //  Perform one QR step.
                //
                case 3:
                {
                    //
                    //  Calculate the shift.
                    //
                    double scale = Math.Max(Math.Abs(s[mn - 1]),
                        Math.Max(Math.Abs(s[mn - 2]),
                            Math.Max(Math.Abs(e[mn - 2]),
                                Math.Max(Math.Abs(s[l - 1]), Math.Abs(e[l - 1])))));

                    double sm = s[mn - 1] / scale;
                    double smm1 = s[mn - 2] / scale;
                    double emm1 = e[mn - 2] / scale;
                    double sl = s[l - 1] / scale;
                    double el = e[l - 1] / scale;
                    double b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.0;
                    double c = sm * emm1 * (sm * emm1);
                    double shift = 0.0;

                    if (b != 0.0 || c != 0.0)
                    {
                        shift = b switch
                        {
                            < 0.0 => -shift,
                            _ => Math.Sqrt(b * b + c)
                        };

                        shift = c / (b + shift);
                    }

                    f = (sl + sm) * (sl - sm) - shift;
                    double g = sl * el;
                    //
                    //  Chase zeros.
                    //
                    mm1 = mn - 1;

                    for (k = l; k <= mm1; k++)
                    {
                        drotg(ref f, ref g, ref cs, ref sn);

                        if (k != l)
                        {
                            e[k - 2] = f;
                        }

                        f = cs * s[k - 1] + sn * e[k - 1];
                        e[k - 1] = cs * e[k - 1] - sn * s[k - 1];
                        g = sn * s[k];
                        s[k] = cs * s[k];

                        switch (wantv)
                        {
                            case true:
                                drot(n, ref v, 1, ref v, 1, cs, sn, + 0 + (k - 1) * ldv, + 0 + k * ldv);
                                break;
                        }

                        drotg(ref f, ref g, ref cs, ref sn);
                        s[k - 1] = f;
                        f = cs * e[k - 1] + sn * s[k];
                        s[k] = -sn * e[k - 1] + cs * s[k];
                        g = sn * e[k];
                        e[k] = cs * e[k];

                        switch (wantu)
                        {
                            case true when k < m:
                                drot(m, ref u, 1, ref u, 1, cs, sn, + 0 + (k - 1) * ldu, + 0 + k * ldu);
                                break;
                        }
                    }

                    e[mn - 2] = f;
                    iter += 1;
                    break;
                }
                //
                //  Convergence.
                //
                case 4:
                {
                    switch (s[l - 1])
                    {
                        //
                        //  Make the singular value nonnegative.
                        //
                        case < 0.0:
                        {
                            s[l - 1] = -s[l - 1];
                            switch (wantv)
                            {
                                case true:
                                    dscal(n, -1.0, ref v, 1, + 0 + (l - 1) * ldv);
                                    break;
                            }

                            break;
                        }
                    }

                    //
                    //  Order the singular value.
                    //
                    for (;;)
                    {
                        if (l == mm)
                        {
                            break;
                        }

                        if (s[l] <= s[l - 1])
                        {
                            break;
                        }

                        t = s[l - 1];
                        s[l - 1] = s[l];
                        s[l] = t;

                        switch (wantv)
                        {
                            case true when l < n:
                                dswap(n, ref v, 1, ref v, 1, + 0 + (l - 1) * ldv, + 0 + l * ldv);
                                break;
                        }

                        switch (wantu)
                        {
                            case true when l < m:
                                dswap(m, ref u, 1, ref u, 1, + 0 + (l - 1) * ldu, + 0 + l * ldu);
                                break;
                        }

                        l += 1;
                    }

                    iter = 0;
                    mn -= 1;
                    break;
                }
            }
        }

        return info;
    }


    public static void dswap(int n, ref double[] x, int incx, ref double[] y, int incy, int xIndex = 0, int yIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSWAP interchanges two vectors.
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
        //    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
        //    David Kincaid, Fred Krogh.
        //    C++ version by John Burkardt.
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
        //    Basic Linear Algebra Subprograms for Fortran Usage,
        //    Algorithm 539,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, September 1979, pages 308-323.
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input/output, double X[*], one of the vectors to swap.
        //
        //    Input, int INCX, the increment between successive entries of X.
        //
        //    Input/output, double Y[*], one of the vectors to swap.
        //
        //    Input, int INCY, the increment between successive elements of Y.
        //
    {
        switch (n)
        {
            case <= 0:
                break;
            default:
            {
                int i;
                double temp;
                switch (incx)
                {
                    case 1 when incy == 1:
                    {
                        int m = n % 3;

                        for (i = 0; i < m; i++)
                        {
                            temp = x[i + xIndex];
                            x[i + xIndex] = y[i + yIndex];
                            y[i + yIndex] = temp;
                        }

                        for (i = m; i < n; i += 3)
                        {
                            temp = x[i + xIndex];
                            x[i + xIndex] = y[i + yIndex];
                            y[i + yIndex] = temp;

                            temp = x[i + 1 + xIndex];
                            x[i + 1 + xIndex] = y[i + 1 + yIndex];
                            y[i + 1 + yIndex] = temp;

                            temp = x[i + 2 + xIndex];
                            x[i + 2 + xIndex] = y[i + 2 + yIndex];
                            y[i + 2 + yIndex] = temp;
                        }

                        break;
                    }
                    default:
                    {
                        int ix = incx switch
                        {
                            >= 0 => 0,
                            _ => (-n + 1) * incx
                        };

                        int iy = incy switch
                        {
                            >= 0 => 0,
                            _ => (-n + 1) * incy
                        };

                        for (i = 0; i < n; i++)
                        {
                            temp = x[ix + xIndex];
                            x[ix + xIndex] = y[iy + yIndex];
                            y[iy + yIndex] = temp;
                            ix += incx;
                            iy += incy;
                        }

                        break;
                    }
                }

                break;
            }
        }
    }

    public static int idamax(int n, double[] dx, int incx, int index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    IDAMAX finds the index of the vector element of maximum absolute value.
        //
        //  Discussion:
        //
        //    WARNING: This index is a 1-based index, not a 0-based index!
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
        //    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
        //    David Kincaid, Fred Krogh.
        //    C++ version by John Burkardt.
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
        //    Basic Linear Algebra Subprograms for Fortran Usage,
        //    Algorithm 539,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, September 1979, pages 308-323.
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double X[*], the vector to be examined.
        //
        //    Input, int INCX, the increment between successive entries of SX.
        //
        //    Output, int IDAMAX, the index of the element of maximum
        //    absolute value.
        //
    {
        double dmax;
        int i;

        int value = 0;

        if (n < 1 || incx <= 0)
        {
            return value;
        }

        value = 1;

        switch (n)
        {
            case 1:
                return value;
        }

        switch (incx)
        {
            case 1:
            {
                dmax = Math.Abs(dx[0 + index]);

                for (i = 1; i < n; i++)
                {
                    if (!(dmax < Math.Abs(dx[i + index])))
                    {
                        continue;
                    }

                    value = i + 1;
                    dmax = Math.Abs(dx[i + index]);
                }

                break;
            }
            default:
            {
                int ix = 0;
                dmax = Math.Abs(dx[0 + index]);
                ix += incx;

                for (i = 1; i < n; i++)
                {
                    if (dmax < Math.Abs(dx[ix + index]))
                    {
                        value = i + 1;
                        dmax = Math.Abs(dx[ix + index]);
                    }

                    ix += incx;
                }

                break;
            }
        }

        return value;
    }

}