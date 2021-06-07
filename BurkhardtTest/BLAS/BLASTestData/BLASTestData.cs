using System;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.BLASTestData
{
    public static class BLASData
    {
        public static double dmach(int job)

//****************************************************************************80
//
//  Purpose:
//
//    DMACH computes machine parameters of double precision real arithmetic.
//
//  Discussion:
//
//    This routine is for testing only.  It is not required by LINPACK.
//
//    If there is trouble with the automatic computation of these quantities,
//    they can be set by direct assignment statements.
//
//    We assume the computer has
//
//      B = base of arithmetic;
//      T = number of base B digits;
//      L = smallest possible exponent;
//      U = largest possible exponent.
//
//    then
//
//      EPS = B^(1-T)
//      TINY = 100.0 * B^(-L+T)
//      HUGE = 0.01 * B^(U-T)
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
//    Input, int JOB:
//    1: requests EPS;
//    2: requests TINY;
//    3: requests HUGE.
//
//    Output, double DMACH, the requested value.
//
        {
            double eps;
            double huge;
            double s;
            double tiny;
            double value;

            eps = 1.0;
            for (;;)
            {
                value = 1.0 + (eps / 2.0);
                if (value <= 1.0)
                {
                    break;
                }

                eps = eps / 2.0;
            }

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
            huge = 1.0 / tiny;

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
                Console.WriteLine("");
                Console.WriteLine("DMACH - Fatal error!");
                Console.WriteLine("  Illegal input value of JOB = " + job + "");
                return 0;
            }

            return value;
        }

        public static float smach(int job)

//****************************************************************************80
//
//  Purpose:
//
//    SMACH computes machine parameters of float arithmetic.
//
//  Discussion:
//
//    This routine is for testing only.  It is not required by LINPACK.
//
//    If there is trouble with the automatic computation of these quantities,
//    they can be set by direct assignment statements.
//
//    We assume the computer has
//
//      B = base of arithmetic;
//      T = number of base B digits;
//      L = smallest possible exponent;
//      U = largest possible exponent.
//
//    then
//
//      EPS = B**(1-T)
//      TINY = 100.0 * B**(-L+T)
//      HUGE = 0.01 * B**(U-T)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
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
//    Input, int JOB:
//    1: requests EPS;
//    2: requests TINY;
//    3: requests HUGE.
//
//    Output, float SMACH, the requested value.
//
        {
            float eps;
            float huge;
            float s;
            float tiny;
            float value;

            eps = 1.0f;
            for (;;)
            {
                value = 1.0f + (eps / 2.0f);
                if (value <= 1.0)
                {
                    break;
                }

                eps = eps / 2.0f;
            }

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
                typeMethods.xerbla("SMACH", 1);
            }

            return value;
        }
        
        
        public static double[] r8mat_test ( char trans, int lda, int m, int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_TEST sets up a test matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //  
            //  Modified:
            //
            //    10 February 2014
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, char TRANS, indicates whether matrix is to be transposed.
            //    'N', no transpose.
            //    'T', transpose the matrix.
            //
            //    Input, int LDA, the leading dimension of the matrix.
            //
            //    Input, int M, N, the number of rows and columns of the matrix.
            //
            //    Output, double R8MAT_TEST[?], the matrix.
            //    if TRANS is 'N', then the matrix is stored in LDA*N entries,
            //    as an M x N matrix;
            //    if TRANS is 'T', then the matrix is stored in LDA*M entries,
            //    as an N x M matrix.
            //
        {
            double[] a;

            if ( trans == 'N' )
            {
                a = new double[lda*n];

                for (int j = 0; j < n; j++ )
                {
                    for (int i = 0; i < m; i++ )
                    {
                        a[i+j*lda] = ( double ) ( 10 * ( i + 1 ) + ( j + 1 ) );
                    }
                }
            }
            else
            {
                a = new double[lda*m];

                for (int j = 0; j < n; j++ )
                {
                    for (int i = 0; i < m; i++ )
                    {
                        a[j+i*lda] = ( double ) ( 10 * ( i + 1 ) + ( j + 1 ) );
                    }
                }
            }
            return a;
        }
        
    }
}