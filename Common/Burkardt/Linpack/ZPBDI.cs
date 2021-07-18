using System.Numerics;

namespace Burkardt.Linpack
{
    public static class ZPBDI
    {
        public static void zpbdi ( Complex[] abd, int lda, int n, int m, double[] det )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPBDI gets the determinant of a hermitian positive definite band matrix.
        //
        //  Discussion:
        //
        //    ZPBDI uses the factors computed by ZPBCO or ZPBFA.
        //
        //    If the inverse is needed, use ZPBSL N times.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2006
        //
        //  Author:
        //
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //
        //  Parameters:
        //
        //    Input, Complex ABD[LDA*N], the output from ZPBCO or ZPBFA.
        //
        //    Input, int LDA, the leading dimension of the array ABD.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int M, the number of diagonals above the main diagonal.
        //
        //    Output, double DET[2], the determinant of the original matrix in the
        //    form determinant = DET(1) * 10.0**DET(2) with 1.0 <= DET(1) < 10.0
        //    or DET(1) == 0.0.
        //
        {
            int i;

            det[0] = 1.0;
            det[1] = 0.0;

            for ( i = 1; i <= n; i++ )
            {
                det[0] = det[0] * ( abd[m+(i-1)*lda].Real ) * ( abd[m+(i-1)*lda].Real );

                if ( det[0] == 0.0 )
                {
                    break;
                }

                while ( det[0] < 1.0 )
                {
                    det[0] = det[0] * 10.0;
                    det[1] = det[1] - 1.0;
                }

                while ( 10.0 <= det[0] )
                {
                    det[0] = det[0] / 10.0;
                    det[1] = det[1] + 1.0;
                }

            }
        }

    }
}