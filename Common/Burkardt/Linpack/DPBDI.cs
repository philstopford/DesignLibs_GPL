namespace Burkardt.Linpack;

public static class DPBDI
{
    public static void dpbdi ( double[] abd, int lda, int n, int m, ref double[] det )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPBDI computes the determinant of a matrix factored by DPBCO or DPBFA.
        //
        //  Discussion:
        //
        //    If the inverse is needed, use DPBSL N times.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2005
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
        //    Pete Stewart.
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
        //    Input, double ABD[LDA*N], the output from DPBCO or DPBFA.
        //
        //    Input, int LDA, the leading dimension of the array ABD.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int M, the number of diagonals above the main diagonal.
        //
        //    Output, double DET[2], the determinant of the original
        //    matrix in the form
        //      determinant = DET[0] * 10.0**DET[1]
        //    with 1.0D+00 <= DET[0] < 10.0D+00 or DET[0] == 0.0D+00.
        //
    {
        int i;
        //
        //  Compute the determinant.
        //
        det[0] = 1.0;
        det[1] = 0.0;
        const double s = 10.0;

        for ( i = 1; i <= n; i++ )
        {
            det[0] = det[0] * abd[m+(i-1)*lda] * abd[m+(i-1)*lda];

            switch (det[0])
            {
                case 0.0:
                    return;
            }

            while ( det[0] < 1.0 )
            {
                det[0] *= s;
                det[1] -= 1.0;
            }

            while ( s <= det[0] )
            {
                det[0] /= s;
                det[1] += 1.0;
            }

        }
    }
}