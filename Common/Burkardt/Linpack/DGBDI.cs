using System;

namespace Burkardt.Linpack;

public static class DGBDI
{
    public static void dgbdi(double[] abd, int lda, int n, int ml, int mu, int[] ipvt,
            ref double[] det )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGBDI computes the determinant of a band matrix factored by DGBCO or DGBFA.
        //
        //  Discussion:
        //
        //    If the inverse is needed, use DGBSL N times.
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
        //    Input, double ABD[LDA*N], the LU factor information from DGBCO 
        //    or DGBFA.
        //
        //    Input, int LDA, the leading dimension of the array ABD.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, MU, the number of diagonals below and above the
        //    main diagonal.  0 <= ML < N, 0 <= MU < N.
        //
        //    Input, int IPVT[N], the pivot vector from DGBCO or DGBFA.
        //
        //    Output, double DET[2], the determinant of the original matrix,
        //    if requested.
        //      determinant = DET[0] * 10.0**DET[1]
        //    with  1.0D+00 <= abs ( DET[0] ) < 10.0D+00 or DET[0] = 0.0D+00.
        //
    {
        int i;

        int m = ml + mu + 1;
        det[0] = 1.0;
        det[1] = 0.0;

        for (i = 1; i <= n; i++)
        {
            if (ipvt[i - 1] != i)
            {
                det[0] = -det[0];
            }

            det[0] *= abd[m - 1 + (i - 1) * lda];

            switch (det[0])
            {
                case 0.0:
                    return;
            }

            while (Math.Abs(det[0]) < 1.0)
            {
                det[0] *= 10.0;
                det[1] -= 1.0;
            }

            while (10.0 <= Math.Abs(det[0]))
            {
                det[0] /= 10.0;
                det[1] += 1.0;
            }
        }
    }
}