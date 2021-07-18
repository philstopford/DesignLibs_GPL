using System.Numerics;
using Burkardt.Types;

namespace Burkardt.Linpack
{
    public static class ZGBDI
    {
        public static void zgbdi(Complex[] abd, int lda, int n, int ml, int mu, int[] ipvt,
            ref Complex[] det )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZGBDI computes the determinant of a band matrix factored by ZGBCO or ZGBFA.
        //
        //  Discussion:
        //
        //    If the inverse is needed, use ZGBSL N times.
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
        //    Input, Complex ABD[LDA*N], the output from ZGBCO or ZGBFA.
        //
        //    Input, int LDA, the leading dimension of A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, the number of diagonals below the main diagonal.
        //
        //    Input, int MU, the number of diagonals above the main diagonal.
        //
        //    Input, int IPVT[N], the pivot vector from ZGBCO or ZGBFA.
        //
        //    Output, Complex DET[2], determinant of original matrix.
        //    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= zabs1 ( DET(1) ) < 10.0
        //    or DET(1) = 0.0.  Also, DET(2) is strictly real.
        //
        {
            int i;
            int m;

            m = ml + mu + 1;
            det[0] = new Complex(1.0, 0.0);
            det[1] = new Complex(0.0, 0.0);

            for (i = 1; i <= n; i++)
            {
                if (ipvt[i - 1] != i)
                {
                    det[0] = -det[0];
                }

                det[0] = det[0] * abd[m - 1 + (i - 1) * lda];

                if (typeMethods.zabs1(det[0]) == 0.0)
                {
                    break;
                }

                while (typeMethods.zabs1(det[0]) < 1.0)
                {
                    det[0] = det[0] * new Complex(10.0, 0.0);
                    det[1] = det[1] - new Complex(1.0, 0.0);
                }

                while (10.0 <= typeMethods.zabs1(det[0]))
                {
                    det[0] = det[0] / new Complex(10.0, 0.0);
                    det[1] = det[1] + new Complex(1.0, 0.0);
                }

            }
        }

    }
}