using Burkardt.Types;

namespace Burkardt.DaubechiesWavelet;

public static class Daub20
{
    public static double[] daub20_transform(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB20_TRANSFORM computes the DAUB20 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vector.
        //    N must be a power of 2 and at least 4.
        //
        //    Input, double X[N], the vector to be transformed. 
        //
        //    Output, double DAUB20_TRANSFORM[N], the transformed vector.
        //
    {
        double[] c =  {
                2.667005790055555E-02,
                1.881768000776914E-01,
                5.272011889317255E-01,
                6.884590394536035E-01,
                2.811723436605774E-01,
                -2.498464243273153E-01,
                -1.959462743773770E-01,
                1.273693403357932E-01,
                9.305736460357235E-02,
                -7.139414716639708E-02,
                -2.945753682187581E-02,
                3.321267405934100E-02,
                3.606553566956169E-03,
                -1.073317548333057E-02,
                1.395351747052901E-03,
                1.992405295185056E-03,
                -6.858566949597116E-04,
                -1.164668551292854E-04,
                9.358867032006959E-05,
                -1.326420289452124E-05
            }
            ;
        const int p = 19;

        double[] y = typeMethods.r8vec_copy_new(n, x);
        double[] z = new double[n];

        int m = n;

        while (4 <= m)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                z[i] = 0.0;
            }

            i = 0;

            int j;
            for (j = 0; j < m - 1; j += 2)
            {
                int k;
                for (k = 0; k < p; k += 2)
                {
                    int j0 = typeMethods.i4_wrap(j + k, 0, m - 1);
                    int j1 = typeMethods.i4_wrap(j + k + 1, 0, m - 1);
                    z[i] = z[i] + c[k] * y[j0] + c[k + 1] * y[j1];
                    z[i + m / 2] = z[i + m / 2] + c[p - k] * y[j0] - c[p - k - 1] * y[j1];
                }

                i += 1;
            }

            for (i = 0; i < m; i++)
            {
                y[i] = z[i];
            }

            m /= 2;
        }

        return y;
    }

    public static double[] daub20_transform_inverse(int n, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB20_TRANSFORM_INVERSE inverts the DAUB20 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vector.
        //    N must be a power of 2 and at least 4.
        //
        //    Input, double Y[N], the transformed vector. 
        //
        //    Output, double DAUB20_TRANSFORM_INVERSE[N], the original vector.
        //
    {
        double[] c =  {
                2.667005790055555E-02,
                1.881768000776914E-01,
                5.272011889317255E-01,
                6.884590394536035E-01,
                2.811723436605774E-01,
                -2.498464243273153E-01,
                -1.959462743773770E-01,
                1.273693403357932E-01,
                9.305736460357235E-02,
                -7.139414716639708E-02,
                -2.945753682187581E-02,
                3.321267405934100E-02,
                3.606553566956169E-03,
                -1.073317548333057E-02,
                1.395351747052901E-03,
                1.992405295185056E-03,
                -6.858566949597116E-04,
                -1.164668551292854E-04,
                9.358867032006959E-05,
                -1.326420289452124E-05
            }
            ;
        const int p = 19;

        double[] x = typeMethods.r8vec_copy_new(n, y);
        double[] z = new double[n];

        int m = 4;
        int q = (p - 1) / 2;

        while (m <= n)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                z[i] = 0.0;
            }

            int j = 0;

            for (i = -q; i < m / 2 - q; i++)
            {
                int k;
                for (k = 0; k < p; k += 2)
                {
                    int i0 = typeMethods.i4_wrap(i + k / 2, 0, m / 2 - 1);
                    int i1 = typeMethods.i4_wrap(i + m / 2 + k / 2, m / 2, m - 1);
                    z[j] = z[j] + c[p - k - 1] * x[i0] + c[k + 1] * x[i1];
                    z[j + 1] = z[j + 1] + c[p - k] * x[i0] - c[k] * x[i1];
                }

                j += 2;
            }

            for (i = 0; i < m; i++)
            {
                x[i] = z[i];
            }

            m *= 2;
        }

        return x;
    }
}