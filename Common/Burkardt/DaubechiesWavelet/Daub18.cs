using Burkardt.Types;

namespace Burkardt.DaubechiesWavelet;

public static class Daub18
{
    public static double[] daub18_transform(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB18_TRANSFORM computes the DAUB18 transform of a vector.
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
        //    Output, double DAUB18_TRANSFORM[N], the transformed vector.
        //
    {
        double[] c =  {
                3.807794736387834E-02,
                2.438346746125903E-01,
                6.048231236901111E-01,
                6.572880780513005E-01,
                1.331973858250075E-01,
                -2.932737832791749E-01,
                -9.684078322297646E-02,
                1.485407493381063E-01,
                3.072568147933337E-02,
                -6.763282906132997E-02,
                2.509471148314519E-04,
                2.236166212367909E-02,
                -4.723204757751397E-03,
                -4.281503682463429E-03,
                1.847646883056226E-03,
                2.303857635231959E-04,
                -2.519631889427101E-04,
                3.934732031627159E-05
            }
            ;
        const int p = 17;

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

    public static double[] daub18_transform_inverse(int n, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB18_TRANSFORM_INVERSE inverts the DAUB18 transform of a vector.
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
        //    Output, double DAUB18_TRANSFORM_INVERSE[N], the original vector.
        //
    {
        double[] c =  {
                3.807794736387834E-02,
                2.438346746125903E-01,
                6.048231236901111E-01,
                6.572880780513005E-01,
                1.331973858250075E-01,
                -2.932737832791749E-01,
                -9.684078322297646E-02,
                1.485407493381063E-01,
                3.072568147933337E-02,
                -6.763282906132997E-02,
                2.509471148314519E-04,
                2.236166212367909E-02,
                -4.723204757751397E-03,
                -4.281503682463429E-03,
                1.847646883056226E-03,
                2.303857635231959E-04,
                -2.519631889427101E-04,
                3.934732031627159E-05
            }
            ;
        const int p = 17;

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