using System;
using Burkardt.Types;

namespace Burkardt.Stroud
{
    public class Parallelipiped
    {
        public static double parallelipiped_volume_3d(double[] x, double[] y, double[] z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PARALLELIPIPED_VOLUME_3D returns the volume of a parallelipiped in 3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X[4], Y[4], Z[4], the coordinates of one corner
            //    of the parallelipiped, and its 3 immediate neighbors.
            //
            //    Output, double PARALLELIPIPED_VOLUME_3D, the volume of
            //    the parallelipiped.
            //
        {
            double volume;

            volume = Math.Abs(
                (z[1] - z[0]) * (y[3] * x[2] - y[2] * x[3]) +
                (z[2] - z[0]) * (x[3] * y[1] - x[1] * y[3]) +
                (z[3] - z[0]) * (x[1] * y[2] - x[2] * y[1]) +
                (z[2] - z[1]) * (y[3] * x[0] - y[0] * x[3]) +
                (z[3] - z[1]) * (x[2] * y[0] - x[0] * y[2]) +
                (z[3] - z[2]) * (x[0] * y[1] - x[1] * y[0]));

            return volume;
        }

        public static double parallelipiped_volume_nd(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PARALLELIPIPED_VOLUME_ND returns the volume of a parallelipiped in ND.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double V[N*(N+1)], the N+1 columns of V contains the 
            //    coordinates of the "corners" of the parallelipiped.
            //
            //    Output, double PARALLELIPIPED_VOLUME_ND, the volume of
            //    the parallelipiped.
            //
        {
            double det;
            int i;
            int info;
            int j;
            int[] pivot;
            double volume;
            double[] w;
            //
            //  Compute the volume of the N-dimensional parallelipiped.
            //
            w = new double[n * n];

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    w[i + j * n] = v[i + (j + 1) * n] - v[i + 0 * n];
                }
            }

            pivot = new int[n];

            info = typeMethods.r8ge_fa(n, ref w, ref pivot);

            if (info != 0)
            {
                volume = 0.0;
            }
            else
            {
                det = typeMethods.r8ge_det(n, w, pivot);

                volume = Math.Abs(det);
            }

            return volume;
        }


    }
}