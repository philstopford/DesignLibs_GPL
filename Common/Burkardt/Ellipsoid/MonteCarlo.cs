using System;
using Burkardt.HyperGeometry.HypersphereNS;
using Burkardt.Types;

namespace Burkardt.Ellipsoid
{
    public static class MonteCarlo
    {
        public static double[] ellipsoid_sample(int m, int n, double[] a, double[] v, double r,
        ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSOID_SAMPLE samples points uniformly from an ellipsoid.
        //
        //  Discussion:
        //
        //    The points X in the ellipsoid are described by a M by M
        //    positive definite symmetric matrix A, a "center" V, and 
        //    a "radius" R, such that
        //
        //      (X-V)' * A * (X-V) <= R * R
        //
        //    The algorithm computes the Cholesky factorization of A:
        //
        //      A = U' * U.
        //
        //    A set of uniformly random points Y is generated, satisfying:
        //
        //      Y' * Y <= R * R.
        //
        //    The appropriate points in the ellipsoid are found by solving
        //
        //      U * X = Y
        //      X = X + V
        //
        //    Thanks to Dr Karl-Heinz Keil for pointing out that the original
        //    coding was actually correct only if A was replaced by its inverse.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation, and Sensitivity
        //    of Queueing Networks,
        //    Krieger, 1992,
        //    ISBN: 0894647644,
        //    LC: QA298.R79.
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double A[M*M], the matrix that describes
        //    the ellipsoid.
        //
        //    Input, double V[M], the "center" of the ellipsoid.
        //
        //    Input, double R, the "radius" of the ellipsoid.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double ELLIPSE_SAMPLE[M*N], the points.
        //
        {
            int i;
            int j;
            double[] t;
            double[] u;
            double[] x;
            //
            //  Get the Cholesky factor U.
            //
            try
            {
                u = typeMethods.r8po_fa(m, a);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("ELLIPSOID_SAMPLE - Fatal error!");
                Console.WriteLine("  R8PO_FA reports that the matrix A");
                Console.WriteLine("  is not positive definite symmetric.");
                return null;
            }

            //
            //  Get the points Y that satisfy Y' * Y <= 1.
            //
            x = Uniform.Sphere.uniform_in_sphere01_map(m, n, ref seed);
            //
            //  Get the points Y that satisfy Y' * Y <= R * R.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    x[i + j * m] = r * x[i + j * m];
                }
            }

            //
            //  Solve U * X = Y.
            //
            for (j = 0; j < n; j++)
            {
                t = typeMethods.r8po_sl(m, u, x, bIndex: + j * m);
                for (i = 0; i < m; i++)
                {
                    x[i + j * m] = t[i];
                }
            }

            //
            //  X = X + V.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    x[i + j * m] = x[i + j * m] + v[i];
                }
            }

            return x;
        }

        public static double ellipsoid_volume(int m, double[] a, double[] v, double r )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSOID_VOLUME returns the volume of an ellipsoid.
        //
        //  Discussion:
        //
        //    The points X in the ellipsoid are described by an M by M
        //    positive definite symmetric matrix A, an M-dimensional point V,
        //    and a "radius" R, such that
        //      (X-V)' * A * (X-V) <= R * R
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, double A[M*M], the matrix that describes
        //    the ellipsoid.  A must be symmetric and positive definite.
        //
        //    Input, double V[M], the "center" of the ellipse.
        //    The value of V is not actually needed by this function.
        //
        //    Input, double R, the "radius" of the ellipse.
        //
        //    Output, double ELLIPSOID_VOLUME, the volume of the ellipsoid.
        //
        {
            int i;
            double sqrt_det;
            double[] u;
            double volume;

            u = typeMethods.r8po_fa(m, a);

            sqrt_det = 1.0;
            for (i = 0; i < m; i++)
            {
                sqrt_det = sqrt_det * u[i + i * m];
            }

            volume = Math.Pow(r, m) * Hypersphere.hypersphere_unit_volume(m) / sqrt_det;
            
            return volume;
        }
    }
}