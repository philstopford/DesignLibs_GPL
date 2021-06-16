using System;
using Burkardt.Types;

namespace Burkardt.Uniform
{
    public static class Ellipsoid
    {
        public static double[] uniform_in_ellipsoid_map(int dim_num, int n, double[] a, double r,
                ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    UNIFORM_IN_ELLIPSOID_MAP maps uniform points into an ellipsoid.
            //
            //  Discussion:
            //
            //    The points X in the ellipsoid are described by a DIM_NUM by DIM_NUM positive
            //    definite symmetric matrix A, and a "radius" R, such that
            //
            //      X' * A * X <= R * R
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
            //    23 May 2005
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
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Input, int N, the number of points.
            //
            //    Input, double A[DIM_NUM*DIM_NUM], the matrix that describes the ellipsoid.
            //
            //    Input, double R, the right hand side of the ellipsoid equation.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double UNIFORM_IN_ELLIPSOID_MAP[DIM_NUM*N], the points.
            //
        {
            int i;
            int info;
            int j;
            double[] u;
            double[] x;
            //
            //  Get the upper triangular Cholesky factor U of A.
            //
            u = new double[dim_num * dim_num];

            for (j = 0; j < dim_num; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    u[i + j * dim_num] = a[i + j * dim_num];
                }
            }

            info = typeMethods.r8po_fa(ref u, dim_num, dim_num);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("UNIFORM_IN_ELLIPSOID_MAP - Fatal error!");
                Console.WriteLine("  R8PO_FA reports that the matrix A");
                Console.WriteLine("  is not positive definite symmetric.");
                return null;
            }

            //
            //  Get the points Y that satisfy Y' * Y <= R * R.
            //
            x = Sphere.uniform_in_sphere01_map(dim_num, n, ref seed);

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    x[i + j * dim_num] = r * x[i + j * dim_num];
                }
            }

            //
            //  Solve U * X = Y.
            //
            for (j = 0; j < n; j++)
            {
                typeMethods.r8po_sl(u, dim_num, dim_num, ref x, bIndex: +j * dim_num);
            }

            return x;
        }

        public static double[] uniform_on_ellipsoid_map(int dim_num, int n, double[] a,
        double r, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_ON_ELLIPSOID_MAP maps uniform points onto an ellipsoid.
        //
        //  Discussion:
        //
        //    The points X on the ellipsoid are described by a DIM_NUM by DIM_NUM
        //    positive definite symmetric matrix A, and a "radius" R, such that
        //
        //      X' * A * X = R * R
        //
        //    The algorithm computes the Cholesky factorization of A:
        //
        //      A = U' * U.
        //
        //    A set of uniformly random points Y is generated, satisfying:
        //
        //      Y' * Y = R * R.
        //
        //    The appropriate points in the ellipsoid are found by solving
        //
        //      U * X = Y
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
        //    23 May 2005
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
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double A[DIM_NUM*DIM_NUM], the matrix that describes the ellipsoid.
        //
        //    Input, double R, the right hand side of the ellipsoid equation.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_ON_ELLIPSOID_MAP[DIM_NUM*N], the points.
        //
        {
            int i;
            int info;
            int j;
            double[] u;
            double[] x;
            //
            //  Get the factor U.
            //
            u = new double[dim_num * dim_num];

            for (j = 0; j < dim_num; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    u[i + j * dim_num] = a[i + j * dim_num];
                }
            }

            info = typeMethods.r8po_fa(ref a, dim_num, dim_num);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("UNIFORM_ON_ELLIPSOID_MAP - Fatal error!");
                Console.WriteLine("  R8PO_FA reports that the matrix A ");
                Console.WriteLine("  is not positive definite symmetric.");
                return null;
            }

            //
            //  Get the points Y that satisfy Y' * Y = R * R.
            //
            x = Sphere.uniform_on_sphere01_map(dim_num, n, ref seed);

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    x[i + j * dim_num] = r * x[i + j * dim_num];
                }
            }

            //
            //  Solve U * X = Y.
            //
            for (j = 0; j < n; j++)
            {
                typeMethods.r8po_sl(u, dim_num, dim_num, ref x, aIndex:0, bIndex: + j * dim_num);
            }

            return x;
        }
    }
}