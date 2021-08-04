using System;
using Burkardt.Types;

namespace Burkardt.Ellipse
{
    public static class MonteCarlo
    {
        public static double ellipse_area1(double[] a, double r )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_AREA1 returns the area of an ellipse defined by a matrix.
        //
        //  Discussion:
        //
        //    The points X in the ellipse are described by a 2 by 2
        //    positive definite symmetric matrix A, and a "radius" R, such that
        //      X' * A * X <= R * R
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2*2], the matrix that describes
        //    the ellipse.  A must be symmetric and positive definite.
        //
        //    Input, double R, the "radius" of the ellipse.
        //
        //    Output, double ELLIPSE_AREA1, the area of the ellipse.
        //
        {
            const double r8_pi = 3.141592653589793;
            double value;

            value = r * r * r8_pi / Math.Sqrt(a[0 + 0 * 2] * a[1 + 1 * 2] - a[1 + 0 * 2] * a[0 + 1 * 2]);

            return value;
        }

        public static double ellipse_area2(double a, double b, double c, double d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPSE_AREA2 returns the area of an ellipse defined by an equation.
            //
            //  Discussion:
            //
            //    The ellipse is described by the formula
            //      a x^2 + b xy + c y^2 = d
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 November 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, B, C, coefficients on the left hand side.
            //
            //    Input, double D, the right hand side.
            //
            //    Output, double ELLIPSE_AREA2, the area of the ellipse.
            //
        {
            const double r8_pi = 3.141592653589793;
            double value;

            value = 2.0 * d * d * r8_pi / Math.Sqrt(4.0 * a * c - b * b);

            return value;
        }

        public static double[] ellipse_sample(int n, double[] a, double r, ref typeMethods.r8vecNormalData data, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_SAMPLE samples points in an ellipse.
        //
        //  Discussion:
        //
        //    The points X in the ellipsoid are described by a 2 by 2 positive
        //    definite symmetric matrix A, and a "radius" R, such that
        //      X' * A * X <= R * R
        //    The algorithm computes the Cholesky factorization of A:
        //      A = U' * U.
        //    A set of uniformly random points Y is generated, satisfying:
        //      Y' * Y <= R * R.
        //    The appropriate points in the ellipsoid are found by solving
        //      U * X = Y
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
        //    Input, int N, the number of points.
        //
        //    Input, double A[2*2], the matrix that describes the ellipse.
        //
        //    Input, double R, the right hand side of the ellipse equation.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double ELLIPSE_SAMPLE[2*N], the points.
        //
        {
            int i;
            int info;
            int j;
            int m = 2;
            double[] u;
            double[] x;
            //
            //  Get the upper triangular Cholesky factor U of A.
            //
            u = new double[m * m];

            for (j = 0; j < m; j++)
            {
                for (i = 0; i < m; i++)
                {
                    u[i + j * m] = a[i + j * m];
                }
            }

            info = typeMethods.r8po_fa(ref u, m, m);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("ELLIPSE_SAMPLE - Fatal error!");
                Console.WriteLine("  R8PO_FA reports that the matrix A");
                Console.WriteLine("  is not positive definite symmetric.");
                return null;
            }

            //
            //  Get the points Y that satisfy Y' * Y <= R * R.
            //
            x = Uniform.Sphere.uniform_in_sphere01_map(m, n, ref data, ref seed);

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
                typeMethods.r8po_sl(u, m, m, ref x, bIndex: + j * m);
            }
            
            return x;
        }
    }
}