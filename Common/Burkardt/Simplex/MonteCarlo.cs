using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SimplexNS
{
    public static class MonteCarlo
    {
        public static double[] simplex_general_sample(int m, int n, double[] t, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_GENERAL_SAMPLE samples a general simplex in M dimensions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 March 2017
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
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of points.
            //
            //    Input, double T[M*(M+1)], the simplex vertices.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double SIMPLEX_GENERAL_SAMPLE[M*N], the points.
            //
        {
            double[] x;
            double[] x1;

            x1 = simplex_unit_sample(m, n, ref seed);

            x = new double[m * n];

            simplex_unit_to_general(m, n, t, x1, x);

            return x;
        }

        public static double simplex_general_volume(int m, double[] t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_GENERAL_VOLUME computes the volume of a simplex in N dimensions.
            //
            //  Discussion:
            //
            //    The formula is: 
            //
            //      volume = 1/M! * det ( B )
            //
            //    where B is the M by M matrix obtained by subtracting one
            //    vector from all the others.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 March 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the dimension of the space.
            //
            //    Input, double T[M*(M+1)], the vertices.
            //
            //    Output, double SIMPLEX_GENERAL_VOLUME, the volume of the simplex.
            //
        {
            double[] b;
            double det;
            int i;
            int j;
            int[] pivot;
            double volume;

            pivot = new int[m];
            b = new double[m * m];

            for (j = 0; j < m; j++)
            {
                for (i = 0; i < m; i++)
                {
                    b[i + j * m] = t[i + j * m] - t[i + m * m];
                }
            }

            typeMethods.r8ge_fa(m, ref b, ref pivot);

            det = typeMethods.r8ge_det(m, b, pivot);

            volume = Math.Abs(det);
            for (i = 1; i <= m; i++)
            {
                volume = volume / (double)(i);
            }

            return volume;
        }

        public static double simplex_unit_monomial_integral(int m, int[] e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_MONOMIAL_INTEGRAL: integrals in the unit simplex in M dimensions.
            //
            //  Discussion:
            //
            //    The monomial is F(X) = product ( 1 <= I <= M ) X(I)^E(I).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int E[M], the exponents.  
            //    Each exponent must be nonnegative.
            //
            //    Output, double SIMPLEX_UNIT_MONOMIAL_INTEGRAL, the integral.
            //
        {
            int i;
            double integral;
            int j;
            int k;

            for (i = 0; i < m; i++)
            {
                if (e[i] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SIMPLEX_UNIT_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    return (1);
                }
            }

            k = 0;
            integral = 1.0;

            for (i = 0; i < m; i++)
            {
                for (j = 1; j <= e[i]; j++)
                {
                    k = k + 1;
                    integral = integral * (double)(j) / (double)(k);
                }
            }

            for (i = 0; i < m; i++)
            {
                k = k + 1;
                integral = integral / (double)(k);
            }

            return integral;
        }

        public static double[] simplex_unit_sample(int m, int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_SAMPLE samples the unit simplex in M dimensions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 January 2015
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
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, int &SEED, a seed for the random 
            //    number generator.
            //
            //    Output, double SIMPLEX_UNIT_SAMPLE_01[M*N], the points.
            //
        {
            double[] e;
            double e_sum;
            int i;
            int j;
            double[] x;

            x = new double[m * n];

            for (j = 0; j < n; j++)
            {
                e = UniformRNG.r8vec_uniform_01_new(m + 1, ref seed);

                for (i = 0; i < m + 1; i++)
                {
                    e[i] = -Math.Log(e[i]);
                }

                e_sum = typeMethods.r8vec_sum(m + 1, e);

                for (i = 0; i < m; i++)
                {
                    x[i + j * m] = e[i] / e_sum;
                }
            }

            return x;
        }

        public static void simplex_unit_to_general(int m, int n, double[] t, double[] ref_,
                double[] phy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_TO_GENERAL maps the unit simplex to a general simplex.
            //
            //  Discussion:
            //
            //    Given that the unit simplex has been mapped to a general simplex
            //    with vertices T, compute the images in T, under the same linear
            //    mapping, of points whose coordinates in the unit simplex are REF.
            //
            //    The vertices of the unit simplex are listed as suggested in the
            //    following:
            //
            //      (0,0,0,...,0)
            //      (1,0,0,...,0)
            //      (0,1,0,...,0)
            //      (0,0,1,...,0)
            //      (...........)
            //      (0,0,0,...,1)
            //
            //    Thanks to Andrei ("spiritualworlds") for pointing out a mistake in the
            //    previous implementation of this routine, 02 March 2008.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of points to transform.
            //
            //    Input, double T[M*(M+1)], the vertices of the
            //    general simplex.
            //
            //    Input, double REF[M*N], points in the
            //    reference triangle.
            //
            //    Output, double PHY[M*N], corresponding points in the physical triangle.
            //
        {
            int dim;
            int point;
            int vertex;
            //
            //  The image of each point is initially the image of the origin.
            //
            //  Insofar as the pre-image differs from the origin in a given vertex
            //  direction, add that proportion of the difference between the images
            //  of the origin and the vertex.
            //
            for (point = 0; point < n; point++)
            {
                for (dim = 0; dim < m; dim++)
                {
                    phy[dim + point * m] = t[dim + 0 * m];

                    for (vertex = 1; vertex < m + 1; vertex++)
                    {
                        phy[dim + point * m] = phy[dim + point * m]
                                               + (t[dim + vertex * m] - t[dim + 0 * m]) * ref_[vertex - 1 + point * m];
                    }
                }
            }

            return;
        }

        public static double simplex_unit_volume(int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_VOLUME returns the volume of the unit simplex in M dimensions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Output, double SIMPLEX_UNIT_VOLUME, the volume.
            //
        {
            int i;
            double volume;

            volume = 1.0;
            for (i = 1; i <= m; i++)
            {
                volume = volume / (double)(i);
            }

            return volume;
        }
    }
}