using System;
using Burkardt.MonomialNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SimplexNS
{
    public static class Simplex
    {

        public static int simplex_num(int m, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_NUM evaluates the N-th Simplex number in M dimensions.
            //
            //  Discussion:
            //
            //     N\M: 1    2    3    4    5
            //    --   --   --   --   --   --
            //     0    0    0    0    0    0
            //     1    1    1    1    1    1
            //     2    2    3    4    5    6
            //     3    3    6   10   15   21
            //     4    4   10   20   35   56
            //     5    5   15   35   70  126
            //     6    6   21   56  126  252
            //     7    7   28   84  210  462
            //     8    8   36  120  330  792
            //     9    9   45  165  495 1287
            //    10   10   55  220  715 2002
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the index of the number.
            //
            //    Output, int SIMPLEX_NUM, the desired value.
            //
        {
            int i;
            int value;

            value = 1;
            for (i = 1; i <= m; i++)
            {
                value = (value * (n + i - 1)) / i;
            }

            return value;
        }

        public static double[] simplex_to_triangle(double[] tvert1, double[] tvert2,
                double[] tvert3, double[] s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_TO_TRIANGLE maps points from the simplex to a triangle.
            //
            //  Discussion:
            //
            //    The simplex has vertices:
            //
            //      (  0, 0 )
            //      (  1, 0 )
            //      (  0, 1 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 June 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double TVERT1[2], TVERT2[2], TVERT3[2], the coordinates
            //    of the vertices of the triangle.  These vertices will be taken
            //    to be the images of (0,0), (1,0) and (0,1) respectively.
            //
            //    Input, double S[2], the coordinates of the point in the simplex.
            //
            //    Output, double SIMPLEX_TO_TRIANGLE[2], the coordinates of the point in
            //    the triangle.
            //
        {
            int i;
            double[] t;

            t = new double[2];

            for (i = 0; i < 2; i++)
            {
                t[i] = tvert1[i] * (1.0 - s[0] - s[1])
                       + tvert2[i] * s[0]
                       + tvert3[i] * s[1];
            }

            return t;
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

        public static double simplex_unit_monomial_integral(int m, int[] expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_MONOMIAL_INTEGRAL integrates a monomial over a simplex.
            //
            //  Discussion:
            //
            //    This routine evaluates a monomial of the form
            //
            //      product ( 1 <= dim <= m ) x(dim)^expon(dim)
            //
            //    where the exponents are nonnegative integers.  Note that
            //    if the combination 0^0 is encountered, it should be treated
            //    as 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int EXPON[M], the exponents.
            //
            //    Output, double SIMPLEX_UNIT_MONOMIAL_INTEGRAL, the value of the integral
            //    of the monomial.
            //
        {
            int dim;
            int i;
            int k;
            double value;

            value = 1.0;

            k = 0;

            for (dim = 0; dim < m; dim++)
            {
                for (i = 1; i <= expon[dim]; i++)
                {
                    k = k + 1;
                    value = value * (double)(i) / (double)(k);
                }
            }

            for (dim = 0; dim < m; dim++)
            {
                k = k + 1;
                value = value / (double)(k);
            }

            return value;
        }

        public static double simplex_unit_monomial_quadrature(int m, int[] expon, int n, double[] x,
                double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_MONOMIAL_QUADRATURE: quadrature of monomials in a unit simplex.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int EXPON[M], the exponents.
            //
            //    Input, int N, the number of points in the rule.
            //
            //    Input, double X[M*N], the quadrature points.
            //
            //    Input, double W[N], the quadrature weights.
            //
            //    Output, double SIMPLEX_UNIT_MONOMIAL_QUADRATURE, the quadrature error.
            //
        {
            double exact = 1.0;
            double quad;
            double quad_error;
            double scale;
            double[] value;
            //
            //  Get the exact value of the integral of the unscaled monomial.
            //
            scale = simplex_unit_monomial_integral(m, expon);
            //
            //  Evaluate the monomial at the quadrature points.
            //
            value = Monomial.monomial_value(m, n, expon, x);
            //
            //  Compute the weighted sum and divide by the exact value.
            //
            quad = typeMethods.r8vec_dot_product(n, w, value) / scale;

            //
            //  Error:
            //
            quad_error = Math.Abs(quad - exact);

            return quad_error;
        }

        public static double[] simplex_unit_sample(int m, int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_SAMPLE returns uniformly random points from a general simplex.
            //
            //  Discussion:
            //
            //    The interior of the unit M dimensional simplex is the set of 
            //    points X(1:M) such that each X(I) is nonnegative, and 
            //    sum(X(1:M)) <= 1.
            //
            //    This routine is valid for any spatial dimension M.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2004
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
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double UNIFORM_IN_SIMPLEX01_MAP[M*N], the points.
            //
        {
            double[] e;
            int i;
            int j;
            double total;
            double[] x;
            //
            //  The construction begins by sampling M+1 points from the
            //  exponential distribution with parameter 1.
            //
            x = new double[m * n];

            for (j = 0; j < n; j++)
            {
                e = UniformRNG.r8vec_uniform_01_new(m + 1, ref seed);

                for (i = 0; i <= m; i++)
                {
                    e[i] = -Math.Log(e[i]);
                }

                total = 0.0;
                for (i = 0; i <= m; i++)
                {
                    total = total + e[i];
                }

                for (i = 0; i < m; i++)
                {
                    x[i + m * j] = e[i] / total;
                }
            }

            return x;
        }

        public static void simplex_unit_to_general(int m, int n, double[] t, double[] ref_,
                ref double[] phy)

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
            int i;
            int j;
            int k;
            //
            //  The image of each point is initially the image of the origin.
            //
            //  Insofar as the pre-image differs from the origin in a given vertex
            //  direction, add that proportion of the difference between the images
            //  of the origin and the vertex.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    phy[i + j * m] = t[i + 0 * m];

                    for (k = 1; k < m + 1; k++)
                    {
                        phy[i + j * m] = phy[i + j * m] + (t[i + k * m] - t[i + 0 * m]) * ref_[k - 1 + j * m];
                    }
                }
            }

        }

        public static double simplex_unit_volume(int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_VOLUME computes the volume of the unit simplex.
            //
            //  Discussion:
            //
            //    The formula is simple: volume = 1/M!.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the dimension of the space.
            //
            //    Output, double SIMPLEX_UNIT_VOLUME, the volume of the cone.
            //
        {
            int i;
            double volume;

            volume = 1.0;
            for (i = 1; i <= m; i++)
            {
                volume = volume / ((double)i);
            }

            return volume;
        }
    }
}