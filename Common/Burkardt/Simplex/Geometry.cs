using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.SimplexNS
{
    public static class Geometry
    {
        public static void simplex_lattice_layer_point_next(int n, int[] c, ref int[] v, ref bool more )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX_LATTICE_LAYER_POINT_NEXT: next simplex lattice layer point.
        //
        //  Discussion:
        //
        //    The simplex lattice layer L is bounded by the lines
        //
        //      0 <= X(1:N),
        //      L - 1 < sum X(1:N) / C(1:N)  <= L.
        //
        //    In particular, layer L = 0 always contains just the origin.
        //
        //    This function returns, one at a time, the points that lie within
        //    a given simplex lattice layer.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int C[N+1], coefficients defining the
        //    lattice layer in entries 1 to N, and the laver index in C[N].
        //    The coefficients should be positive, and C[N] must be nonnegative.
        //
        //    Input/output, int V[N].  On first call for a given layer,
        //    the input value of V is not important.  On a repeated call for the same
        //    layer, the input value of V should be the output value from the previous
        //    call.  On output, V contains the next lattice layer point.
        //
        //    Input/output, bool *MORE.  On input, set MORE to FALSE to indicate
        //    that this is the first call for a given layer.  Thereafter, the input
        //    value should be the output value from the previous call.  On output,
        //    MORE is TRUE if the returned value V is a new point.
        //    If the output value is FALSE, then no more points were found,
        //    and V was reset to 0, and the lattice layer has been exhausted.
        //
        {
            int c1n;
            int i;
            int j;
            int lhs;
            int rhs1;
            int rhs2;
            //
            //  Treat layer C[N] = 0 specially.
            //
            if (c[n] == 0)
            {
                if (!(more))
                {
                    for (j = 0; j < n; j++)
                    {
                        v[j] = 0;
                    }

                    more = true;
                }
                else
                {
                    more = false;
                }

                return;
            }

            //
            //  Compute the first point.
            //
            if (!(more))
            {
                v[0] = (c[n] - 1) * c[0] + 1;
                for (j = 1; j < n; j++)
                {
                    v[j] = 0;
                }

                more = true;
            }
            else
            {
                c1n = typeMethods.i4vec_lcm(n, c);

                rhs1 = c1n * (c[n] - 1);
                rhs2 = c1n * c[n];
                //
                //  Try to increment component I.
                //
                for (i = 0; i < n; i++)
                {
                    v[i] = v[i] + 1;

                    for (j = 0; j < i; j++)
                    {
                        v[j] = 0;
                    }

                    if (0 < i)
                    {
                        v[0] = rhs1;
                        for (j = 1; j < n; j++)
                        {
                            v[0] = v[0] - (c1n / c[j]) * v[j];
                        }

                        v[0] = (c[0] * v[0]) / c1n;
                        v[0] = Math.Max(v[0], 0);
                    }

                    lhs = 0;
                    for (j = 0; j < n; j++)
                    {
                        lhs = lhs + (c1n / c[j]) * v[j];
                    }

                    if (lhs <= rhs1)
                    {
                        v[0] = v[0] + 1;
                        lhs = lhs + c1n / c[0];
                    }

                    if (lhs <= rhs2)
                    {
                        return;
                    }
                }

                for (j = 0; j < n; j++)
                {
                    v[j] = 0;
                }

                more = false;
            }
        }

        public static void simplex_lattice_point_next(int n, int[] c, int[] v, ref bool more )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX_LATTICE_POINT_NEXT returns the next simplex lattice point.
        //
        //  Discussion:
        //
        //    The lattice simplex is defined by the vertices:
        //
        //      (0,0,...,0), (C[N]/C[0],0,...,0), (0,C[N]/C[1],...,0) ...
        //      (0,0,...C(N]/C[N-1])
        //
        //    The lattice simplex is bounded by the lines
        //
        //      0 <= V[0:N-1],
        //      V[0] / C[0] + V[1] / C[1] + ... + V[N-1] / C[N-1] <= C[N]
        //
        //    Lattice points are listed one at a time, starting at the origin,
        //    with V[0] increasing first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, int C[N+1], coefficients defining the
        //    lattice simplex.  These should be positive.
        //
        //    Input/output, int V[N].  On first call, the input
        //    value is not important.  On a repeated call, the input value should
        //    be the output value from the previous call.  On output, V contains
        //    the next lattice point.
        //
        //    Input/output, bool *MORE.  On input, set MORE to FALSE to indicate
        //    that this is the first call for a given simplex.  Thereafter, the input
        //    value should be the output value from the previous call.  On output,
        //    MORE is TRUE if not only is the returned value V a lattice point,
        //    but the routine can be called again for another lattice point.
        //    If the output value is FALSE, then no more lattice points were found,
        //    and V was reset to 0, and the routine should not be called further
        //    for this simplex.
        //
        {
            int c1n;
            int i;
            int j;
            int lhs;
            int rhs;
            int term;

            if (!(more))
            {
                typeMethods.i4vec_zero(n, ref v);
                more = true;
            }
            else
            {
                c1n = typeMethods.i4vec_lcm(n, c);
                rhs = c1n * c[n];

                lhs = 0;
                for (i = 0; i < n; i++)
                {
                    term = 1;
                    for (j = 0; j < n; j++)
                    {
                        if (i == j)
                        {
                            term = term * v[j];
                        }
                        else
                        {
                            term = term * c[j];
                        }
                    }

                    lhs = lhs + term;
                }

                for (i = 0; i < n; i++)
                {
                    if (lhs + c1n / c[i] <= rhs)
                    {
                        v[i] = v[i] + 1;
                        more = true;
                        return;
                    }

                    lhs = lhs - c1n * v[i] / c[i];
                    v[i] = 0;
                }

                more = false;
            }
        }

        public static int simplex_unit_lattice_point_nd(int d, int s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_LATTICE_POINT_ND: count lattice points.
            //
            //  Discussion:
            //
            //    The simplex is assumed to be the unit D-dimensional simplex:
            //
            //    ( (0,0,...,0), (1,0,...,0), (0,1,...,0), ... (0,,0,...,1) )
            //
            //    or a copy of this simplex scaled by an integer S:
            //
            //    ( (0,0,...,0), (S,0,...,0), (0,S,...,0), ... (0,,0,...,S) )
            //
            //    The routine returns the number of integer lattice points that appear
            //    inside the simplex or on its boundary.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 July 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Matthias Beck, Sinai Robins,
            //    Computing the Continuous Discretely,
            //    Springer, 2006,
            //    ISBN13: 978-0387291390,
            //    LC: QA640.7.B43.
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Input, int S, the scale factor.
            //
            //    Output, int SIMPLEX_UNIT_LATTICE_POINT_ND, the number of lattice points.
            //
        {
            int i;
            int n;

            n = 1;
            for (i = 1; i <= d; i++)
            {
                n = (n * (s + i)) / i;
            }

            return n;
        }

        public static double simplex_unit_volume_nd(int dim_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_VOLUME_ND computes the volume of the unit simplex in ND.
            //
            //  Discussion:
            //
            //    The formula is simple: volume = 1/DIM_NUM!.
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
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Output, double SIMPLEX_UNIT_VOLUME_ND, the volume of the cone.
            //
        {
            int i;
            double volume;

            volume = 1.0;
            for (i = 1; i <= dim_num; i++)
            {
                volume = volume / ((double) i);
            }

            return volume;
        }

        public static double simplex_volume_nd(int dim_num, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_VOLUME_ND computes the volume of a simplex in ND.
            //
            //  Discussion:
            //
            //    The formula is:
            //
            //      volume = 1/DIM_NUM! * det ( A )
            //
            //    where A is the DIM_NUM by DIM_NUM matrix obtained by subtracting one
            //    vector from all the others.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Input, double A[DIM_NUM*(DIM_NUM+1)], the points that define the simplex.
            //
            //    Output, double SIMPLEX_VOLUME_ND, the volume of the simplex.
            //
        {
            double[] b;
            double det;
            int i;
            int info;
            int j;
            int[] pivot;
            double volume;

            b = new double [dim_num * dim_num];
            pivot = new int [dim_num];

            for (j = 0; j < dim_num; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    b[i + j * dim_num] = a[i + j * dim_num] - a[i + dim_num * dim_num];
                }
            }

            info = Matrix.dge_fa(dim_num, ref b, ref pivot);

            if (info != 0)
            {
                volume = -1.0;
            }
            else
            {
                det = Matrix.dge_det(dim_num, b, pivot);

                volume = Math.Abs(det);
                for (i = 1; i <= dim_num; i++)
                {
                    volume = volume / ((double) i);
                }
            }
            
            return volume;
        }

    }
}