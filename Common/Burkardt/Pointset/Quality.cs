using System;
using Burkardt.MatrixNS;
using Burkardt.Sampling;
using Burkardt.Stroud;
using Burkardt.Types;

namespace Burkardt.Pointset;

public static class Quality
{
    public static double alpha_measure(int n, double[] z, int triangle_order, int triangle_num,
            int[] triangle_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ALPHA_MEASURE determines the triangulated pointset quality measure ALPHA.
        //
        //  Discusion:
        //
        //    The ALPHA measure evaluates the uniformity of the shapes of the triangles
        //    defined by a triangulated pointset.
        //
        //    We compute the minimum angle among all the triangles in the triangulated
        //    dataset and divide by the maximum possible value (which, in degrees,
        //    is 60).  The best possible value is 1, and the worst 0.  A good
        //    triangulation should have an ALPHA score close to 1.
        //
        //    The code has been modified to 'allow' 6-node triangulations.
        //    However, no effort is made to actually process the midside nodes.
        //    Only information from the vertices is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, real ( kind = 8 ) Z(2,N), the points.
        //
        //    Input, int TRIANGLE_ORDER, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
        //    the triangulation.
        //
        //    Output, double ALPHA_MEASURE, the ALPHA quality measure.
        //
    {
        int triangle;

        double alpha = typeMethods.r8_huge();

        for (triangle = 0; triangle < triangle_num; triangle++)
        {
            int a_index = triangle_node[0 + triangle * 3];
            int b_index = triangle_node[1 + triangle * 3];
            int c_index = triangle_node[2 + triangle * 3];

            double a_x = z[0 + (a_index - 1) * 2];
            double a_y = z[1 + (a_index - 1) * 2];
            double b_x = z[0 + (b_index - 1) * 2];
            double b_y = z[1 + (b_index - 1) * 2];
            double c_x = z[0 + (c_index - 1) * 2];
            double c_y = z[1 + (c_index - 1) * 2];

            double ab_len = Math.Sqrt(Math.Pow(a_x - b_x, 2) + Math.Pow(a_y - b_y, 2));
            double bc_len = Math.Sqrt(Math.Pow(b_x - c_x, 2) + Math.Pow(b_y - c_y, 2));
            double ca_len = Math.Sqrt(Math.Pow(c_x - a_x, 2) + Math.Pow(c_y - a_y, 2));
            double c_angle;
            double a_angle;
            double b_angle;
            switch (ab_len)
            {
                //
                //  Take care of a ridiculous special case.
                //
                case 0.0 when bc_len == 0.0 && ca_len == 0.0:
                    a_angle = 2.0 * Math.PI / 3.0;
                    b_angle = 2.0 * Math.PI / 3.0;
                    c_angle = 2.0 * Math.PI / 3.0;
                    break;
                default:
                {
                    if (ca_len == 0.0 || ab_len == 0.0)
                    {
                        a_angle = Math.PI;
                    }
                    else
                    {
                        a_angle = Helpers.arc_cosine(
                            (ca_len * ca_len + ab_len * ab_len - bc_len * bc_len)
                            / (2.0 * ca_len * ab_len));
                    }

                    if (ab_len == 0.0 || bc_len == 0.0)
                    {
                        b_angle = Math.PI;
                    }
                    else
                    {
                        b_angle = Helpers.arc_cosine(
                            (ab_len * ab_len + bc_len * bc_len - ca_len * ca_len)
                            / (2.0 * ab_len * bc_len));
                    }

                    if (bc_len == 0.0 || ca_len == 0.0)
                    {
                        c_angle = Math.PI;
                    }
                    else
                    {
                        c_angle = Helpers.arc_cosine(
                            (bc_len * bc_len + ca_len * ca_len - ab_len * ab_len)
                            / (2.0 * bc_len * ca_len));
                    }

                    break;
                }
            }

            alpha = Math.Min(alpha, a_angle);
            alpha = Math.Min(alpha, b_angle);
            alpha = Math.Min(alpha, c_angle);
        }

        //
        //  Normalize angles from [0,60] degrees into qualities in [0,1].
        //
        double value = alpha * 3.0 / Math.PI;

        return value;
    }
        
    public static double area_measure ( int n, double[] z, int triangle_order, int triangle_num,
            int[] triangle_node )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AREA_MEASURE determines the area ratio quality measure.
        //
        //  Discusion:
        //
        //    This measure computes the area of every triangle, and returns
        //    the ratio of the minimum to the maximum triangle.  A value of
        //    1 is "perfect", indicating that all triangles have the same area.
        //    A value of 0 is the worst possible result.
        //
        //    The code has been modified to 'allow' 6-node triangulations.
        //    However, no effort is made to actually process the midside nodes.
        //    Only information from the vertices is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[2*N], the points.
        //
        //    Input, int TRIANGLE_ORDER, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
        //    the triangulation.
        //
        //    Output, double AREA_MEASURE, the AREA quality measure.
        //
    {
        int triangle;

        double area_max = 0.0;
        double area_min = typeMethods.r8_huge();

        for (triangle = 0; triangle < triangle_num; triangle++)
        {
            double x1 = z[0 + (triangle_node[0 + triangle * 3] - 1) * 2];
            double y1 = z[1 + (triangle_node[0 + triangle * 3] - 1) * 2];
            double x2 = z[0 + (triangle_node[1 + triangle * 3] - 1) * 2];
            double y2 = z[1 + (triangle_node[1 + triangle * 3] - 1) * 2];
            double x3 = z[0 + (triangle_node[2 + triangle * 3] - 1) * 2];
            double y3 = z[1 + (triangle_node[2 + triangle * 3] - 1) * 2];

            double area = 0.5 * Math.Abs(x1 * (y2 - y3)
                                         + x2 * (y3 - y1)
                                         + x3 * (y1 - y2));

            area_min = Math.Min(area_min, area);
            area_max = Math.Max(area_max, area);
        }

        double value = area_max switch
        {
            > 0.0 => area_min / area_max,
            _ => 0.0
        };

        return value;
    }

    public static double beta_measure(int dim_num, int n, double[] z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_MEASURE determines the pointset quality measure BETA.
        //
        //  Discussion:
        //
        //    The BETA measure of point distribution quality for a set Z of
        //    N points in an DIM_NUM dimensional region is defined as follows:
        //
        //    For each point Z(I), determine the nearest distinct element of
        //    the pointset by
        //
        //      GAMMA(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
        //
        //    Let GAMMA_AVE be the average of GAMMA(1:N).
        //
        //    Let GAMMA_STD be the standard deviation of the GAMMA's:
        //
        //      GAMMA_STD = Math.Sqrt ( 1 / ( N - 1 )
        //        * sum ( 1 <= I <= N ) ( GAMMA(I) - GAMMA_AVE )**2 ) )
        //
        //    Then BETA is the standard deviation normalized by the average:
        //
        //      BETA = GAMMA_STD / GAMMA_AVE.
        //
        //    For an ideally regular mesh, the GAMMA(I)'s will be equal and
        //    BETA will be zero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Output, double BETA_MEASURE, the BETA quality measure.
        //
    {
        double gamma_std;
        int i;

        double[] gamma = Spacing.pointset_spacing(dim_num, n, z);

        double gamma_ave = 0.0;
        for (i = 0; i < n; i++)
        {
            gamma_ave += gamma[i];
        }

        gamma_ave /= n;

        switch (n)
        {
            case > 1:
            {
                gamma_std = 0.0;
                for (i = 0; i < n; i++)
                {
                    gamma_std += Math.Pow(gamma[i] - gamma_ave, 2);
                }

                gamma_std = Math.Sqrt(gamma_std / (n - 1));
                break;
            }
            default:
                gamma_std = 0.0;
                break;
        }

        double value = gamma_ave switch
        {
            > 0.0 => gamma_std / gamma_ave,
            _ => 0.0
        };

        return value;
    }

    public static double chi_measure(int dim_num, int n, double[] z, int ns,
            Func<int, int, int, GeometrySampleResult> sample_routine,
            int seed_init)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_MEASURE determines the pointset quality measure CHI.
        //
        //  Discussion:
        //
        //    The CHI measure of point distribution quality for a set Z of
        //    N points in an DIM_NUM-dimensional region is defined as follows:
        //
        //    Assign every point X in the region to the nearest element
        //    Z(I) of the point set.  For each Z(I), let H(I) be the maximum
        //    distance between Z(I) and any point assigned to it by this process.
        //
        //    For each point Z(I), we determine the nearest distinct element of
        //    the pointset by
        //
        //      GAMMA(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
        //
        //    Then
        //
        //      CHI(I) = 2 * H(I) / GAMMA(I)
        //
        //    and
        //
        //      CHI = maximum ( 1 <= I <= N ) CHI(I)
        //
        //    This quantity can be estimated by using sampling to pick a large
        //    number of points in the region, rather than all of them.
        //
        //    For an ideally regular mesh, all the CHI(I)'s will be equal.
        //    Any deviation from regularity increases the value of some entries
        //    of CHI; thus, given two meshes, the one with a lower value of
        //    CHI is to be recommended.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Input, int NS, the number of sample points.
        //
        //    Input, double *SAMPLE_ROUTINE, the name of a routine which
        //    is used to produce an DIM_NUM by N array of sample points in the region,
        //    of the form:
        //      double *sample_routine ( int dim_num, int n, int *seed )
        //
        //    Input, int SEED_INIT, the initial value of the random number seed.
        //
        //    Output, double CHI_MEASURE, the CHI quality measure.
        //
    {
        int[] closest = new int[1];
        int j;
        int k;

        int seed = seed_init;

        double[] chi_vec = new double[n];
        double[] h = new double[n];

        for (j = 0; j < n; j++)
        {
            h[j] = 0.0;
        }

        for (k = 1; k <= ns; k++)
        {
            GeometrySampleResult result = sample_routine(dim_num, 1, seed);
            double[] x = result.result;
            seed = result.seed;

            find_closest(dim_num, n, 1, x, z, closest);

            double dist = 0.0;
            int i;
            for (i = 0; i < dim_num; i++)
            {
                dist += Math.Pow(x[i] - z[i + closest[0] * dim_num], 2);
            }

            h[closest[0]] = Math.Max(h[closest[0]], dist);

        }

        double[] gamma = Spacing.pointset_spacing(dim_num, n, z);

        double chi = 0.0;

        for (j = 0; j < n; j++)
        {
            chi_vec[j] = 2.0 * Math.Sqrt(h[j]) / gamma[j];
            chi = Math.Max(chi, chi_vec[j]);
        }

        return chi;
    }

    public static double d_measure(int dim_num, int n, double[] z, int ns,
            Func<int, int, int, GeometrySampleResult> sample_routine,
            int seed_init)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    D_MEASURE determines the pointset quality measure D.
        //
        //  Discussion:
        //
        //    The D measure of point distribution quality for a set Z of
        //    N points in an DIM_NUM-dimensional region is defined as follows:
        //
        //    For each point Z(I) in the pointset, let V(I) be the subregion
        //    defined by the intersection of the region with the Voronoi
        //    region associated with Z(I).
        //
        //    Let D(I) be the determinant of the deviatoric tensor associated with
        //    the region V(I).
        //
        //    Then D = maximum ( 1 <= I <= N ) D(I).
        //
        //    This quantity can be estimated using sampling.  A given number of
        //    sample points are generated in the region, assigned to the nearest
        //    element of the pointset, and used to approximate the Voronoi regions
        //    and the deviatoric tensors.
        //
        //    In an ideally regular mesh, each deviatoric tensor would have a
        //    zero determinant, and hence D would be zero.  In general, the smaller
        //    D, the better.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Input, int NS, the number of sample points.
        //
        //    Input, double *SAMPLE_ROUTINE, the name of a routine which
        //    is used to produce an DIM_NUM by N array of sample points in the region,
        //    of the form:
        //      double *sample_routine ( int dim_num, int n, int *seed )
        //
        //    Input, int SEED_INIT, the initial value of the random number seed.
        //
        //    Output, double D_MEASURE, the D quality measure.
        //
    {
        int[] closest = new int[1];
        int i;
        int i1;
        int i2;
        int j;
        int k;

        double[] a = new double[dim_num * dim_num];
        double[] centroid = new double[dim_num * n];
        int[] hit = new int[n];
        double[] moment = new double[dim_num * dim_num * n];
        double[] tri = new double[n];

        int seed = seed_init;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < dim_num; i++)
            {
                centroid[i + j * dim_num] = 0.0;
            }
        }

        for (j = 0; j < n; j++)
        {
            hit[j] = 0;
        }

        for (j = 0; j < n; j++)
        {
            for (i2 = 0; i2 < dim_num; i2++)
            {
                for (i1 = 0; i1 < dim_num; i1++)
                {
                    moment[i1 + i2 * dim_num + j * dim_num * dim_num] = 0.0;
                }
            }
        }

        for (k = 1; k <= ns; k++)
        {
            GeometrySampleResult result = sample_routine(dim_num, 1, seed);
            double[] x = result.result;
            seed = result.seed;

            find_closest(dim_num, n, 1, x, z, closest);

            hit[closest[0]] += 1;

            for (i = 0; i < dim_num; i++)
            {
                centroid[i + closest[0] * dim_num] += x[i];
            }

            for (i1 = 0; i1 < dim_num; i1++)
            {
                for (i2 = 0; i2 < dim_num; i2++)
                {
                    moment[i1 + i2 * dim_num + closest[0] * dim_num * dim_num] += x[i1] * x[i2];
                }
            }
        }

        for (j = 0; j < n; j++)
        {
            switch (hit[j])
            {
                case > 0:
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        centroid[i + j * dim_num] /= hit[j];
                    }

                    for (i1 = 0; i1 < dim_num; i1++)
                    {
                        for (i2 = 0; i2 < dim_num; i2++)
                        {
                            moment[i1 + i2 * dim_num + j * dim_num * dim_num] /= hit[j];
                        }
                    }

                    for (i1 = 0; i1 < dim_num; i1++)
                    {
                        for (i2 = 0; i2 < dim_num; i2++)
                        {
                            moment[i1 + i2 * dim_num + j * dim_num * dim_num] -= centroid[i1 + j * dim_num] * centroid[i2 + j * dim_num];
                        }
                    }

                    break;
                }
            }
        }

        for (j = 0; j < n; j++)
        {
            tri[j] = 0.0;
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < dim_num; i++)
            {
                tri[j] += moment[i + i * dim_num + j * dim_num * dim_num];
            }
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < dim_num; i++)
            {
                moment[i + i * dim_num + j * dim_num * dim_num] -= tri[j] / dim_num;
            }
        }

        double d = 0.0;

        for (j = 0; j < n; j++)
        {
            for (i2 = 0; i2 < dim_num; i2++)
            {
                for (i1 = 0; i1 < dim_num; i1++)
                {
                    a[i1 + i2 * dim_num] = moment[i1 + i2 * dim_num + j * dim_num * dim_num];
                }
            }

            double di = Matrix.dge_det(dim_num, ref a);

            d = Math.Max(d, di);

        }

        return d;
    }

    public static double e_measure(int dim_num, int n, double[] z, int ns,
            Func<int, int, int, GeometrySampleResult> sample_routine,
            int seed_init)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    E_MEASURE determines the pointset quality measure E.
        //
        //  Discussion:
        //
        //    The E measure of point distribution quality for a set Z of
        //    N points in an DIM_NUM dimensional region is defined as follows:
        //
        //    Assign every point X in the region to the nearest element
        //    Z(I) of the point set.  For each point Z(I), let E_VEC(I) be the
        //    integral of the distance between Z(I) and all the points assigned to
        //    it:
        //
        //      E_VEC(I) = Integral ( all X nearest to Z(I) ) distance ( X, Z(I) )
        //
        //    If we let VOLUME be the volume of the region, then we define E by:
        //
        //      E = sum ( 1 <= I <= N ) E_VEC(I) / VOLUME
        //
        //    This quantity can be estimated by using sampling to pick a large
        //    number of points in the region, rather than all of them.
        //
        //    The E measure is minimized by a centroidal Voronoi tessellation.
        //
        //    Given two meshes, the one with a lower value of E is to be recommended.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Input, int NS, the number of sample points.
        //
        //    Input, double *SAMPLE_ROUTINE, the name of a routine which
        //    is used to produce an DIM_NUM by N array of sample points in the region,
        //    of the form:
        //      double *sample_routine ( int dim_num, int n, int *seed )
        //
        //    Input, int SEED_INIT, the initial value of the random number seed.
        //
        //    Output, double E_MEASURE, the E quality measure.
        //
    {
        int[] closest = new int[1];
        int j;
        int k;

        int seed = seed_init;

        double[] e_vec = new double[n];

        for (j = 0; j < n; j++)
        {
            e_vec[j] = 0.0;
        }

        for (k = 1; k <= ns; k++)
        {
            GeometrySampleResult result = sample_routine(dim_num, 1, seed);
            double[] x = result.result;
            seed = result.seed;

            find_closest(dim_num, n, 1, x, z, closest);

            double dist = 0.0;
            int i;
            for (i = 0; i < dim_num; i++)
            {
                dist += Math.Pow(x[i] - z[i + closest[0] * dim_num], 2);
            }

            e_vec[closest[0]] += dist;

        }

        double e = 0.0;
        for (j = 0; j < n; j++)
        {
            e += e_vec[j];
        }

        e /= ns;

        return e;
    }

    public static void find_closest ( int dim_num, int n, int sample_num, double[] s, double[] r,
            int[] nearest )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FIND_CLOSEST finds the nearest R point to each S point.
        //
        //  Discussion:
        //
        //    This routine finds the closest Voronoi cell generator by checking every
        //    one.  For problems with many cells, this process can take the bulk
        //    of the CPU time.  Other approaches, which group the cell generators into
        //    bins, can run faster by a large factor.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of cell generators.
        //
        //    Input, int SAMPLE_NUM, the number of sample points.
        //
        //    Input, double S[DIM_NUM*SAMPLE_NUM], the points to be checked.
        //
        //    Input, double R[DIM_NUM*N], the cell generators.
        //
        //    Output, int NEAREST[SAMPLE_NUM], the (0-based) index of the nearest
        //    cell generator.
        //
    {
        int js;

        for ( js = 0; js < sample_num; js++ )
        {
            double dist_sq_min = typeMethods.r8_huge ( );
            nearest[js] = -1;

            int jr;
            for ( jr = 0; jr < n; jr++ )
            {
                double dist_sq = 0.0;
                int i;
                for ( i = 0; i < dim_num; i++ )
                {
                    dist_sq += ( s[i+js*dim_num] - r[i+jr*dim_num] )
                               * ( s[i+js*dim_num] - r[i+jr*dim_num] );
                }

                if ( jr == 0 || dist_sq < dist_sq_min )
                {
                    dist_sq_min = dist_sq;
                    nearest[js] = jr;
                }
            }
        }
    }

    public static double gamma_measure(int dim_num, int n, double[] z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_MEASURE determines the pointset quality measure GAMMA.
        //
        //  Discussion:
        //
        //    The GAMMA measure of point distribution quality for a set Z of
        //    N points in an DIM_NUM-dimensional region is defined as follows:
        //
        //      GAMMA = ( GAMMA_MAX / GAMMA_MIN ),
        //
        //    where
        //
        //      GAMMA_MAX = maximum ( 1 <= I <= N ) DIST_MIN(I)
        //      GAMMA_MIN = minimum ( 1 <= I <= N ) DIST_MIN(I)
        //
        //    and
        //
        //      DIST_MIN(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
        //
        //
        //    Note that, in this code, the variable DIST_SQ_MIN is actually the square
        //    of the minimum point distance, and so when we compute GAMMA, we must
        //    take the square root of the given ratio.
        //
        //    GAMMA must be at least 1.  For an ideally regular mesh, GAMMA would
        //    be equal to one.  Given two meshes, this measure recommends the one
        //    with the smaller value of GAMMA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Output, double GAMMA_MEASURE, the GAMMA quality measure.
        //
        //  Local parameters:
        //
        //    Local, double GAMMA_SQ_MAX, the maximum, over all points,
        //    of the minimum squared distance to a distinct point.
        //
        //    Local, double GAMMA_SQ_MIN, the minimum, over all points,
        //    of the minimum squared distance to a distinct point.
        //
    {
        int j1;
        double gamma;
        switch (n)
        {
            //
            //  Take care of ridiculous cases.
            //
            case <= 1:
                gamma = 0.0;
                return gamma;
        }

        double gamma_sq_max = 0.0;
        double gamma_sq_min = typeMethods.r8_huge();

        for (j1 = 0; j1 < n; j1++)
        {
            double dist_sq_min = typeMethods.r8_huge();

            int j2;
            for (j2 = 0; j2 < n; j2++)
            {
                if (j2 != j1)
                {

                    double dist_sq = 0.0;
                    int i;
                    for (i = 0; i < dim_num; i++)
                    {
                        dist_sq += Math.Pow(z[i + j1 * dim_num] - z[i + j2 * dim_num], 2);
                    }

                    if (dist_sq < dist_sq_min)
                    {
                        dist_sq_min = dist_sq;
                    }
                }

            }

            gamma_sq_max = Math.Max(gamma_sq_max, dist_sq_min);
            gamma_sq_min = Math.Min(gamma_sq_min, dist_sq_min);

        }

        gamma = gamma_sq_min switch
        {
            <= 0.0 => typeMethods.r8_huge(),
            _ => Math.Sqrt(gamma_sq_max / gamma_sq_min)
        };

        return gamma;
    }

    public static double h_measure(int dim_num, int n, double[] z, int ns,
            Func<int, int, int, GeometrySampleResult> sample_routine,
            int seed_init)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    H_MEASURE determines the pointset quality measure H.
        //
        //  Discussion:
        //
        //    The H measure of dispersion for a set of N points in an DIM_NUM-dimensional
        //    region is the maximum distance between a point in the region and some
        //    point in the set.
        //
        //    To compute this quantity exactly, for every point X in the region,
        //    find the nearest element Z of the point set and compute the distance.
        //    H is then the maximum of all these distances.
        //
        //    To ESTIMATE this quantity, carry out the same process, but only for
        //    NS sample points in the region.
        //
        //    Under this measure, a mesh with a smaller value of H is preferable.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Input, int NS, the number of sample points.
        //
        //    Input, double *SAMPLE_ROUTINE, the name of a routine which
        //    is used to produce an DIM_NUM by N array of sample points in the region,
        //    of the form:
        //      double *sample_routine ( int dim_num, int n, int *seed )
        //
        //    Input, int SEED_INIT, the initial value of the random number seed.
        //
        //    Output, double H_MEASURE, the H quality measure.
        //
    {
        int[] closest = new int[1];
        int k;

        int seed = seed_init;
        double h = 0.0;

        for (k = 1; k <= ns; k++)
        {
            GeometrySampleResult result = sample_routine(dim_num, 1, seed);
            double[] x = result.result;
            seed = result.seed;

            find_closest(dim_num, n, 1, x, z, closest);

            double dist = 0.0;
            int i;
            for (i = 0; i < dim_num; i++)
            {
                dist += Math.Pow(x[i] - z[i + closest[0] * dim_num], 2);
            }

            h = Math.Max(h, dist);

        }

        h = Math.Sqrt(h);

        return h;
    }

    public static double lambda_measure(int dim_num, int n, double[] z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAMBDA_MEASURE determines the pointset quality measure LAMBDA.
        //
        //  Discussion:
        //
        //    The LAMBDA measure of point distribution quality for a set Z of
        //    N points in an DIM_NUM-dimensional region is defined as follows:
        //
        //    Let
        //
        //      GAMMA(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
        //
        //    and let
        //
        //      GAMMA_AVE = sum ( 1 <= I <= N ) GAMMA(I) / N
        //
        //    then
        //
        //      LAMBDA = Math.Sqrt ( sum ( 1 <= I <= N ) ( GAMMA(I) - GAMMA_AVE )**2 / N )
        //        / GAMMA_AVE
        //
        //    An ideally regular mesh would have GAMMA(I) = GAMMA_AVE for all I,
        //    so that LAMBDA would be 0.  Under this measure, the mesh with the
        //    smaller value of LAMBDA is to be preferred.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Output, double LAMBDA_MEASURE, the LAMBDA quality measure.
        //
        //  Local parameters:
        //
        //    Local, double GAMMA_MAX, the maximum, over all points,
        //    of the minimum distance to a distinct point.
        //
        //    Local, double GAMMA_MIN, the minimum, over all points,
        //    of the minimum distance to a distinct point.
        //
    {
        int j;
        double lambda;
        switch (n)
        {
            //
            //  Take care of ridiculous cases.
            //
            case <= 1:
                lambda = 0.0;
                return lambda;
        }

        //
        //  Compute the minimum spacing between distinct points of the set.
        //
        double[] gamma = Spacing.pointset_spacing(dim_num, n, z);
        //
        //  Average the minimum spacing.
        //
        double gamma_ave = 0.0;
        for (j = 0; j < n; j++)
        {
            gamma_ave += gamma[j];
        }

        gamma_ave /= n;
        switch (gamma_ave)
        {
            //
            //  Compute a weighted variance.
            //
            case <= 0.0:
                lambda = typeMethods.r8_huge();
                break;
            default:
            {
                lambda = 0.0;
                for (j = 0; j < n; j++)
                {
                    lambda += Math.Pow(gamma[j] - gamma_ave, 2);
                }

                lambda = Math.Sqrt(lambda / n);
                lambda /= gamma_ave;
                break;
            }
        }

        return lambda;
    }

    public static double mu_measure(int dim_num, int n, double[] z, int ns,
            Func<int, int, int, GeometrySampleResult> sample_routine,
            int seed_init)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MU_MEASURE determines the pointset quality measure MU.
        //
        //  Discussion:
        //
        //    The MU measure of dispersion for a set of N points in an DIM_NUM-dimensional
        //    region takes the ratio of the largest and smallest half-diameters
        //    of the Voronoi cells defined by a pointset.
        //
        //    To compute this quantity exactly, for every point X in the region,
        //    find the nearest element Z of the point set and compute the distance.
        //
        //    Then, for each element Z(I) of the point set, define H(I) to be the
        //    maximum of these distances.
        //
        //    MU is then the ratio of the maximum and minimum values of H.
        //
        //    To ESTIMATE this quantity, carry out the same process, but only for
        //    NS sample points in the region.
        //
        //    In an ideally regular mesh, MU would be 1.  MU must be at least 1.
        //    Under this measure, the mesh with the smaller value of MU is to be
        //    preferred.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Input, int NS, the number of sample points.
        //
        //    Input, double *SAMPLE_ROUTINE, the name of a routine which
        //    is used to produce an DIM_NUM by N array of sample points in the region,
        //    of the form:
        //      double *sample_routine ( int dim_num, int n, int *seed )
        //
        //    Input, int SEED_INIT, the initial value of the random number seed.
        //
        //    Output, double MU_MEASURE, the MU quality measure.
        //
    {
        int k;
        int[] closest = new int[1];
        int j;

        double[] h = new double[n];

        int seed = seed_init;

        for (j = 0; j < n; j++)
        {
            h[j] = 0.0;
        }

        for (k = 1; k <= ns; k++)
        {
            GeometrySampleResult result = sample_routine(dim_num, 1, seed);
            double[] x = result.result;
            seed = result.seed;

            find_closest(dim_num, n, 1, x, z, closest);

            double dist = 0.0;
            int i;
            for (i = 0; i < dim_num; i++)
            {
                dist += Math.Pow(x[i] - z[i + closest[0] * dim_num], 2);
            }

            h[closest[0]] = Math.Max(h[closest[0]], dist);

        }

        double h_max = h[0];
        for (j = 1; j < n; j++)
        {
            h_max = Math.Max(h_max, h[j]);
        }

        h_max = Math.Sqrt(h_max);

        double h_min = h[0];
        for (j = 1; j < n; j++)
        {
            h_min = Math.Min(h_min, h[j]);
        }

        h_min = Math.Sqrt(h_min);

        double mu = h_min switch
        {
            0.0 => typeMethods.r8_huge(),
            _ => h_max / h_min
        };

        return mu;
    }

    public static double nu_measure(int dim_num, int n, double[] z, int ns,
            Func<int, int, int, GeometrySampleResult> sample_routine,
            int seed_init)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NU_MEASURE determines the pointset quality measure NU.
        //
        //  Discussion:
        //
        //    The NU measure of dispersion for a set of N points in an DIM_NUM-dimensional
        //    region is defined as follows:
        //
        //    For each element Z(I) of the pointset, let VOLUME(I) be the volume
        //    of the corresponding Voronoi subregion, restricted to the region.
        //
        //    Then
        //
        //      NU = max ( 1 <= I <= N ) VOLUME(I) / min ( 1 <= I <= N ) VOLUME(I)
        //
        //    This quantity can be estimated by using a large number of sampling
        //    points to estimate the Voronoi volumes.
        //
        //    For an ideally uniform pointset, the Voronoi volumes would be equal,
        //    so that NU would be 1.  In any case, NU must be 1 or greater.  In
        //    comparing two meshes, the one with the lower value of NU would be
        //    preferred.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Input, int NS, the number of sample points.
        //
        //    Input, double *SAMPLE_ROUTINE, the name of a routine which
        //    is used to produce an DIM_NUM by N array of sample points in the region,
        //    of the form:
        //      double *sample_routine ( int dim_num, int n, int *seed )
        //
        //    Input, int SEED_INIT, the initial value of the random number seed.
        //
        //    Output, double NU_MEASURE, the NU quality measure.
        //
    {
        int[] closest = new int[1];
        int j;
        int k;

        int[] hit = new int[n];
        double[] volume = new double[n];

        int seed = seed_init;

        for (j = 0; j < n; j++)
        {
            hit[j] = 0;
        }

        for (k = 1; k <= ns; k++)
        {
            GeometrySampleResult result = sample_routine(dim_num, 1, seed);
            double[] x = result.result;
            seed = result.seed;

            find_closest(dim_num, n, 1, x, z, closest);

            hit[closest[0]] += 1;

        }

        for (j = 0; j < n; j++)
        {
            volume[j] = hit[j] / (double) ns;
        }

        double volume_max = 0.0;
        for (j = 0; j < n; j++)
        {
            volume_max = Math.Max(volume_max, volume[j]);
        }

        double volume_min = typeMethods.r8_huge();
        for (j = 0; j < n; j++)
        {
            volume_min = Math.Min(volume_min, volume[j]);
        }

        double nu = volume_min switch
        {
            0.0 => typeMethods.r8_huge(),
            _ => volume_max / volume_min
        };

        return nu;
    }

    public static double q_measure(int n, double[] z, int triangle_order, int triangle_num,
            int[] triangle_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Q_MEASURE determines the triangulated pointset quality measure Q.
        //
        //  Discussion:
        //
        //    The Q measure evaluates the uniformity of the shapes of the triangles
        //    defined by a triangulated pointset.
        //
        //    For a single triangle T, the value of Q(T) is defined as follows:
        //
        //      TAU_IN = radius of the inscribed circle,
        //      TAU_OUT = radius of the circumscribed circle,
        //
        //      Q(T) = 2 * TAU_IN / TAU_OUT
        //        = ( B + C - A ) * ( C + A - B ) * ( A + B - C ) / ( A * B * C )
        //
        //    where A, B and C are the lengths of the sides of the triangle T.
        //
        //    The Q measure computes the value of Q(T) for every triangle T in the
        //    triangulation, and then computes the minimum of this
        //    set of values:
        //
        //      Q_MEASURE = min ( all T in triangulation ) Q(T)
        //
        //    In an ideally regular mesh, all triangles would have the same
        //    equilateral shape, for which Q = 1.  A good mesh would have
        //    0.5 < Q.
        //
        //    Given the 2D coordinates of a set of N nodes, stored as Z(1:2,1:N),
        //    a triangulation is a list of TRIANGLE_NUM triples of node indices that form
        //    triangles.  Generally, a maximal triangulation is expected, namely,
        //    a triangulation whose image is a planar graph, but for which the
        //    addition of any new triangle would mean the graph was no longer planar.
        //    A Delaunay triangulation is a maximal triangulation which maximizes
        //    the minimum angle that occurs in any triangle.
        //
        //    The code has been modified to 'allow' 6-node triangulations.
        //    However, no effort is made to actually process the midside nodes.
        //    Only information from the vertices is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //    Per-Olof Persson and Gilbert Strang,
        //    A Simple Mesh Generator in MATLAB,
        //    SIAM Review,
        //    Volume 46, Number 2, pages 329-345, June 2004.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[2*N], the points.
        //
        //    Input, int TRIANGLE_ORDER, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
        //    the triangulation.
        //
        //    Output, double Q_MEASURE, the Q quality measure.
        //
    {
        int triangle;
        double value;

        switch (triangle_num)
        {
            case < 1:
                value = -1.0;
                return value;
        }

        double q_min = typeMethods.r8_huge();

        for (triangle = 0; triangle < triangle_num; triangle++)
        {
            int a_index = triangle_node[0 + triangle * 3];
            int b_index = triangle_node[1 + triangle * 3];
            int c_index = triangle_node[2 + triangle * 3];

            double ab_length = Math.Sqrt(
                Math.Pow(z[0 + (a_index - 1) * 2] - z[0 + (b_index - 1) * 2], 2)
                + Math.Pow(z[1 + (a_index - 1) * 2] - z[1 + (b_index - 1) * 2], 2));

            double bc_length = Math.Sqrt(
                Math.Pow(z[0 + (b_index - 1) * 2] - z[0 + (c_index - 1) * 2], 2)
                + Math.Pow(z[1 + (b_index - 1) * 2] - z[1 + (c_index - 1) * 2], 2));

            double ca_length = Math.Sqrt(
                Math.Pow(z[0 + (c_index - 1) * 2] - z[0 + (a_index - 1) * 2], 2)
                + Math.Pow(z[1 + (c_index - 1) * 2] - z[1 + (a_index - 1) * 2], 2));

            double q = (bc_length + ca_length - ab_length)
                       * (ca_length + ab_length - bc_length)
                       * (ab_length + bc_length - ca_length)
                       / (ab_length * bc_length * ca_length);

            q_min = Math.Min(q_min, q);
        }

        value = q_min;

        return value;
    }

    public static double r0_measure(int dim_num, int n, double[] z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R0_MEASURE determines the pointset quality measure R0.
        //
        //  Discussion:
        //
        //    The R0 measure of point distribution quality for a set Z of
        //    N points in an DIM_NUM-dimensional region is defined as follows:
        //
        //      R0 = sum ( 1 <= I /= J <= N ) log ( 1 / distance ( Z(I), Z(J) ) )
        //         / ( N * ( N - 1 ) )
        //
        //    The divisor of ( N * ( N - 1 ) ) means that R0 is essentially an
        //
        //    R0 is undefined if N < 2 or if any two points are equal.
        //
        //    R0 is known as the Riesz S-energy for S = 0.
        //
        //    Given two meshes, this measure recommends the one with the smaller
        //    value of R0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    D P Hardin and E B Saff,
        //    Discretizing Manifolds via Minimum Energy Points,
        //    Notices of the AMS,
        //    Volume 51, Number 10, November 2004, pages 1186-1194.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Output, double R0_MEASURE, the R0 quality measure.
        //
    {
        int j1;
        double value;
        switch (n)
        {
            //
            //  Take care of ridiculous cases.
            //
            case <= 1:
                value = typeMethods.r8_huge();
                return value;
        }

        value = 0.0;

        for (j1 = 0; j1 < n; j1++)
        {
            int j2;
            for (j2 = 0; j2 < n; j2++)
            {
                if (j2 != j1)
                {
                    double dist = 0.0;
                    int i;
                    for (i = 0; i < dim_num; i++)
                    {
                        dist += Math.Pow(z[i + j1 * dim_num] - z[i + j2 * dim_num], 2);
                    }

                    dist = Math.Sqrt(dist);

                    switch (dist)
                    {
                        case 0.0:
                            value = typeMethods.r8_huge();
                            return value;
                        default:
                            value += Math.Log(1.0 / dist);
                            break;
                    }
                }
            }
        }

        value /= n * (n - 1);

        return value;
    }

    public static double[] radius_maximus(int dim_num, int n, double[] z, bool walls)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RADIUS_MAXIMUS finds the biggest possible nonintersecting sphere.
        //
        //  Discussion:
        //
        //    We are given a set of N points in DIM_NUM space.  We imagine that
        //    at each point simultaneously, a sphere begins to expand.
        //    Each sphere stops expanding as soon as it touches another sphere.
        //    The radius of these spheres is to be computed.
        //
        //    If WALLS is true, then the spheres must not extend outside the
        //    "walls" of the unit hypersquare.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of spatial dimensions.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the point coordinates.
        //    If WALLS is TRUE, these values must be between 0 and 1.
        //
        //    Input, logical WALLS, is TRUE if the spheres must not extend
        //    outside the unit hypercube.  If WALLS is FALSE, then this
        //    restriction is not imposed.
        //
        //    Output, double RADIUS(N), the radius of the
        //    maximal nonintersecting sphere around each point.
        //
    {
        const int FIXED = 0;
        const int FREE = 1;
        int i;
        int j;

        double[] radius = new double[n];
        int[] status = new int[n];

        switch (walls)
        {
            case true:
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        switch (z[i + j * dim_num])
                        {
                            case < 0.0:
                                Console.WriteLine("");
                                Console.WriteLine("RADIUS_MAXIMUS - Fatal error!");
                                Console.WriteLine("  Some coordinate is less than 0.");
                                return null;
                            case > 1.0:
                                Console.WriteLine("");
                                Console.WriteLine("RADIUS_MAXIMUS - Fatal error!");
                                Console.WriteLine("  Some coordinate is greater than 1.");
                                return null;
                        }
                    }
                }

                break;
            }
        }

        //
        //  Initially, all points are "free".
        //
        for (j = 0; j < n; j++)
        {
            radius[j] = 0.0;
        }

        for (j = 0; j < n; j++)
        {
            status[j] = FREE;
        }

        for (;;)
        {
            //
            //  If all points are fixed, we're done.
            //
            bool done = true;

            for (j = 0; j < n; j++)
            {
                if (status[j] == FIXED)
                {
                    continue;
                }

                done = false;
                break;
            }

            if (done)
            {
                break;
            }

            //
            //  Look at all the free points.
            //  Imagine an expanding sphere at each free point, and determine
            //  which such sphere will first have to stop expanding.
            //
            int next = -1;
            double radius_min = typeMethods.r8_huge();

            int j1;
            for (j1 = 0; j1 < n; j1++)
            {
                if (status[j1] != FREE)
                {
                    continue;
                }

                double radius_i;
                switch (walls)
                {
                    case true:
                    {
                        radius_i = typeMethods.r8_huge();
                        for (i = 0; i < dim_num; i++)
                        {
                            radius_i = Math.Min(radius_i, z[i + j1 * dim_num]);
                        }

                        for (i = 0; i < dim_num; i++)
                        {
                            radius_i = Math.Min(radius_i, 1.0 - z[i + j1 * dim_num]);
                        }

                        break;
                    }
                    default:
                        radius_i = typeMethods.r8_huge();
                        break;
                }

                int j2;
                for (j2 = 0; j2 < n; j2++)
                {
                    if (j2 != j1)
                    {
                        double distance_j = 0.0;
                        for (i = 0; i < dim_num; i++)
                        {
                            distance_j += Math.Pow(z[i + j1 * dim_num] - z[i + j2 * dim_num], 2);
                        }

                        distance_j = Math.Sqrt(distance_j);

                        radius_i = status[j2] == FREE ? Math.Min(radius_i, distance_j / 2.0) : Math.Min(radius_i, distance_j - radius[j2]);
                    }
                }

                if (!(radius_i < radius_min))
                {
                    continue;
                }

                next = j1;
                radius_min = radius_i;

            }

            switch (next)
            {
                case -1:
                    Console.WriteLine("");
                    Console.WriteLine("RADIUS_MAXIMUS - Fatal error!");
                    Console.WriteLine("  There were points left to handle, but could");
                    Console.WriteLine("  not choose the next one to work on.");
                    return null;
                default:
                    radius[next] = radius_min;
                    status[next] = FIXED;
                    break;
            }
        }

        return radius;
    }

    public static double sphere_measure(int dim_num, int n, double[] z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_MEASURE determines the pointset quality measure S.
        //
        //  Discussion:
        //
        //    This routine computes a measure of even spacing for a set of N points
        //    in the DIM_NUM-dimensional unit hypercube.  We will discuss the program
        //    as though the space is 2-dimensional and the spheres are circles, but
        //    the program may be used for general DIM_NUM-dimensional data.
        //
        //    The points are assumed to lie in the unit square.
        //
        //    The program makes a circle-packing measurement on the points
        //    by assuming that, at each point, a circle is centered; all
        //    the circles start out with zero radius, and then expand
        //    together at the same rate.  A circle stops expanding as soon
        //    as it touches any other circle.
        //
        //    The amount of area covered by the circles is compared to the
        //    area of the unit square.  This measurement has a certain amount
        //    of boundary effect: some circles will naturally extend outside
        //    the unit hypercube.  If this is a concern, is possible to restrict
        //    the circles to remain inside the unit hypercube.  In any case,
        //    this problem generally goes away as the number of points increases.
        //
        //    Since we are interested in the coverage of the unit hypercube,
        //    it is probably best if the circles are restricted.  This way,
        //    computing the area of the circles gives a measure of the even
        //    coverage of the region, relative to the presumably best possible
        //    covering, by the same number of circles, but of equal radius.
        //
        //    In the limit, the maximum relative packing density of a 2D
        //    region with equal-sized circles is 0.9069.  In 3D, a density
        //    of at least 0.74 can be achieved, and it is known that no
        //    greater than 0.7796 is possible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double Z[DIM_NUM*N], the points.
        //
        //    Output, double SPHERE_MEASURE, the amount of volume taken up
        //    by the nonintersecting spheres of maximum radius around each
        //    point.  Ignoring boundary effects, the "ideal" value would be
        //    1 (achievable only in 1 dimension), and the maximum value
        //    possible is the sphere packing density in the given spatial
        //    dimension.  If boundary effects can be ignored, the value of
        //    SPHERE_VOLUME reports how closely the given set of points
        //    behaves like a set of close-packed spheres.
        //
        //  Local Parameters:
        //
        //    Local, logical WALLS, is TRUE if the spheres are restricted
        //    to lie within the unit hypercube.
        //
    {
        int i;
        const bool verbose = false;
        const bool walls = true;

        if (!typeMethods.r8mat_in_01(dim_num, n, z))
        {
            Console.WriteLine("");
            Console.WriteLine("SPHERE_MEASURE - Fatal error!");
            Console.WriteLine("  Some of the data is not inside the unit hypercube.");
            return typeMethods.r8_huge();
        }

        double[] radius = radius_maximus(dim_num, n, z, walls);

        double sphere = 0.0;
        for (i = 0; i < n; i++)
        {
            double volume = Sphere.sphere_volume_nd(dim_num, radius[i]);
            sphere += volume;
        }

        switch (verbose)
        {
            case true:
            {
                double radius_ave = 0.0;
                double radius_min = typeMethods.r8_huge();
                double radius_max = 0.0;
                int j;
                for (j = 0; j < n; j++)
                {
                    radius_ave += radius[j];
                    radius_min = Math.Min(radius_min, radius[j]);
                    radius_max = Math.Max(radius_max, radius[j]);
                }

                Console.WriteLine("");
                Console.WriteLine("  Number of dimensions is " + dim_num + "");
                Console.WriteLine("  Number of points is " + n + "");
                switch (walls)
                {
                    case true:
                        Console.WriteLine("  Spheres are required to stay in the unit hypercube.");
                        break;
                    default:
                        Console.WriteLine("  Spheres are NOT required to stay in the unit hypercube.");
                        break;
                }

                Console.WriteLine("");
                Console.WriteLine("  Average radius = " + radius_ave + "");
                Console.WriteLine("  Minimum radius = " + radius_min + "");
                Console.WriteLine("  Maximum radius = " + radius_max + "");
                Console.WriteLine("  Sphere volume =  " + sphere + "");
                break;
            }
        }

        return sphere;
    }

    public static double tau_measure(int dim_num, int n, double[] z, int ns,
            Func<int, int, int, GeometrySampleResult> sample_routine,
            int seed_init)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TAU_MEASURE determines the pointset quality measure TAU.
        //
        //  Discussion:
        //
        //    The TAU measure of point distribution quality for a set Z of
        //    N points in an DIM_NUM-dimensional region is defined as follows:
        //
        //    For each point Z(I) in the pointset, let V(I) be the subregion
        //    defined by the intersection of the region with the Voronoi
        //    region associated with Z(I).
        //
        //    Let T(I) be the trace of the second moment tensor about the point
        //    Z(I), associated with the subregion V(I).  Let T_BAR be the average
        //    of the values of T(1:N).
        //
        //    Then TAU = maximum ( 1 <= I <= N ) abs ( T(I) - TBAR ).
        //
        //    This quantity can be estimated using sampling.  A given number of
        //    sample points are generated in the region, assigned to the nearest
        //    element of the pointset, and used to approximate the Voronoi regions
        //    and the second moment tensors.
        //
        //    In an ideally regular mesh, the values of T would be equal, and so
        //    TAU would be zero.  In general, the smaller TAU, the better.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[DIM_NUM*N], the point distribution.
        //
        //    Input, int NS, the number of sample points.
        //
        //    Input, double *SAMPLE_ROUTINE, the name of a routine which
        //    is used to produce an DIM_NUM by N array of sample points in the region,
        //    of the form:
        //      double *sample_routine ( int dim_num, int n, int *seed )
        //
        //    Input, int SEED_INIT, the initial value of the random number seed.
        //
        //    Output, double TAU_MEASURE, a quality measure.
        //
    {
        int[] closest = new int[1];
        int i;
        int i1;
        int i2;
        int j;
        int k;

        double[] centroid = new double[dim_num * n];
        int[] hit = new int[n];
        double[] moment = new double[dim_num * dim_num * n];
        double[] t = new double[n];

        int seed = seed_init;

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < dim_num; i++)
            {
                centroid[i + j * dim_num] = 0.0;
            }
        }

        for (j = 0; j < n; j++)
        {
            hit[j] = 0;
        }

        for (j = 0; j < n; j++)
        {
            for (i2 = 0; i2 < dim_num; i2++)
            {
                for (i1 = 0; i1 < dim_num; i1++)
                {
                    moment[i1 + i2 * dim_num + j * dim_num * dim_num] = 0.0;
                }
            }
        }

        for (k = 1; k <= ns; k++)
        {
            GeometrySampleResult result = sample_routine(dim_num, 1, seed);
            double[] x = result.result;
            seed = result.seed;

            find_closest(dim_num, n, 1, x, z, closest);

            hit[closest[0]] += 1;

            for (i = 0; i < dim_num; i++)
            {
                centroid[i + closest[0] * dim_num] += x[i];
            }

            for (i1 = 0; i1 < dim_num; i1++)
            {
                for (i2 = 0; i2 < dim_num; i2++)
                {
                    moment[i1 + i2 * dim_num + closest[0] * dim_num * dim_num] += x[i1] * x[i2];
                }
            }
        }

        for (j = 0; j < n; j++)
        {
            switch (hit[j])
            {
                case > 0:
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        centroid[i + j * dim_num] /= hit[j];
                    }

                    for (i1 = 0; i1 < dim_num; i1++)
                    {
                        for (i2 = 0; i2 < dim_num; i2++)
                        {
                            moment[i1 + i2 * dim_num + j * dim_num * dim_num] /= hit[j];
                        }
                    }

                    for (i1 = 0; i1 < dim_num; i1++)
                    {
                        for (i2 = 0; i2 < dim_num; i2++)
                        {
                            moment[i1 + i2 * dim_num + j * dim_num * dim_num] -= centroid[i1 + j * dim_num] * centroid[i2 + j * dim_num];
                        }
                    }

                    break;
                }
            }
        }

        for (j = 0; j < n; j++)
        {
            t[j] = 0.0;
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < dim_num; i++)
            {
                t[j] += moment[i + i * dim_num + j * dim_num * dim_num];
            }
        }

        double t_bar = 0.0;

        for (j = 0; j < n; j++)
        {
            t_bar += t[j];
        }

        t_bar /= n;

        double tau = 0.0;
        for (j = 0; j < n; j++)
        {
            tau = Math.Max(tau, Math.Abs(t[j] - t_bar));
        }

        return tau;
    }

}