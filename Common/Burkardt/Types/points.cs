using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static int point_unique_count ( int m, int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_UNIQUE_COUNT counts the unique points.
        //
        //  Discussion:
        //
        //    The input data is an M x N array A, representing the M-dimensional
        //    coordinates of N points.
        //
        //    The algorithm relies on the fact that, in a sorted list, points that
        //    are exactly equal must occur consecutively.
        //
        //    The output is the number of unique points in the list.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, double A[M*N], the array of N columns of data.
        //
        //    Output, int POINT_UNIQUE_COUNT, the number of unique points.
        //
    {
        int j;
        int unique_num;

        switch (n)
        {
            case <= 0:
                unique_num = 0;
                return unique_num;
        }
        //
        //  Implicitly sort the array.
        //
        int[] indx = r8col_sort_heap_index_a ( m, n, a );
        //
        //  Two points are considered equal only if they exactly match.
        //  In that case, equal points can only occur as consecutive items
        //  in the sorted list.   This makes counting easy.
        //
        unique_num = 1;
        int unique_index = indx[0];

        for ( j = 1; j < n; j++ )
        {
            int i;
            for ( i = 0; i < m; i++ )
            {
                if ( Math.Abs(a[i+unique_index*m] - a[i+indx[j]*m]) > double.Epsilon )
                {
                    unique_num += 1;
                    unique_index = indx[j];
                }
            }
        }
            
        return unique_num;
    }
        
    public static void points_plot(string file_name, int node_num, double[] node_xy,
            bool node_label)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_PLOT plots a pointset.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string FILE_NAME, the name of the file to create.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the nodes.
        //
        //    Input, bool NODE_LABEL, is TRUE if the nodes are to be labeled.
        //
        //  Local parameters:
        //
        //    int CIRCLE_SIZE, controls the size of the circles depicting
        //    the nodes.  Currently set to 5.  3 is pretty small, and 1 is
        //    barely visible.
        //
    {
        int circle_size = 3;
        int delta;
        List<string> file_unit = new();
        int node;
        int x_ps;
        int x_ps_max = 576;
        int x_ps_max_clip = 594;
        int x_ps_min = 36;
        int x_ps_min_clip = 18;
        int y_ps;
        int y_ps_max = 666;
        int y_ps_max_clip = 684;
        int y_ps_min = 126;
        int y_ps_min_clip = 108;

        //
        //  We need to do some figuring here, so that we can determine
        //  the range of the data, and hence the height and width
        //  of the piece of paper.
        //
        double x_max = -r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (x_max < node_xy[0 + node * 2])
            {
                x_max = node_xy[0 + node * 2];
            }
        }

        double x_min = r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (node_xy[0 + node * 2] < x_min)
            {
                x_min = node_xy[0 + node * 2];
            }
        }

        double x_scale = x_max - x_min;

        x_max += 0.05 * x_scale;
        x_min -= 0.05 * x_scale;
        x_scale = x_max - x_min;

        double y_max = -r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (y_max < node_xy[1 + node * 2])
            {
                y_max = node_xy[1 + node * 2];
            }
        }

        double y_min = r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (node_xy[1 + node * 2] < y_min)
            {
                y_min = node_xy[1 + node * 2];
            }
        }

        double y_scale = y_max - y_min;

        y_max += 0.05 * y_scale;
        y_min -= 0.05 * y_scale;
        y_scale = y_max - y_min;

        if (x_scale < y_scale)
        {
            delta = (int) ((x_ps_max - x_ps_min)
                * (y_scale - x_scale) / (2.0 * y_scale));

            x_ps_max -= delta;
            x_ps_min += delta;

            x_ps_max_clip -= delta;
            x_ps_min_clip += delta;

            x_scale = y_scale;
        }
        else if (y_scale < x_scale)
        {
            delta = (int) ((y_ps_max - y_ps_min)
                * (x_scale - y_scale) / (2.0 * x_scale));

            y_ps_max -= delta;
            y_ps_min += delta;

            y_ps_max_clip -= delta;
            y_ps_min_clip += delta;

            y_scale = x_scale;
        }

        file_unit.Add("%!PS-Adobe-3.0 EPSF-3.0");
        file_unit.Add("%%Creator: points_plot.C");
        file_unit.Add("%%Title: " + file_name + "");

        file_unit.Add("%%Pages: 1");
        file_unit.Add("%%BoundingBox:  "
                      + x_ps_min + "  "
                      + y_ps_min + "  "
                      + x_ps_max + "  "
                      + y_ps_max + "");
        file_unit.Add("%%Document-Fonts: Times-Roman");
        file_unit.Add("%%LanguageLevel: 1");
        file_unit.Add("%%EndComments");
        file_unit.Add("%%BeginProlog");
        file_unit.Add("/inch {72 mul} def");
        file_unit.Add("%%EndProlog");
        file_unit.Add("%%Page:      1     1");
        file_unit.Add("save");
        file_unit.Add("%");
        file_unit.Add("% Set the RGB line color to very light gray.");
        file_unit.Add("%");
        file_unit.Add(" 0.9000 0.9000 0.9000 setrgbcolor");
        file_unit.Add("%");
        file_unit.Add("% Draw a gray border around the page.");
        file_unit.Add("%");
        file_unit.Add("newpath");
        file_unit.Add(x_ps_min + "  "
                               + y_ps_min + "  moveto");
        file_unit.Add(x_ps_max + "  "
                               + y_ps_min + "  lineto");
        file_unit.Add(x_ps_max + "  "
                               + y_ps_max + "  lineto");
        file_unit.Add(x_ps_min + "  "
                               + y_ps_max + "  lineto");
        file_unit.Add(x_ps_min + "  "
                               + y_ps_min + "  lineto");
        file_unit.Add("stroke");
        file_unit.Add("%");
        file_unit.Add("% Set RGB line color to black.");
        file_unit.Add("%");
        file_unit.Add(" 0.0000 0.0000 0.0000 setrgbcolor");
        file_unit.Add("%");
        file_unit.Add("%  Set the font and its size:");
        file_unit.Add("%");
        file_unit.Add("/Times-Roman findfont");
        file_unit.Add("0.50 inch scalefont");
        file_unit.Add("setfont");
        file_unit.Add("%");
        file_unit.Add("%  Print a title:");
        file_unit.Add("%");
        file_unit.Add("%  210  702 moveto");
        file_unit.Add("%(Pointset) show");
        file_unit.Add("%");
        file_unit.Add("% Define a clipping polygon");
        file_unit.Add("%");
        file_unit.Add("newpath");
        file_unit.Add(x_ps_min_clip + "  "
                                    + y_ps_min_clip + "  moveto");
        file_unit.Add(x_ps_max_clip + "  "
                                    + y_ps_min_clip + "  lineto");
        file_unit.Add(x_ps_max_clip + "  "
                                    + y_ps_max_clip + "  lineto");
        file_unit.Add(x_ps_min_clip + "  "
                                    + y_ps_max_clip + "  lineto");
        file_unit.Add(x_ps_min_clip + "  "
                                    + y_ps_min_clip + "  lineto");
        file_unit.Add("clip newpath");
        //
        //  Draw the nodes.
        //
        file_unit.Add("%");
        file_unit.Add("%  Draw filled dots at each node:");
        file_unit.Add("%");
        file_unit.Add("%  Set the color to blue:");
        file_unit.Add("%");
        file_unit.Add("0.000  0.150  0.750  setrgbcolor");
        file_unit.Add("%");

        for (node = 0; node < node_num; node++)
        {
            x_ps = (int) (
                ((x_max - node_xy[0 + node * 2]) * x_ps_min
                 + (+node_xy[0 + node * 2] - x_min) * x_ps_max)
                / (x_max - x_min));

            y_ps = (int) (
                ((y_max - node_xy[1 + node * 2]) * y_ps_min
                 + (node_xy[1 + node * 2] - y_min) * y_ps_max)
                / (y_max - y_min));

            file_unit.Add("newpath  "
                          + x_ps + "  "
                          + y_ps + "  "
                          + circle_size + " 0 360 arc closepath fill");
        }

        //
        //  Label the nodes.
        //
        file_unit.Add("%");
        file_unit.Add("%  Label the nodes:");
        file_unit.Add("%");
        file_unit.Add("%  Set the color to darker blue:");
        file_unit.Add("%");
        file_unit.Add("0.000  0.250  0.850  setrgbcolor");
        file_unit.Add("/Times-Roman findfont");
        file_unit.Add("0.20 inch scalefont");
        file_unit.Add("setfont");

        file_unit.Add("%");

        for (node = 0; node < node_num; node++)
        {
            x_ps = (int) (
                ((x_max - node_xy[0 + node * 2]) * x_ps_min
                 + (+node_xy[0 + node * 2] - x_min) * x_ps_max)
                / (x_max - x_min));

            y_ps = (int) (
                ((y_max - node_xy[1 + node * 2]) * y_ps_min
                 + (node_xy[1 + node * 2] - y_min) * y_ps_max)
                / (y_max - y_min));

            file_unit.Add("newpath  "
                          + x_ps + "  "
                          + (y_ps + 5) + "  moveto ("
                          + (node + 1) + ") show");
        }

        file_unit.Add("%");
        file_unit.Add("restore showpage");
        file_unit.Add("%");
        file_unit.Add("% End of page");
        file_unit.Add("%");
        file_unit.Add("%%Trailer");
        file_unit.Add("%%EOF");

        try
        {
            File.WriteAllLines(file_name, file_unit);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("POINTS_PLOT - Fatal error!");
            Console.WriteLine("  Could not open the output EPS file.");
        }
    }

    public static int point_radial_tol_unique_count(int m, int n, double[] a, double tol,
            ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_RADIAL_TOL_UNIQUE_COUNT counts the tolerably unique points.
        //
        //  Discussion:
        //
        //    The input data is an M x N array A, representing the M-dimensional
        //    coordinates of N points.
        //
        //    The output is the number of tolerably unique points in the list.
        //
        //    This program performs the same task as POINT_TOL_UNIQUE_COUNT.
        //    But that program is guaranteed to use N^2 comparisons.
        //
        //    It is hoped that this function, on the other hand, will tend
        //    to use O(N) comparisons after an O(NLog(N)) sort.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, double A[M*N], the array of N columns of data.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Input/output, int *SEED, a seed for the random
        //    number generator.
        //
        //    Output, int POINT_RADIAL_TOL_UNIQUE_COUNT, the number of tolerably
        //    unique points.
        //
    {
        int i;
        int j;
        int unique_num;

        switch (n)
        {
            case <= 0:
                unique_num = 0;
                return unique_num;
        }

        //
        //  Assign a base point Z randomly in the convex hull.
        //
        double[] w = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        double w_sum = r8vec_sum(n, w);
        for (j = 0; j < n; j++)
        {
            w[j] /= w_sum;
        }

        double[] z = new double[m];
        for (i = 0; i < m; i++)
        {
            z[i] = 0.0;
            for (j = 0; j < n; j++)
            {
                z[i] += a[i + j * m] * w[j];
            }
        }

        //
        //  Compute the radial distance R of each point to Z.
        //
        double[] r = new double[n];

        for (j = 0; j < n; j++)
        {
            r[j] = 0.0;
            for (i = 0; i < m; i++)
            {
                r[j] += Math.Pow(a[i + j * m] - z[i], 2);
            }

            r[j] = Math.Sqrt(r[j]);
        }

        //
        //  Implicitly sort the R array.
        //
        int[] indx = r8vec_sort_heap_index_a_new(n, r);
        //
        //  To determine if a point I is tolerably unique, we only have to check
        //  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
        //
        unique_num = 0;

        bool[] unique = new bool[n];
        for (i = 0; i < n; i++)
        {
            unique[i] = true;
        }

        for (i = 0; i < n; i++)
        {
            switch (unique[indx[i]])
            {
                case true:
                {
                    //
                    //  Point INDX(I) is unique, in that no earlier point is near it.
                    //
                    unique_num += 1;
                    //
                    //  Look for later points which are close to point INDX(I)
                    //  in terms of R.
                    //
                    int hi = i;

                    while (hi < n - 1)
                    {
                        if (r[indx[i]] + tol < r[indx[hi + 1]])
                        {
                            break;
                        }

                        hi += 1;
                    }

                    //
                    //  Points INDX(I+1) through INDX(HI) have an R value close to
                    //  point INDX(I).  Are they truly close to point INDEX(I)?
                    //
                    for (j = i + 1; j <= hi; j++)
                    {
                        switch (unique[indx[j]])
                        {
                            case true:
                            {
                                double dist = 0.0;
                                int k;
                                for (k = 0; k < m; k++)
                                {
                                    dist += Math.Pow(a[k + indx[i] * m] - a[k + indx[j] * m], 2);
                                }

                                dist = Math.Sqrt(dist);

                                if (dist <= tol)
                                {
                                    unique[indx[j]] = false;
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }

        return unique_num;
    }

    public static void point_radial_tol_unique_count_inc1(int m, int n1, double[] a1,
            double tol, ref int seed, ref double[] z, ref double[] r1, ref int[] indx1, ref bool[] unique1,
            ref int unique_num1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_RADIAL_TOL_UNIQUE_COUNT_INC1 counts the tolerably unique points.
        //
        //  Discussion:
        //
        //    The input data includes an M x N1 array A1 of a set of N1
        //    "permanent" points and N2 "temporary" points.
        //
        //    This is a two step version of POINT_RADIAL_TOL_UNIQUE_COUNT_INC.
        //
        //    This means that we want to identify the tolerably unique points
        //    among the permanent points before processing the temporary points.
        //
        //    If many sets of temporary data are considered, this function will
        //    do a lot of unnecessary work resorting the permanent data; it would
        //    be possible to avoid repetitions of that work at the expense of saving
        //    various work vectors.  This function accepts the overhead of the
        //    repeated calculations for the benefit of only having to "remember"
        //    the number of unique points discovered.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N1, the number of permanent points.
        //
        //    Input, double A1[M*N1], the permanent points.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Input/output, int *SEED, a seed for the random
        //    number generator.
        //
        //    Output, double Z[M], a random base vector used to
        //    linearly sort the data.
        //
        //    Output, double R1[N1], the scalar values assigned to
        //    the data for sorting.
        //
        //    Output, int INDX1[N1], the ascending sort index
        //    for A1.
        //
        //    Output, bool UNIQUE1[N1], is TRUE for each unique permanent point.
        //
        //    Output, int *UNIQUE_NUM1, the number of tolerably
        //    unique permanent points.
        //
    {
        int i;
        int j1;
        //
        //  Assign a base point Z randomly in the convex hull of the permanent points.
        //
        double[] w = UniformRNG.r8vec_uniform_01_new(n1, ref seed);
        double w_sum = r8vec_sum(n1, w);
        for (j1 = 0; j1 < n1; j1++)
        {
            w[j1] /= w_sum;
        }

        for (i = 0; i < m; i++)
        {
            z[i] = 0.0;
            for (j1 = 0; j1 < n1; j1++)
            {
                z[i] += a1[i + j1 * m] * w[j1];
            }
        }

        //
        //  Initialize the permanent point data.
        //
        for (j1 = 0; j1 < n1; j1++)
        {
            r1[j1] = 0.0;
            for (i = 0; i < m; i++)
            {
                r1[j1] += Math.Pow(a1[i + j1 * m] - z[i], 2);
            }

            r1[j1] = Math.Sqrt(r1[j1]);
        }

        indx1 = r8vec_sort_heap_index_a(n1, r1);

        unique_num1 = 0;
        for (j1 = 0; j1 < n1; j1++)
        {
            unique1[j1] = true;
        }

        //
        //  STEP 1:
        //  Compare PERMANENT POINTS to PERMANENT POINTS.
        //
        for (j1 = 0; j1 < n1; j1++)
        {
            switch (unique1[indx1[j1]])
            {
                case true:
                {
                    unique_num1 += 1;

                    int hi = j1;

                    while (hi < n1 - 1)
                    {
                        if (r1[indx1[j1]] + tol < r1[indx1[hi + 1]])
                        {
                            break;
                        }

                        hi += 1;
                    }

                    int k1;
                    for (k1 = j1 + 1; k1 <= hi; k1++)
                    {
                        switch (unique1[indx1[k1]])
                        {
                            case true:
                            {
                                double dist = 0.0;
                                for (i = 0; i < m; i++)
                                {
                                    dist += Math.Pow(a1[i + indx1[j1] * m] - a1[i + indx1[k1] * m], 2);
                                }

                                dist = Math.Sqrt(dist);

                                if (dist <= tol)
                                {
                                    unique1[indx1[k1]] = false;
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }
    }

    public static void point_radial_tol_unique_count_inc2(int m, int n1, double[] a1, int n2,
            double[] a2, double tol, double[] z, double[] r1, int[] indx1, bool[] unique1,
            ref int unique_num2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_RADIAL_TOL_UNIQUE_COUNT_INC2 counts the tolerably unique points.
        //
        //  Discussion:
        //
        //    The input data includes an M x N1 array A1 and an M x N2 array A2,
        //    representing the M-dimensional coordinates of a set of N1
        //    "permanent" points and N2 "temporary" points.
        //
        //    This is an "incremental" version of POINT_RADIAL_TOL_UNIQUE_COUNT.
        //
        //    This means that we want to identify the tolerably unique points
        //    among the permanent points before processing the temporary points.
        //
        //    If many sets of temporary data are considered, this function will
        //    do a lot of unnecessary work resorting the permanent data; it would
        //    be possible to avoid repetitions of that work at the expense of saving
        //    various work vectors.  This function accepts the overhead of the
        //    repeated calculations for the benefit of only having to "remember"
        //    the number of unique points discovered.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N1, the number of permanent points.
        //
        //    Input, double A1[M*N1], the permanent points.
        //
        //    Input, int N2, the number of temporary points.
        //
        //    Input, double A2[M*N2], the temporary points.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Input, double Z[M], a random base vector used to
        //    linearly sort the data.
        //
        //    Input, double R1[N1], the scalar values assigned to
        //    the data for sorting.
        //
        //    Input, int INDX1[N1], the ascending sort index
        //    for A1.
        //
        //    Input, bool UNIQUE1[N1], is TRUE for each unique permanent point.
        //
        //    Output, int *UNIQUE_NUM2, the number of additional
        //    tolerably unique points if the temporary points are included.
        //
    {
        double dist;
        int i;
        int j1;
        int j2;
        int j2_hi = 0;
        int j2_lo = 0;
        //
        //  Initialize the temporary point data.
        //
        double[] r2 = new double[n2];
        for (j2 = 0; j2 < n2; j2++)
        {
            r2[j2] = 0.0;
            for (i = 0; i < m; i++)
            {
                r2[j2] += Math.Pow(a2[i + j2 * m] - z[i], 2);
            }

            r2[j2] = Math.Sqrt(r2[j2]);
        }

        int[] indx2 = r8vec_sort_heap_index_a(n2, r2);

        bool[] unique2 = new bool[n2];
        for (j2 = 0; j2 < n2; j2++)
        {
            unique2[j2] = true;
        }

        unique_num2 = 0;
        //
        //  STEP 2:
        //  Use PERMANENT points to eliminate TEMPORARY points.
        //
        for (j1 = 0; j1 < n1; j1++)
        {
            switch (unique1[indx1[j1]])
            {
                case true:
                {
                    double r_lo = r1[indx1[j1]] - tol;
                    double r_hi = r1[indx1[j1]] + tol;

                    r8vec_index_sorted_range(n2, r2, indx2, r_lo, r_hi,
                        ref j2_lo, ref j2_hi);

                    for (j2 = j2_lo; j2 <= j2_hi; j2++)
                    {
                        switch (unique2[indx2[j2]])
                        {
                            case true:
                            {
                                dist = 0.0;
                                for (i = 0; i < m; i++)
                                {
                                    dist += Math.Pow(a1[i + indx1[j1] * m]
                                                     - a2[i + indx2[j2] * m], 2);
                                }

                                dist = Math.Sqrt(dist);
                                if (dist <= tol)
                                {
                                    unique2[indx2[j2]] = false;
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }

        //
        //  STEP 3:
        //  Use TEMPORARY points to eliminate TEMPORARY points.
        //
        for (j2 = 0; j2 < n2; j2++)
        {
            switch (unique2[indx2[j2]])
            {
                case true:
                {
                    unique_num2 += 1;

                    int hi = j2;

                    while (hi < n2 - 1)
                    {
                        if (r2[indx2[j2]] + tol < r2[indx2[hi + 1]])
                        {
                            break;
                        }

                        hi += 1;
                    }

                    int k2;
                    for (k2 = j2 + 1; k2 <= hi; k2++)
                    {
                        switch (unique2[indx2[k2]])
                        {
                            case true:
                            {
                                dist = 0.0;
                                for (i = 0; i < m; i++)
                                {
                                    dist += Math.Pow(a2[i + indx2[j2] * m] - a2[i + indx2[k2] * m], 2);
                                }

                                dist = Math.Sqrt(dist);

                                if (dist <= tol)
                                {
                                    unique2[indx2[k2]] = false;
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }
    }

    public static int point_radial_unique_count(int m, int n, double[] a, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_RADIAL_UNIQUE_COUNT counts the unique points.
        //
        //  Discussion:
        //
        //    The input data is an M x N array A, representing the M-dimensional
        //    coordinates of N points.
        //
        //    The output is the number of unique points in the list.
        //
        //    This program performs the same task as POINT_UNIQUE_COUNT, and
        //    carries out more work.  Hence, it is not a substitute for
        //    POINT_UNIQUE_COUNT.  Instead, it is intended to be a starting point
        //    for a similar program which includes a tolerance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, double A[M*N], the array of N columns of data.
        //
        //    Input/output, int *SEED, a seed for the random
        //    number generator.
        //
        //    Output, int POINT_RADIAL_UNIQUE_COUNT, the number of unique points.
        //
    {
        bool equal = false;
        int i;
        int j;
        int unique_num;

        switch (n)
        {
            case <= 0:
                unique_num = 0;
                return unique_num;
        }

        //
        //  Assign a base point Z randomly in the convex hull.
        //
        double[] w = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        double w_sum = r8vec_sum(n, w);
        for (j = 0; j < n; j++)
        {
            w[j] /= w_sum;
        }

        double[] z = new double[m];
        for (i = 0; i < m; i++)
        {
            z[i] = 0.0;
            for (j = 0; j < n; j++)
            {
                z[i] += a[i + j * m] * w[j];
            }
        }

        //
        //  Compute the radial distance R of each point to Z.
        //
        double[] r = new double[n];

        for (j = 0; j < n; j++)
        {
            r[j] = 0.0;
            for (i = 0; i < m; i++)
            {
                r[j] += Math.Pow(a[i + j * m] - z[i], 2);
            }

            r[j] = Math.Sqrt(r[j]);
        }

        //
        //  Implicitly sort the R array.
        //
        int[] indx = r8vec_sort_heap_index_a(n, r);
        //
        //  To determine if a point is unique, we only have to check
        //  whether it is distinct from all points with the same
        //  R value and lower ordering.
        //
        unique_num = 0;
        int hi = -1;

        while (hi < n - 1)
        {
            //
            //  Advance LO.
            //
            int lo = hi + 1;
            //
            //  Extend HI.
            //
            hi = lo;

            while (hi < n - 1)
            {
                if (Math.Abs(r[indx[hi + 1]] - r[indx[lo]]) <= double.Epsilon)
                {
                    hi += 1;
                }
                else
                {
                    break;
                }
            }

            //
            //  Points INDX(LO) through INDX(HI) have same R value.
            //
            //  Find the unique ones.
            //
            unique_num += 1;

            int j1;
            for (j1 = lo + 1; j1 <= hi; j1++)
            {
                int j2;
                for (j2 = lo; j2 < j1; j2++)
                {
                    equal = true;
                    for (i = 0; i < m; i++)
                    {
                        if (!(Math.Abs(a[i + indx[j2] * m] - a[i + indx[j1] * m]) > double.Epsilon))
                        {
                            continue;
                        }

                        equal = false;
                        break;
                    }

                    if (equal)
                    {
                        break;
                    }
                }

                switch (equal)
                {
                    case false:
                        unique_num += 1;
                        break;
                }
            }
        }

        return unique_num;
    }

    public static int point_tol_unique_count(int m, int n, double[] a, double tol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_TOL_UNIQUE_COUNT counts the tolerably unique points.
        //
        //  Discussion:
        //
        //    The input data is an M x N array A, representing the M-dimensional
        //    coordinates of N points.
        //
        //    This function uses a simple but expensive approach.  The first point
        //    is accepted as unique.  Each subsequent point is accepted as unique
        //    only if it is at least a tolerance away from all accepted unique points.
        //    This means the expected amount of work is O(N^2).
        //
        //    The output is the number of unique points in the list.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, double A[M*N], the array of N columns of data.
        //
        //    Input, double TOL, a tolerance.
        //
        //    Output, int POINT_TOL_UNIQUE_COUNT, the number of unique points.
        //
    {
        int i;

        bool[] unique = new bool[n];

        for (i = 0; i < n; i++)
        {
            unique[i] = true;
        }

        int unique_num = n;

        for (i = 1; i < n; i++)
        {
            int j;
            for (j = 0; j < i; j++)
            {
                if (!unique[j])
                {
                    continue;
                }

                double dist = 0.0;
                int k;
                for (k = 0; k < m; k++)
                {
                    dist += Math.Pow(a[k + i * m] - a[k + j * m], 2);
                }

                dist = Math.Sqrt(dist);
                if (!(dist <= tol))
                {
                    continue;
                }

                unique[i] = false;
                unique_num -= 1;
                break;
            }
        }

        return unique_num;
    }

    public static int point_tol_unique_index ( int m, int n, double[] a, double tol, ref int[] xdnu )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_TOL_UNIQUE_INDEX indexes the tolerably unique points.
        //
        //  Discussion:
        //
        //    This routine uses an algorithm that is O(N^2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the data values.
        //
        //    Input, int N, the number of data values.
        //
        //    Input, double A[M*N], the data values.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int XDNU[N], the index, in A, of the tolerably unique
        //    point that "represents" this point.
        //
        //    Output, int POINT_TOL_UNIQUE_INDEX, the number of tolerably
        //    unique points.
        //
    {
        int i;

        bool[] unique = new bool[n];

        for ( i = 0; i < n; i++ )
        {
            unique[i] = true;
        }
        for ( i = 0; i < n; i++ )
        {
            xdnu[i] = i;
        }
        int unique_num = n;

        i = 0;
        xdnu[0] = 0;

        for ( i = 1; i < n; i++ )
        {
            int j;
            for ( j = 0; j < i; j++ )
            {
                if (!unique[j])
                {
                    continue;
                }

                double dist = 0.0;
                int k;
                for ( k = 0; k < m; k++ )
                {
                    dist += Math.Pow ( a[k+i*m] - a[k+j*m], 2 );
                }
                dist = Math.Sqrt ( dist );
                if (!(dist <= tol))
                {
                    continue;
                }

                unique[i] = false;
                unique_num -= 1;
                xdnu[i] = j;
                break;
            }
        }
            
        return unique_num;
    }
        
    public static int point_radial_tol_unique_index(int m, int n, double[] a, double tol,
            ref int seed, ref int[] undx, ref int[] xdnu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_RADIAL_TOL_UNIQUE_INDEX indexes the tolerably unique points.
        //
        //  Discussion:
        //
        //    The input data is an M x N array A, representing the M-dimensional
        //    coordinates of N points.
        //
        //    The output is:
        //    * the number of tolerably unique points in the list;
        //    * the index, in the list of unique items, of the representatives
        //      of each point;
        //    * the index, in A, of the tolerably unique representatives.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, double A[M*N], the array of N columns of data.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Input/output, int SEED, a seed for the random
        //    number generator.
        //
        //    Output, int UNDX[UNIQUE_NUM], the index, in A, of the
        //    tolerably unique points.
        //
        //    Output, int XDNU[N], the index, in UNDX, of the
        //    tolerably unique point that "represents" this point.
        //
        //    Output, int POINT_RADIAL_TOL_UNIQUE_INDEX, the number of tolerably
        //    unique points.
        //
    {
        int i;
        int j;
        int unique_num;

        switch (n)
        {
            case <= 0:
                unique_num = 0;
                return unique_num;
        }

        //
        //  Assign a base point Z randomly in the convex hull.
        //
        double[] w = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        double w_sum = r8vec_sum(n, w);
        for (j = 0; j < n; j++)
        {
            w[j] /= w_sum;
        }

        double[] z = new double[m];
        for (i = 0; i < m; i++)
        {
            z[i] = 0.0;
            for (j = 0; j < n; j++)
            {
                z[i] += a[i + j * m] * w[j];
            }
        }

        //
        //  Compute the radial distance R of each point to Z.
        //
        double[] r = new double[n];

        for (j = 0; j < n; j++)
        {
            r[j] = 0.0;
            for (i = 0; i < m; i++)
            {
                r[j] += Math.Pow(a[i + j * m] - z[i], 2);
            }

            r[j] = Math.Sqrt(r[j]);
        }

        //
        //  Implicitly sort the R array.
        //
        int[] indx = r8vec_sort_heap_index_a_new(n, r);
        //
        //  To determine if a point I is tolerably unique, we only have to check
        //  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
        //
        unique_num = 0;

        bool[] unique = new bool[n];
        for (i = 0; i < n; i++)
        {
            unique[i] = true;
        }

        for (i = 0; i < n; i++)
        {
            switch (unique[indx[i]])
            {
                case true:
                {
                    //
                    //  Point INDX(I) is unique, in that no earlier point is near it.
                    //
                    xdnu[indx[i]] = unique_num;
                    undx[unique_num] = indx[i];
                    unique_num += 1;
                    //
                    //  Look for later points which are close to point INDX(I)
                    //  in terms of R.
                    //
                    int hi = i;

                    while (hi < n - 1)
                    {
                        if (r[indx[i]] + tol < r[indx[hi + 1]])
                        {
                            break;
                        }

                        hi += 1;
                    }

                    //
                    //  Points INDX(I+1) through INDX(HI) have an R value close to
                    //  point INDX(I).  Are they truly close to point INDEX(I)?
                    //
                    for (j = i + 1; j <= hi; j++)
                    {
                        switch (unique[indx[j]])
                        {
                            case true:
                            {
                                double dist = 0.0;
                                int k;
                                for (k = 0; k < m; k++)
                                {
                                    dist += Math.Pow(a[k + indx[i] * m] - a[k + indx[j] * m], 2);
                                }

                                dist = Math.Sqrt(dist);

                                if (dist <= tol)
                                {
                                    unique[indx[j]] = false;
                                    xdnu[indx[j]] = xdnu[indx[i]];
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }

        return unique_num;
    }

    public static void point_radial_tol_unique_index_inc1(int m, int n1, double[] a1,
            double tol, ref int seed, ref double[] z, ref double[] r1, ref int[] indx1, ref bool[] unique1,
            ref int unique_num1, ref int[] undx1, int[] xdnu1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_RADIAL_TOL_UNIQUE_INDEX_INC1 indexes the tolerably unique points.
        //
        //  Discussion:
        //
        //    The input data includes an M x N1 array A1 of
        //    "permanent" points.
        //
        //    This is a two step version of POINT_RADIAL_TOL_UNIQUE_INDEX_INC.
        //
        //    The output is:
        //    * the number of tolerably unique points in the list;
        //    * the index, in the list of unique items, of the representatives
        //      of each point;
        //    * the index, in A1, of the tolerably unique representatives.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N1, the number of permanent points.
        //
        //    Input, double A1[M*N1], the permanent points.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Input/output, int *SEED, a seed for the random
        //    number generator.
        //
        //    Output, double Z[M], a random base vector used to
        //    linearly sort the data.
        //
        //    Output, double R1[N1], the scalar values assigned to
        //    the data for sorting.
        //
        //    Output, int INDX1[N1], the ascending sort index for A1.
        //
        //    Output, bool UNIQUE1[N1], is TRUE for unique permanent points.
        //
        //    Output, int *UNIQUE_NUM1, the number of tolerably unique points
        //    with just the permanent points.
        //
        //    Output, int UNDX1[UNIQUE_NUM1], the index, in A1, of the tolerably
        //    unique points.
        //
        //    Output, int XDNU1[N1], the index, in UNDX1, of the tolerably unique
        //    point that "represents" this point.
        //
    {
        int i;
        int j1;
        //
        //  Assign a base point Z randomly in the convex hull of the permanent points.
        //
        double[] w = UniformRNG.r8vec_uniform_01_new(n1, ref seed);
        double w_sum = r8vec_sum(n1, w);
        for (j1 = 0; j1 < n1; j1++)
        {
            w[j1] /= w_sum;
        }

        for (i = 0; i < m; i++)
        {
            z[i] = 0.0;
            for (j1 = 0; j1 < n1; j1++)
            {
                z[i] += a1[i + j1 * m] * w[j1];
            }
        }

        //
        //  Initialize the permanent point data.
        //
        for (j1 = 0; j1 < n1; j1++)
        {
            r1[j1] = 0.0;
            for (i = 0; i < m; i++)
            {
                r1[j1] += Math.Pow(a1[i + j1 * m] - z[i], 2);
            }

            r1[j1] = Math.Sqrt(r1[j1]);
        }

        indx1 = r8vec_sort_heap_index_a(n1, r1);

        unique_num1 = 0;
        for (j1 = 0; j1 < n1; j1++)
        {
            unique1[j1] = true;
        }

        //
        //  STEP 1:
        //  Compare PERMANENT POINTS to PERMANENT POINTS.
        //
        for (j1 = 0; j1 < n1; j1++)
        {
            switch (unique1[indx1[j1]])
            {
                case true:
                {
                    xdnu1[indx1[j1]] = unique_num1;
                    undx1[unique_num1] = indx1[j1];
                    unique_num1 += 1;

                    int hi = j1;

                    while (hi < n1 - 1)
                    {
                        if (r1[indx1[j1]] + tol < r1[indx1[hi + 1]])
                        {
                            break;
                        }

                        hi += 1;
                    }

                    int k1;
                    for (k1 = j1 + 1; k1 <= hi; k1++)
                    {
                        switch (unique1[indx1[k1]])
                        {
                            case true:
                            {
                                double dist = 0.0;
                                for (i = 0; i < m; i++)
                                {
                                    dist += Math.Pow(a1[i + indx1[j1] * m] - a1[i + indx1[k1] * m], 2);
                                }

                                dist = Math.Sqrt(dist);

                                if (dist <= tol)
                                {
                                    unique1[indx1[k1]] = false;
                                    xdnu1[indx1[k1]] = xdnu1[indx1[j1]];
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }
    }

    public static void point_radial_tol_unique_index_inc2(int m, int n1, double[] a1, int n2,
            double[] a2, double tol, double[] z, double[] r1, int[] indx1, bool[] unique1,
            int unique_num1, int[] undx1, int[] xdnu1, ref double[] r2,
            ref int[] indx2, ref bool[] unique2, ref int unique_num2, ref int[] undx2, ref int[] xdnu2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_RADIAL_TOL_UNIQUE_INDEX_INC2 indexes unique temporary points.
        //
        //  Discussion:
        //
        //    The input data includes an M x N1 array A1 and an M x N2 array A2,
        //    representing the M-dimensional coordinates of a set of N1
        //    "permanent" points and N2 "temporary" points.
        //
        //    For notation, we use "A" to describe the M x (N1+N2) array that would be
        //    formed by starting with A1 and appending A2.
        //
        //    The output is:
        //    * the number of tolerably unique points in the list;
        //    * the index, in the list of unique items, of the representatives
        //      of each point;
        //    * the index, in A, of the tolerably unique representatives.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N1, the number of permanent points.
        //
        //    Input, double A1[M*N1], the permanent points.
        //
        //    Input, int N2, the number of temporary points.
        //
        //    Input, double A2[M*N2], the temporary points.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Input, double Z[M], a random base vector used to
        //    linearly sort the data.
        //
        //    Input, double R1[N1], the scalar values assigned to
        //    A1 for sorting.
        //
        //    Input, int INDX1[N1], the ascending sort index for A1.
        //
        //    Input, bool UNIQUE1[N1], is TRUE for unique permanent points.
        //
        //    Input, int UNIQUE_NUM1, the number of tolerably unique permanent points.
        //
        //    Input, int UNDX1[UNIQUE_NUM1],
        //    the index in A1 of the tolerably unique permanent points.
        //
        //    Input, int XDNU1[N1], the index in UNDX1
        //    of the tolerably unique permanent point that "represents" this point.
        //
        //    Output, double R2[N2], the scalar values assigned to
        //    A2 for sorting.
        //
        //    Output, int INDX2[N2], the ascending sort index for A2.
        //
        //    Output, bool UNIQUE2[N2], is TRUE for unique temporary points.
        //
        //    Output, int *UNIQUE_NUM2, the number
        //    of tolerably unique temporary points.
        //
        //    Output, int UNDX2[UNIQUE_NUM2],
        //    the index in A2 of the tolerably unique points, incremented by N1.
        //
        //    Output, int XDNU2[N2], the index, in UNDX1
        //    or UNDX2, of the tolerably unique point that "represents" this
        //    temporary point.  If the value represents an index in UNDX2, this
        //    can be inferred by the fact that its value is greater than or
        //    equal to UNIQUE_NUM1.  To reference UNDX2, the value should then be
        //    decremented by UNIQUE_NUM1.
        //
    {
        double dist;
        int i;
        int j1;
        int j2;
        int j2_hi = 0;
        int j2_lo = 0;
        //
        //  Initialize the temporary point data.
        //
        for (j2 = 0; j2 < n2; j2++)
        {
            r2[j2] = 0.0;
            for (i = 0; i < m; i++)
            {
                r2[j2] += Math.Pow(a2[i + j2 * m] - z[i], 2);
            }

            r2[j2] = Math.Sqrt(r2[j2]);
        }

        indx2 = r8vec_sort_heap_index_a(n2, r2);

        for (j2 = 0; j2 < n2; j2++)
        {
            unique2[j2] = true;
        }

        unique_num2 = 0;
        //
        //  STEP 2:
        //  Use PERMANENT points to eliminate TEMPORARY points.
        //
        for (j1 = 0; j1 < n1; j1++)
        {
            switch (unique1[indx1[j1]])
            {
                case true:
                {
                    double r_lo = r1[indx1[j1]] - tol;
                    double r_hi = r1[indx1[j1]] + tol;

                    r8vec_index_sorted_range(n2, r2, indx2, r_lo, r_hi,
                        ref j2_lo, ref j2_hi);

                    for (j2 = j2_lo; j2 <= j2_hi; j2++)
                    {
                        switch (unique2[indx2[j2]])
                        {
                            case true:
                            {
                                dist = 0.0;
                                for (i = 0; i < m; i++)
                                {
                                    dist += Math.Pow(a1[i + indx1[j1] * m]
                                                     - a2[i + indx2[j2] * m], 2);
                                }

                                dist = Math.Sqrt(dist);
                                if (dist <= tol)
                                {
                                    unique2[indx2[j2]] = false;
                                    xdnu2[indx2[j2]] = xdnu1[indx1[j1]];
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }

        //
        //  STEP 3:
        //  Use TEMPORARY points to eliminate TEMPORARY points.
        //
        for (j2 = 0; j2 < n2; j2++)
        {
            switch (unique2[indx2[j2]])
            {
                case true:
                {
                    xdnu2[indx2[j2]] = unique_num1 + unique_num2;
                    undx2[unique_num2] = indx2[j2] + n1;
                    unique_num2 += 1;

                    int hi = j2;

                    while (hi < n2 - 1)
                    {
                        if (r2[indx2[j2]] + tol < r2[indx2[hi + 1]])
                        {
                            break;
                        }

                        hi += 1;
                    }

                    int k2;
                    for (k2 = j2 + 1; k2 <= hi; k2++)
                    {
                        switch (unique2[indx2[k2]])
                        {
                            case true:
                            {
                                dist = 0.0;
                                for (i = 0; i < m; i++)
                                {
                                    dist += Math.Pow(a2[i + indx2[j2] * m] - a2[i + indx2[k2] * m], 2);
                                }

                                dist = Math.Sqrt(dist);

                                if (dist <= tol)
                                {
                                    unique2[indx2[k2]] = false;
                                    xdnu2[indx2[k2]] = xdnu2[indx2[j2]];
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }

    }

    public static void point_radial_tol_unique_index_inc3(int m, int n1, double[] a1,
            double[] r1, int[] indx1, bool[] unique1, int unique_num1, int[] undx1,
            int[] xdnu1, int n2, double[] a2, double[] r2, int[] indx2, bool[] unique2,
            int unique_num2, int[] undx2, int[] xdnu2, ref int n3, ref double[] a3, ref double[] r3,
            ref int[] indx3, ref bool[] unique3, ref int unique_num3, ref int[] undx3, ref int[] xdnu3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_RADIAL_TOL_UNIQUE_INDEX_INC3 merges index data.
        //
        //  Discussion:
        //
        //    This function may be called after *INDEX_INC1 has created index
        //    information for the permanent data, and *INDEX_INC2 has created
        //    augmenting information for a set of temporary data which now is
        //    to be merged with the permanent data.
        //
        //    The function merges the data and index information to create a
        //    new "permanent" data set.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N1, the number of permanent points.
        //
        //    Input, double A1[M*N1], the permanent points.
        //
        //    Input, double R1[N1], the scalar values assigned to
        //    the data for sorting.
        //
        //    Input, int INDX1[N1], the ascending sort index
        //    for A1.
        //
        //    Input, bool UNIQUE1[N1], is TRUE for each unique permanent point.
        //
        //    Input, int UNIQUE_NUM1, the number
        //    of tolerably unique points with just the permanent points.
        //
        //    Input, int UNDX1[UNIQUE_NUM1],
        //    the index in A1 of the tolerably unique points.
        //
        //    Input, int XDNU1[N1], the index in UNDX1
        //    of the tolerably unique point that "represents" this point.
        //
        //    Input, int N2, the number of temporary points.
        //
        //    Input, double A2[M,N2], the temporary points.
        //
        //    Input, double R2[N2], the scalar values assigned to
        //    the data for sorting.
        //
        //    Input, int INDX2[N2], the ascending sort index
        //    for A2.
        //
        //    Input, bool UNIQUE2[N2], is TRUE for each unique temporary point.
        //
        //    Input, int UNIQUE_NUM2, the number
        //    of tolerably unique temporary points.
        //
        //    Input, int UNDX2[UNIQUE_NUM2],
        //    the index in A2 of the tolerably unique points, incremented by UNIQUE_NUM1.
        //
        //    Input, int XDNU2[N2], the index in UNDX1 or UNDX2
        //    of the tolerably unique point that "represents" this point.
        //
        //    Output, int *N3, the number of permanent points.
        //
        //    Output, double A3[M,N3], the permanent points.
        //
        //    Output, double R3[N3], the scalar values assigned to
        //    the data for sorting.
        //
        //    Output, int INDX3[N3], the ascending sort index
        //    for A3.
        //
        //    Output, bool UNIQUE3[N3], is TRUE for each unique permanent point.
        //
        //    Output, int *UNIQUE_NUM3, the number
        //    of tolerably unique points.
        //
        //    Output, int UNDX3[UNIQUE_NUM3],
        //    the index in A3 of the tolerably unique points.
        //
        //    Output, int XDNU3[N3], the index in UNDX3
        //    of the tolerably unique point that "represents" this point.
        //
    {
        int i;
        int i1;
        int i2;
        int i3;

        n3 = n1 + n2;

        for (i1 = 0; i1 < n1; i1++)
        {
            for (i = 0; i < m; i++)
            {
                a3[i + i1 * m] = a1[i + i1 * m];
            }
        }

        for (i2 = 0; i2 < n2; i2++)
        {
            i3 = n1 + i2;
            for (i = 0; i < m; i++)
            {
                a3[i + i3 * m] = a2[i + i2 * m];
            }
        }

        for (i1 = 0; i1 < n1; i1++)
        {
            r3[i1] = r1[i1];
        }

        for (i2 = 0; i2 < n2; i2++)
        {
            i3 = n1 + i2;
            r3[i3] = r2[i2];
        }

        //
        //  Interleave the two INDX arrays so that INDX3 presents the entries
        //  of A3 in ascending R3 order.
        //
        i1 = 0;
        i2 = 0;

        for (i3 = 0; i3 < n3; i3++)
        {
            double v1 = i1 < n1 ? r1[indx1[i1]] : r8_huge();

            double v2 = i2 < n2 ? r2[indx2[i2]] : r8_huge();

            if (v1 <= v2)
            {
                indx3[i3] = indx1[i1];
                i1 += 1;
            }
            else
            {
                indx3[i3] = indx2[i2] + n1;
                i2 += 1;
            }
        }

        unique_num3 = unique_num1 + unique_num2;

        for (i1 = 0; i1 < n1; i1++)
        {
            unique3[i1] = unique1[i1];
        }

        for (i2 = 0; i2 < n2; i2++)
        {
            i3 = n1 + i2;
            unique3[i3] = unique2[i2];
        }

        //
        //  The entries in UNDX2 were already incremented by N2 if they pointed
        //  to an entry of A2, so all entries in UNDX2 correctly index A3.
        //
        for (i1 = 0; i1 < unique_num1; i1++)
        {
            undx3[i1] = undx1[i1];
        }

        for (i2 = 0; i2 < unique_num2; i2++)
        {
            i3 = unique_num1 + i2;
            undx3[i3] = undx2[i2];
        }

        //
        //  Note that the entries of XDNU2 were already incremented by N2
        //  so that they correctly index A3, not A2.
        //
        for (i1 = 0; i1 < n1; i1++)
        {
            xdnu3[i1] = xdnu1[i1];
        }

        for (i2 = 0; i2 < n2; i2++)
        {
            i3 = n1 + i2;
            xdnu3[i3] = xdnu2[i2];
        }

    }

    public static void point_unique_index(int m, int n, double[] a, int unique_num, ref int[] undx,
            ref int[] xdnu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_UNIQUE_INDEX indexes unique points.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The goal of this routine is to determine a vector UNDX,
        //    which points to the unique elements of A, in sorted order,
        //    and a vector XDNU, which identifies, for each entry of A, the index of
        //    the unique sorted element of A.
        //
        //    This is all done with index vectors, so that the elements of
        //    A are never moved.
        //
        //    The first step of the algorithm requires the indexed sorting
        //    of A, which creates arrays INDX and XDNI.  (If all the entries
        //    of A are unique, then these arrays are the same as UNDX and XDNU.)
        //
        //    We then use INDX to examine the entries of A in sorted order,
        //    noting the unique entries, creating the entries of XDNU and
        //    UNDX as we go.
        //
        //    Once this process has been completed, the vector A could be
        //    replaced by a compressed vector XU, containing the unique entries
        //    of A in sorted order, using the formula
        //
        //      XU(*) = A(UNDX(*)).
        //
        //    We could then, if we wished, reconstruct the entire vector A, or
        //    any element of it, by index, as follows:
        //
        //      A(I) = XU(XDNU(I)).
        //
        //    We could then replace A by the combination of XU and XDNU.
        //
        //    Later, when we need the I-th entry of A, we can locate it as
        //    the XDNU(I)-th entry of XU.
        //
        //    Here is an example of a vector A, the sort and inverse sort
        //    index vectors, and the unique sort and inverse unique sort vectors
        //    and the compressed unique sorted vector.
        //
        //      I     A  Indx  Xdni       XU  Undx  Xdnu
        //    ----+-----+-----+-----+--------+-----+-----+
        //      0 | 11.     0     0 |    11.     0     0
        //      1 | 22.     2     4 |    22.     1     1
        //      2 | 11.     5     1 |    33.     3     0
        //      3 | 33.     8     7 |    55.     4     2
        //      4 | 55.     1     8 |                  3
        //      5 | 11.     6     2 |                  0
        //      6 | 22.     7     5 |                  1
        //      7 | 22.     3     6 |                  1
        //      8 | 11.     4     3 |                  0
        //
        //    INDX(2) = 3 means that sorted item(2) is A(3).
        //    XDNI(2) = 5 means that A(2) is sorted item(5).
        //
        //    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
        //    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
        //
        //    XU(XDNU(I))) = A(I).
        //    XU(I)        = A(UNDX(I)).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the data values.
        //
        //    Input, int N, the number of data values,
        //
        //    Input, double A[M*N], the data values.
        //
        //    Input, int UNIQUE_NUM, the number of unique values in A.
        //    This value is only required for languages in which the size of
        //    UNDX must be known in advance.
        //
        //    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
        //
        //    Output, int XDNU[N], the XDNU vector.
        //
    {
        //
        //  Implicitly sort the array.
        //
        int[] indx = r8col_sort_heap_index_a(m, n, a);
        //
        //  Walk through the implicitly sorted array.
        //
        int i = 0;

        int j = 0;
        undx[j] = indx[i];

        xdnu[indx[i]] = j;

        for (i = 1; i < n; i++)
        {
            double diff = 0.0;
            int k;
            for (k = 0; k < m; k++)
            {
                diff = Math.Max(diff,
                    Math.Abs(a[k + indx[i] * m] - a[k + undx[j] * m]));
            }

            switch (diff)
            {
                case > 0.0:
                    j += 1;
                    undx[j] = indx[i];
                    break;
            }

            xdnu[indx[i]] = j;
        }
    }

}