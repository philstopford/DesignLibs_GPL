using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class NEWS
    {
        public static int[] news(int m, int n, int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NEWS demonstrates the NEWS stencil for edge detection.
            //
            //  Discussion:
            //
            //    Given a black and white image A, which we regard as an M by N array
            //    of pixels, we want to produce an array E of the same shape, which
            //    contains information describing the location of edges.
            //
            //    A simple algorithm for trying to detect edges in an array that
            //    represents an image is the NEWS scheme.  For each pixel A(C),
            //    we consider its North, East, West, and South pixel neighbors.  The
            //    indexing of arrays and images do not correspond, so we will use
            //    these directions instead:
            //
            //             A(N)
            //              |
            //              |
            //      A(W)---A(C)---A(E)
            //              |
            //              |
            //             A(S)
            //
            //    Entry E(C) of the edge array will be computed by
            //
            //      E(C) = abs ( A(N) - A(S) ) + abs ( A(E) - A(W) )
            //
            //    Pixels of A that represent edges will tend to have high values
            //    of E, while pixels that are interior to a region of roughly the
            //    same shade will tend to have low values.
            //
            //    Thus, an edge detection scheme would use the NEWS stencil to
            //    compute the E array, determine E_MAX, the maximum entry in E,
            //    choose some threshold value E_THRESH, and declare pixel A(I,J)
            //    to be associated with an edge whenever E(I,J) is greater than E_THRESH.
            //
            //    In this program, we demonstrate the NEWS stencil using a PGM
            //    grayscale image of coins.  At the end, we use the edge information
            //    to produce a color image in which the edges of the coins have been
            //    outlined in red.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the image.
            //
            //    Input, int A[M*N], the gray scale image data, presumably
            //    integers between 0 and 255.
            //
            //    Output, int NEWS[M*N], is 1 for each pixel that is part of an
            //    edge, and 0 otherwise.
            //
        {
            int[] b;
            int[] e;
            int e_max;
            int i;
            int j;
            int thresh;
            //
            //  For neatness, we add a border of zeros to the image,
            //  then fill in the border by copying the nearby original values.
            //  This will be our M+2 by N+2 data array B.
            //
            b = new int[(m + 2) * (n + 2)];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    b[i + 1 + (j + 1) * (m + 2)] = a[i + j * m];
                }
            }

            for (j = 1; j < n + 1; j++)
            {
                b[0 + j * (m + 2)] = b[1 + j * (m + 2)];
                b[m + 1 + j * (m + 2)] = b[m + j * (m + 2)];
            }

            for (i = 1; i < m + 1; i++)
            {
                b[i + 0 * (m + 2)] = b[i + 1 * (m + 2)];
                b[i + (n + 1) * (m + 2)] = b[i + n * (m + 2)];
            }

            b[0 + 0 * (m + 2)] = (b[0 + 1 * (m + 2)] + b[1 + 0 * (m + 2)]) / 2;
            b[m + 1 + 0 * (m + 2)] = (b[m + 1 + 1 * (m + 2)] + b[m + 0 + 0 * (m + 2)]) / 2;
            b[0 + (n + 1) * (m + 2)] = (b[0 + (n + 0) * (m + 2)] + b[1 + (n + 1) * (m + 2)]) / 2;
            b[m + 1 + (n + 1) * (m + 2)] = (b[m + 1 + (n + 0) * (m + 2)] + b[m + 0 + (n + 1) * (m + 2)]) / 2;
            //
            //  Apply the NEWS Operator.  We do not process the boundary pixels.
            //
            //  The picture is:
            //
            //   |  0 +1  0 |     |  0  0   0 |
            //   |  0  0  0 |  +  | -1  0  +1 |
            //   |  0 -1  0 |     |  0  0   0 |
            //
            e = new int[m * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    e[i + j * m] = Math.Abs(-b[i + (j + 1) * (m + 2)] + b[i + 2 + (j + 1) * (m + 2)])
                                   + Math.Abs(-b[i + 1 + j * (m + 2)] + b[i + 1 + (j + 2) * (m + 2)]);
                }
            }

            //
            //  Remap E so the largest value is 255.
            //
            e_max = typeMethods.i4mat_max(m, n, e);
            //
            //  Threshold the data.  Set the threshold to give enough detail
            //  to guess the coin denominations.
            //
            thresh = e_max / 5;

            Console.WriteLine("");
            Console.WriteLine("NEWS:");
            Console.WriteLine("  E_MAX = " + e_max + "");
            Console.WriteLine("  Using threshold value THRESH = " + thresh + "");

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (e[i + j * m] < thresh)
                    {
                        e[i + j * m] = 0;
                    }
                    else
                    {
                        e[i + j * m] = 1;
                    }
                }
            }

            return e;
        }
        

        public static int[] gray_median_news(int m, int n, int[] gray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_MEDIAN_NEWS uses a median NEWS filter on a gray scale image to remove noise.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 July 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns of pixels.
            //
            //    Input, int GRAY[M*N], the noisy grayscale data.
            //
            //    Output, int GRAY_MEDIAN_NEWS[M*N], the grayscale data for the filtered image.
            //
        {
            int[] gray2;
            int i;
            int j;
            int[] p = new int[5];

            gray2 = new int[m * n];
            //
            //  Process the main part of the image:
            //
            for (i = 1; i < m - 1; i++)
            {
                for (j = 1; j < n - 1; j++)
                {
                    p[0] = gray[i - 1 + j * m];
                    p[1] = gray[i + 1 + j * m];
                    p[2] = gray[i + (j + 1) * m];
                    p[3] = gray[i + (j - 1) * m];
                    p[4] = gray[i + j * m];

                    gray2[i + j * m] = typeMethods.i4vec_median(5, ref p);
                }
            }

            //
            //  Process the four borders.
            //  Get an odd number of data points, 
            //
            for (i = 1; i < m - 1; i++)
            {
                j = 0;
                p[0] = gray[i - 1 + j * m];
                p[1] = gray[i + 1 + j * m];
                p[2] = gray[i + j * m];
                p[3] = gray[i + (j + 1) * m];
                p[4] = gray[i + (j + 2) * m];
                gray2[i + j * m] = typeMethods.i4vec_median(5, ref p);

                j = n - 1;
                p[0] = gray[i - 1 + j * m];
                p[1] = gray[i + 1 + j * m];
                p[2] = gray[i + (j - 2) * m];
                p[3] = gray[i + (j - 1) * m];
                p[4] = gray[i + j * m];
                gray2[i + j * m] = typeMethods.i4vec_median(5, ref p);
            }

            for (j = 1; j < n - 1; j++)
            {
                i = 0;
                p[0] = gray[i + j * m];
                p[1] = gray[i + 1 + j * m];
                p[2] = gray[i + 2 + j * m];
                p[3] = gray[i + (j - 1) * m];
                p[4] = gray[i + (j + 1) * m];
                gray2[i + j * m] = typeMethods.i4vec_median(5, ref p);

                i = m - 1;
                p[0] = gray[i - 2 + j * m];
                p[1] = gray[i - 1 + j * m];
                p[2] = gray[i + j * m];
                p[3] = gray[i + (j - 1) * m];
                p[4] = gray[i + (j + 1) * m];
                gray2[i + j * m] = typeMethods.i4vec_median(5, ref p);
            }

            //
            //  Process the four corners.
            //
            i = 0;
            j = 0;
            p[0] = gray[i + 1 + j * m];
            p[1] = gray[i + j * m];
            p[2] = gray[i + (j + 1) * m];
            gray2[i + j * m] = typeMethods.i4vec_median(3, ref p);

            i = 0;
            j = n - 1;
            p[0] = gray[i + 1 + j * m];
            p[1] = gray[i + j * m];
            p[2] = gray[i + (j - 1) * m];
            gray2[i + j * m] = typeMethods.i4vec_median(3, ref p);

            i = m - 1;
            j = 0;
            p[0] = gray[i - 1 + j * m];
            p[1] = gray[i + j * m];
            p[2] = gray[i + (j + 1) * m];
            gray2[i + j * m] = typeMethods.i4vec_median(3, ref p);

            i = m - 1;
            j = n - 1;
            p[0] = gray[i - 1 + j * m];
            p[1] = gray[i + j * m];
            p[2] = gray[i + (j - 1) * m];
            gray2[i + j * m] = typeMethods.i4vec_median(3, ref p);

            return gray2;
        }
    }
}