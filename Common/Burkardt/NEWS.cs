using Burkardt.Types;

namespace Burkardt
{
    public static class NEWS
    {
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