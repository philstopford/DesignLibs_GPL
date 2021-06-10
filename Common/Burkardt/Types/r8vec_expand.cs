namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8vec_expand_linear(int n, double[] x, int fat )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_EXPAND_LINEAR linearly interpolates new data into an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of input data values.
        //
        //    Input, double X[N], the original data.
        //
        //    Input, int FAT, the number of data values to interpolate
        //    between each pair of original data values.
        //
        //    Output, double R8VEC_EXPAND_LINEAR[(N-1)*(FAT+1)+1], the "fattened" data.
        //
        {
            int i;
            int j;
            int k;
            double[] xfat;

            xfat = new double[(n - 1) * (fat + 1) + 1];

            k = 0;

            for (i = 0; i < n - 1; i++)
            {
                xfat[k] = x[i];
                k = k + 1;

                for (j = 1; j <= fat; j++)
                {
                    xfat[k] = ((double) (fat - j + 1) * x[i]
                               + (double) (j) * x[i + 1])
                              / (double) (fat + 1);
                    k = k + 1;
                }
            }

            xfat[k] = x[n - 1];
            k = k + 1;

            return xfat;
        }

        public static double[] r8vec_expand_linear2(int n, double[] x, int before, int fat,
        int after )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_EXPAND_LINEAR2 linearly interpolates new data into an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    This routine starts with a vector of data.
        //
        //    The intent is to "fatten" the data, that is, to insert more points
        //    between successive values of the original data.
        //
        //    There will also be extra points placed BEFORE the first original
        //    value and AFTER that last original value.
        //
        //    The "fattened" data is equally spaced between the original points.
        //
        //    The BEFORE data uses the spacing of the first original interval,
        //    and the AFTER data uses the spacing of the last original interval.
        //
        //  Example:
        //
        //    N = 3
        //    BEFORE = 3
        //    FAT = 2
        //    AFTER = 1
        //
        //    X    = (/                   0.0,           6.0,             7.0       /)
        //    XFAT = (/ -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 6.33, 6.66, 7.0, 7.66 /)
        //            3 "BEFORE's"        Old  2 "FATS"  Old    2 "FATS"  Old  1 "AFTER"
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of input data values.
        //    N must be at least 2.
        //
        //    Input, double X[N], the original data.
        //
        //    Input, int BEFORE, the number of "before" values.
        //
        //    Input, int FAT, the number of data values to interpolate
        //    between each pair of original data values.
        //
        //    Input, int AFTER, the number of "after" values.
        //
        //    Output, double R8VEC_EXPAND_LINEAR2[BEFORE+(N-1)*(FAT+1)+1+AFTER], the
        //    "fattened" data.
        //
        {
            int i;
            int j;
            int k;
            double[] xfat;

            xfat = new double[before + (n - 1) * (fat + 1) + 1 + after];

            k = 0;
            //
            //  Points BEFORE.
            //
            for (j = 1 - before + fat; j <= fat; j++)
            {
                xfat[k] = ((double) (fat - j + 1) * (x[0] - (x[1] - x[0]))
                           + (double) (j) * x[0])
                          / (double) (fat + 1);
                k = k + 1;
            }

            //
            //  Original points and FAT points.
            //
            for (i = 0; i < n - 1; i++)
            {
                xfat[k] = x[0];
                k = k + 1;
                for (j = 1; j <= fat; j++)
                {
                    xfat[k] = ((double) (fat - j + 1) * x[i]
                               + (double) (j) * x[i + 1])
                              / (double) (fat + 1);
                    k = k + 1;
                }
            }

            xfat[k] = x[n - 1];
            k = k + 1;
            //
            //  Points AFTER.
            //
            for (j = 1; j <= after; j++)
            {
                xfat[k] = ((double) (fat - j + 1) * x[n - 1]
                           + (double) (j) * (x[n - 1] + (x[n - 1] - x[n - 2])))
                          / (double) (fat + 1);
                k = k + 1;
            }

            return xfat;
        }
    }
}