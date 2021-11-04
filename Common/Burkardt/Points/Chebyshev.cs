using System;

namespace Burkardt.PointsNS
{
    public static partial class Points
    {
                public static double[] chebyshev1(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBYSHEV1 returns the Type 1 Chebyshev points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2018
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double CHEBYSHEV1[N], the points.
            //
        {
            double angle;
            int i;
            
            double[] x;

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                angle = Math.PI * (double) (2 * i + 1) / (double) (2 * n);
                x[i] = Math.Cos(angle);
            }

            return x;
        }

        public static double[] chebyshev2(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBYSHEV2 returns the Type 2 Chebyshev points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2018
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double CHEBYSHEV2[N], the points.
            //
        {
            double angle;
            int i;
            
            double[] x;

            x = new double[n];

            if (n == 1)
            {
                x[0] = 0.0;
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    angle = Math.PI * (double) (n - i - 1) / (double) (n - 1);
                    x[i] = Math.Cos(angle);
                }
            }

            return x;
        }

        public static double[] chebyshev3(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBYSHEV3 returns the Type 3 Chebyshev points.
            //
            //  Discussion:
            //
            //    Note that this point set is NOT symmetric in [-1,+1].
            //    It is sometimes augmented by the value -1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2018
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double CHEBYSHEV3[N], the points.
            //
        {
            double angle;
            int i;
            
            double[] x;

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                angle = Math.PI * (double) (2 * n - 2 * i - 1)
                        / (double) (2 * n + 1);
                x[i] = Math.Cos(angle);
            }

            return x;
        }

        public static double[] chebyshev4(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBYSHEV4 returns the Type 4 Chebyshev points.
            //
            //  Discussion:
            //
            //    Note that this point set is NOT symmetric in [-1,+1].
            //    It is sometimes augmented by the value +1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 March 2018
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double CHEBYSHEV4[N], the points.
            //
        {
            double angle;
            int i;
            
            double[] x;

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                angle = Math.PI * (double) (2 * n - 2 * i)
                        / (double) (2 * n + 1);
                x[i] = Math.Cos(angle);
            }

            return x;
        }

    }
}