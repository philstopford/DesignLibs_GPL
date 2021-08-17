using System;

namespace StroudTest
{
    public static class mono
    {
        public static double mono_000_3d(int unused, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_000_3D evaluates X^0 Y^0 Z^0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the spatial dimension (which is 3 here).
            //
            //    Input, double X[N], the evaluation point.
            //
            //    Output, double MONO_000_3D, the value of the monomial at X.
            //
        {
            double value;

            value = 1.0;

            return value;
        }

        public static double mono_111_3d(int unused, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_111_3D evaluates X^1 Y^1 Z^1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the spatial dimension (which is 3 here).
            //
            //    Input, double X[N], the evaluation point.
            //
            //    Output, double MONO_111_3D, the value of the monomial at X.
            //
        {
            double value;

            value = Math.Pow(x[0], 1)
                    * Math.Pow(x[1], 1)
                    * Math.Pow(x[2], 1);

            return value;
        }

        public static double mono_202_3d(int unused, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_202_3D evaluates X^2 Y^0 Z^2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the spatial dimension (which is 3 here).
            //
            //    Input, double X[N], the evaluation point.
            //
            //    Output, double MONO_202_3D, the value of the monomial at X.
            //
        {
            double value;

            value = Math.Pow(x[0], 2)
                    * Math.Pow(x[1], 0)
                    * Math.Pow(x[2], 2);

            return value;
        }

        public static double mono_422_3d(int unused, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONO_422_3D evaluates X^4 Y^2 Z^2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the spatial dimension (which is 3 here).
            //
            //    Input, double X[N], the evaluation point.
            //
            //    Output, double MONO_422_3D, the value of the monomial at X.
            //
        {
            double value;

            value = Math.Pow(x[0], 4)
                    * Math.Pow(x[1], 2)
                    * Math.Pow(x[2], 2);

            return value;
        }
    }
}