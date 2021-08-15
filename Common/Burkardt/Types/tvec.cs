namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] tvec_even(int nt)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TVEC_EVEN computes an evenly spaced set of angles between 0 and 2*PI.
            //
            //  Discussion:
            //
            //    The computation realizes that 0 = 2 * PI.
            //
            //  Example:
            //
            //    NT = 4
            //
            //    T = ( 0, PI/2, PI, 3*PI/2 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NT, the number of values to compute.
            //
            //    Output, double TVEC[NT], the evenly spaced angles, in radians.
            //
        {
            int i;
            double pi = 3.141592653589793;
            double[] t;

            if (nt < 1)
            {
                return null;
            }

            t = new double[nt];

            for (i = 1; i <= nt; i++)
            {
                t[i - 1] = (double)(2 * (i - 1)) * pi / (double)(nt);
            }

            return t;
        }

        public static double[] tvec_even2(int nt)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TVEC_EVEN2 computes evenly spaced angles between 0 and 2*PI.
            //
            //  Discussion:
            //
            //    The computation realizes that 0 = 2 * PI.  The values are equally
            //    spaced in the circle, do not include 0, and are symmetric about 0.
            //
            //  Example:
            //
            //    NT = 4
            //
            //    T = ( PI/4, 3*PI/4, 5*PI/4, 7*PI/4 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NT, the number of values to compute.
            //
            //    Output, double TVEC[NT], the evenly spaced angles, in radians.
            //
        {
            int i;
            double pi = 3.141592653589793;
            double[] t;

            if (nt < 1)
            {
                return null;
            }

            t = new double[nt];

            for (i = 1; i <= nt; i++)
            {
                t[i - 1] = (double)(2 * i - 1) * pi / (double)(nt);
            }

            return t;
        }

        public static double[] tvec_even3(int nt)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TVEC_EVEN3 computes an evenly spaced set of angles between 0 and 2*PI.
            //
            //  Discussion:
            //
            //    The angles begin with 0 and end with 2*PI.
            //
            //  Example:
            //
            //    NT = 4
            //
            //    T = ( 0, 2*PI/3, 4*PI/3 2*PI )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NT, the number of values to compute.
            //
            //    Output, double TVEC[NT], the evenly spaced angles, in radians.
            //
        {
            int i;
            double pi = 3.141592653589793;
            double[] t;

            if (nt < 1)
            {
                return null;
            }

            t = new double[nt];

            if (nt == 1)
            {
                t[0] = pi;
            }
            else
            {
                for (i = 1; i <= nt; i++)
                {
                    t[i - 1] = (double)(2 * (i - 1)) * pi / (double)(nt - 1);
                }
            }

            return t;
        }

        public static double[] tvec_even_bracket(int nt, double theta1, double theta2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TVEC_EVEN_BRACKET computes evenly spaced angles between THETA1 and THETA2.
            //
            //  Example:
            //
            //    NT = 4
            //    THETA1 = 30
            //    THETA2 = 90
            //
            //    T = ( 30, 50, 70, 90 )
            //
            //  Discussion:
            //
            //    The interval between THETA1 and THETA2 is divided into NT-1 subintervals.
            //
            //    The angles returned are the breakpoints of these subintervals,
            //    including THETA1 and THETA2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NT, the number of values to compute.
            //
            //    Input, double THETA1, THETA2, the limiting angles.
            //
            //    Output, double TVEC_EVEN_BRACKET[NT], the evenly spaced angles.
            //
        {
            int i;
            double[] t;

            if (nt < 1)
            {
                return null;
            }

            t = new double[nt];

            if (nt == 1)
            {
                t[0] = (theta1 + theta2) / 2.0;
            }
            else
            {
                for (i = 1; i <= nt; i++)
                {
                    t[i - 1] = ((double)(nt - i) * theta1
                                + (double)(i - 1) * theta2)
                               / (double)(nt - 1);
                }
            }

            return t;
        }

        public static double[] tvec_even_bracket2(int nt, double theta1, double theta2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TVEC_EVEN_BRACKET2 computes evenly spaced angles from THETA1 to THETA2.
            //
            //  Discussion:
            //
            //    The interval between THETA1 and THETA2 is divided into NT+1 subintervals.
            //
            //    The angles returned are the internal NT breakpoints of the subintervals.
            //
            //  Example:
            //
            //    NT = 5
            //    THETA1 = 30
            //    THETA2 = 90
            //
            //    T = ( 40, 50, 60, 70, 80 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NT, the number of values to compute.
            //
            //    Input, double THETA1, THETA2, the limiting angles.
            //
            //    Output, double TVEC_EVEN_BRACKET2[NT], the evenly spaced angles.
            //
        {
            int i;
            double[] t;

            if (nt < 1)
            {
                return null;
            }

            t = new double[nt];

            for (i = 1; i <= nt; i++)
            {
                t[i - 1] = ((double)(nt + 1 - i) * theta1
                            + (double)(i) * theta2)
                           / (double)(nt + 1);
            }

            return t;
        }

        public static double[] tvec_even_bracket3(int nt, double theta1, double theta2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TVEC_EVEN_BRACKET3 computes evenly spaced angles from THETA1 to THETA2.
            //
            //  Discussion:
            //
            //    The interval between THETA1 and THETA2 is divided into NT subintervals.
            //
            //    The angles returned are the midpoints of each subinterval.
            //
            //  Example:
            //
            //    NT = 3
            //    THETA1 = 30
            //    THETA2 = 90
            //
            //    T = ( 40, 60, 80 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NT, the number of values to compute.
            //
            //    Input, double THETA1, THETA2, the limiting angles.
            //
            //    Output, double TVEC_EVEN_BRACKET3[NT], the evenly spaced angles.
            //
        {
            int i;
            double[] t;

            t = new double[nt];

            for (i = 1; i <= nt; i++)
            {
                t[i - 1] = ((double)(2 * nt - 2 * i + 1) * theta1
                            + (double)(2 * i - 1) * theta2)
                           / (double)(2 * nt);
            }

            return t;
        }

    }
}