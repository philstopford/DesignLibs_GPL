using System;

namespace BlendTest;

internal static partial class Program
{
    private static double identity_r ( double r, int i )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    IDENTITY_R returns a data component given (R).
        //
        //  Discussion:
        //
        //    This is a dummy routine, which simply returns (R).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the coordinate of a point that lies on the
        //    boundary of the cube.
        //
        //    Input, int I, the component of the data which is to be returned.
        //
        //    Output, double *XI, the I-th component of the data vector X, evaluated
        //    at the point (R), which is on an endpoint of the unit line segment.
        //
    {
        double xi;
        switch (i)
        {
            case 0:
                xi = r;
                break;
            default:
                Console.WriteLine();;
                Console.WriteLine("IDENTITY_R - Fatal error!");
                Console.WriteLine("  Illegal component index I = " + i);
                xi = 0.0;
                break;
        }

        return xi;
    }

    private static double cubic_rs ( double r, double s, int i )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBIC_RS evaluates a function of R and S used for some tests.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double xi = 20.0 * ( r * r * s * s * s );

        return xi;
    }

    private static double quad_rst ( double r, double s, double t, int i )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUAD_RST evaluates a function of R, S and T used for some tests.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double xi = 18.0 * ( r * r + s + t );

        return xi;
    }


    private static double identity_rs ( double r, double s, int i )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    IDENTITY_RS returns a data component given (R,S).
        //
        //  Discussion:
        //
        //    This is a dummy routine, which simply returns (R,S).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, S, the coordinates of a point that lies on the
        //    boundary of the square.
        //
        //    Input, int I, the component of the data which is to be returned.
        //
        //    Output, double *XI, the I-th component of the data vector X, evaluated
        //    at the point (R,S), which is on a corner, or edge, of the unit square.
        //
    {
        double xi;
        switch (i)
        {
            case 0:
                xi = r;
                break;
            case 1:
                xi = s;
                break;
            default:
                Console.WriteLine();
                Console.WriteLine("IDENTITY_RS - Fatal error!");
                Console.WriteLine("  Illegal component index I = " + i);
                xi = 0.0;
                break;
        }

        return xi;
    }

    private static double identity_rst ( double r, double s, double t, int i )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    IDENTITY_RST returns a data component given (R,S,T).
        //
        //  Discussion:
        //
        //    This is a dummy routine, which simply returns (R,S,T).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, S, T, the coordinates of a point that lies on the
        //    boundary of the cube.
        //
        //    Input, int I, the component of the data which is to be returned.
        //
        //    Output, double *XI, the I-th component of the data vector X, evaluated
        //    at the point (R,S), which is on a corner, edge or face of the unit cube.
        //
    {
        double xi;
        switch (i)
        {
            case 0:
                xi = r;
                break;
            case 1:
                xi = s;
                break;
            case 2:
                xi = t;
                break;
            default:
                Console.WriteLine();
                Console.WriteLine("IDENTITY_RST - Fatal error!");
                Console.WriteLine("  Illegal component index I = " + i);
                xi = 0.0;
                break;
        }

        return xi;
    }


    private static double stretch_r ( double r, int i )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRETCH_R returns a data component given (R).
        //
        //  Discussion:
        //
        //    This routine shifts by 1 and stretches by 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the coordinate of a point that lies on the
        //    boundary of the cube.
        //
        //    Input, int I, the component of the data which is to be returned.
        //
        //    Output, double *XI, the I-th component of the data vector X, evaluated
        //    at the point (R), which is on an endpoint of the unit line segment.
        //
    {
        double xi;
        switch (i)
        {
            case 0:
                xi = 2.0 * r + 1.0;
                break;
            default:
                Console.WriteLine();
                Console.WriteLine("STRETCH_R - Fatal error!");
                Console.WriteLine("  Illegal component index I = " + i);
                xi = 0.0;
                break;
        }

        return xi;
    }


    private static double stretch_rs ( double r, double s, int i )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRETCH_RS returns a data component given (R,S).
        //
        //  Discussion:
        //
        //    This routine shifts by (1,2) and stretches by (3,4).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, S, the coordinates of a point that lies on the
        //    boundary of the square.
        //
        //    Input, int I, the component of the data which is to be returned.
        //
        //    Output, double *XI, the I-th component of the data vector X, evaluated
        //    at the point (R,S), which is on a corner, or edge, of the unit square.
        //
    {
        double xi;
        switch (i)
        {
            case 0:
                xi = 3.0 * r + 1.0;
                break;
            case 1:
                xi = 4.0 * s + 2.0;
                break;
            default:
                Console.WriteLine();
                Console.WriteLine("STRETCH_RS - Fatal error!");
                Console.WriteLine("  Illegal component index I = " + i);
                xi = 0.0;
                break;
        }

        return xi;
    }


    private static double stretch_rst ( double r, double s, double t, int i )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRETCH_RST returns a data component given (R,S,T).
        //
        //  Discussion:
        //
        //    This routine shifts by (1,2,3) and stretches by (4,5,6)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, S, T, the coordinates of a point that lies on the
        //    boundary of the cube.
        //
        //    Input, int I, the component of the data which is to be returned.
        //
        //    Output, double *XI, the I-th component of the data vector X, evaluated
        //    at the point (R,S), which is on a corner, edge or face of the unit cube.
        //
    {
        double xi;
        switch (i)
        {
            case 0:
                xi = 4.0 * r + 1.0;
                break;
            case 1:
                xi = 5.0 * s + 2.0;
                break;
            case 2:
                xi = 6.0 * t + 3.0;
                break;
            default:
                Console.WriteLine();
                Console.WriteLine("STRETCH_RST - Fatal error!");
                Console.WriteLine("  Illegal component index I = " + i);
                xi = 0.0;
                break;
        }

        return xi;
    }
}