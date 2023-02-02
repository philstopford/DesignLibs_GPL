using System;
using Burkardt.Types;

namespace Burkardt.Geometry;

public static class LocalMinimum
{
    public static void minabs(double x1, double y1, double x2, double y2, double x3, double y3,
            ref double xmin, ref double ymin)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MINABS finds a local minimum of F(X) = A * abs ( X ) + B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, double X1, Y1, X2, Y2, X3, Y3, are three sets of data
        //    of the form ( X, F(X) ).  The three X values must be distinct.
        //    On output, the data has been sorted so that X1 < X2 < X3,
        //    and the Y values have been rearranged accordingly.
        //
        //    Output, double *XMIN, *YMIN.  XMIN is a point within the interval
        //    spanned by X1, X2 and X3, at which F takes its local minimum
        //    value YMIN.
        //
    {
        //
        //  Refuse to deal with coincident data.
        //
        if (Math.Abs(x1 - x2) <= typeMethods.r8_epsilon() || Math.Abs(x2 - x3) <= typeMethods.r8_epsilon() || Math.Abs(x3 - x1) <= typeMethods.r8_epsilon())
        {
            Console.WriteLine("");
            Console.WriteLine("MINABS - Fatal error!");
            Console.WriteLine("  X values are equal.");
            return;
        }

        //
        //  Sort the data.
        //
        if (x2 < x1)
        {
            typeMethods.r8_swap(ref x1, ref x2);
            typeMethods.r8_swap(ref y1, ref y2);
        }

        if (x3 < x1)
        {
            typeMethods.r8_swap(ref x1, ref x3);
            typeMethods.r8_swap(ref y1, ref y3);
        }

        if (x3 < x2)
        {
            typeMethods.r8_swap(ref x2, ref x3);
            typeMethods.r8_swap(ref y2, ref y3);
        }

        //
        //  Now determine the slopes.
        //
        double slope12 = (y2 - y1) / (x2 - x1);
        double slope23 = (y3 - y2) / (x3 - x2);
        double slope13 = (y3 - y1) / (x3 - x1);
        //
        //  Case 1: Minimum must be at an endpoint.
        //
        if (slope13 <= slope12 || 0.0 <= slope12)
        {
            if (y1 < y3)
            {
                xmin = x1;
                ymin = y1;
            }
            else
            {
                xmin = x3;
                ymin = y3;
            }
        }
        //
        //  Case 2: The curve decreases, and decreases faster than the line
        //  joining the endpoints.
        //
        //  Whichever of SLOPE12 and SLOPE23 is the greater in magnitude
        //  represents the actual slope of the underlying function.
        //  Find where two lines of that slope, passing through the
        //  endpoint data, intersect.
        //
        else
        {
            double slope = Math.Max(Math.Abs(slope12), slope23);
            xmin = 0.5 * (x1 + x3 + (y1 - y3) / slope);
            ymin = y1 - slope * (xmin - x1);
        }
    }

    public static bool minquad(double x1, double y1, double x2, double y2, double x3, double y3,
            ref double xmin, ref double ymin)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MINQUAD finds a local minimum of F(X) = A * X^2 + B * X + C.
        //
        //  Discussion:
        //
        //    MINQUAD is primarily intended as a utility routine for use by
        //    DISLSLS3.  The square of the distance function between a point
        //    and a line segment has the form of F(X).  Hence, we can seek
        //    the line on the second segment which minimizes the square of
        //    the distance to the other line segment.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 1998
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X1, Y1, X2, Y2, X3, Y3, are three sets of data
        //    of the form ( X, F(X) ).  The three X values must be distinct.
        //
        //    Output, double *XMIN, *YMIN.  XMIN is a point within the interval
        //    spanned by X1, X2 and X3, at which F takes its local minimum
        //    value YMIN.
        //
        //    Output, bool MINQUAD,
        //    true if no error,
        //    false if error because X values are not distinct.
        //
    {
        double x = 0;
        double y = 0;

        xmin = 0.0;
        ymin = 0.0;
        //
        //  Refuse to deal with coincident data.
        //
        if (Math.Abs(x1 - x2) <= typeMethods.r8_epsilon() || Math.Abs(x2 - x3) <= typeMethods.r8_epsilon() || Math.Abs(x3 - x1) <= typeMethods.r8_epsilon())
        {
            return false;
        }

        //
        //  Find the interval endpoints.
        //
        double xleft = x1;
        if (x2 < xleft)
        {
            xleft = x2;
        }

        if (x3 < xleft)
        {
            xleft = x3;
        }

        double xrite = x1;
        if (xrite < x2)
        {
            xrite = x2;
        }

        if (xrite < x3)
        {
            xrite = x3;
        }

        //
        //  Find the minimizer and its function value over the three input points.
        //
        if (y1 <= y2 && y1 <= y3)
        {
            xmin = x1;
            ymin = y1;
        }
        else if (y2 <= y1 && y2 <= y3)
        {
            xmin = x2;
            ymin = y2;
        }
        else if (y3 <= y1 && y3 <= y2)
        {
            xmin = x3;
            ymin = y3;
        }

        //
        //  Find the minimizer and its function value over the real line.
        //
        int ierror = ParabolaNS.Geometry.parabola_ex(x1, y1, x2, y2, x3, y3, ref x, ref y);

        if (ierror == 2 || !(y < ymin) || !(xleft < x) || !(x < xrite))
        {
            return true;
        }

        xmin = x;
        ymin = y;

        return true;
    }

}