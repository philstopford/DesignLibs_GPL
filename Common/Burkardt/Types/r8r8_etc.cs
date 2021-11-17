using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static int r8r8_compare(double x1, double y1, double x2, double y2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8R8_COMPARE compares two R8R8's.
        //
        //  Discussion:
        //
        //    An R8R8 is simply a pair of R8 values, stored separately.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X1, Y1, the first vector.
        //
        //    Input, double X2, Y2, the second vector.
        //
        //    Output, int R8R8_COMPARE:
        //    -1, (X1,Y1) < (X2,Y2);
        //     0, (X1,Y1) = (X2,Y2);
        //    +1, (X1,Y1) > (X2,Y2).
        //
    {
        int value;

        if (x1 < x2)
        {
            value = -1;
        }
        else if (x2 < x1)
        {
            value = +1;
        }
        else if (y1 < y2)
        {
            value = -1;
        }
        else if (y2 < y1)
        {
            value = +1;
        }
        else
        {
            value = 0;
        }

        return value;
    }

    public static void r8r8_print(double a1, double a2, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8R8_PRINT prints an R8R8.
        //
        //  Discussion:
        //
        //    An R8R8 is a pair of R8 values, regarded as a single item.
        //
        //    A format is used which suggests a coordinate pair:
        //
        //  Example:
        //
        //    Center : ( 1.23, 7.45 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A1, A2, the coordinates of the vector.
        //
        //    Input, string TITLE, a title.
        //
    {
        string cout = "  " + title + " : ";
        cout += "  ( " + a1.ToString().PadLeft(12) + 
                ", " + a2.ToString().PadLeft(12) + " )";

        Console.WriteLine(cout);
    }

    public static int r8r8r8_compare(double x1, double y1, double z1, double x2, double y2,
            double z2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8R8R8_COMPARE compares two R8R8R8's.
        //
        //  Discussion:
        //
        //    An R8R8R8 is simply 3 R8 values, stored as scalars.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X1, Y1, Z1, the first vector.
        //
        //    Input, double X2, Y2, Z2, the second vector.
        //
        //    Output, int R8R8R8_COMPARE:
        //    -1, (X1,Y1,Z1) < (X2,Y2,Z2);
        //     0, (X1,Y1,Z1) = (X2,Y2,Z2);
        //    +1, (X1,Y1,Z1) > (X2,Y2,Z2).
        //
    {
        int value;

        if (x1 < x2)
        {
            value = -1;
        }
        else if (x2 < x1)
        {
            value = +1;
        }
        else if (y1 < y2)
        {
            value = -1;
        }
        else if (y2 < y1)
        {
            value = +1;
        }
        else if (z1 < z2)
        {
            value = -1;
        }
        else if (z2 < z1)
        {
            value = +1;
        }
        else
        {
            value = 0;
        }

        return value;
    }
}