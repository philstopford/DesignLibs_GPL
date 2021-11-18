using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r82poly2_print(double a, double b, double c, double d, double e,
            double f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82POLY2_PRINT prints a second order polynomial in two variables.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, D, E, F, the coefficients.
        //
    {
        Console.WriteLine("  " + a.ToString().PadLeft(8)
                               + " * x^2 + " + b.ToString().PadLeft(8)
                               + " * y^2 + " + c.ToString().PadLeft(8)
                               + " * xy  + " + "");
        Console.WriteLine("  " + d.ToString().PadLeft(8)
                               + " * x   + " + e.ToString().PadLeft(8)
                               + " * y   + " + f.ToString().PadLeft(8) + "");

    }

    public static int r82poly2_type(double a, double b, double c, double d, double e, double f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82POLY2_TYPE analyzes a second order polynomial in two variables.
        //
        //  Discussion:
        //
        //    The polynomial has the form
        //
        //      A x^2 + B y^2 + C xy + Dx + Ey + F = 0
        //
        //    The possible types of the solution set are:
        //
        //     1: a hyperbola;
        //        9x^2 -  4y^2       -36x - 24y -  36 = 0
        //     2: a parabola;
        //        4x^2 +  1y^2 - 4xy + 3x -  4y +   1 = 0;
        //     3: an ellipse;
        //        9x^2 + 16y^2       +36x - 32y -  92 = 0;
        //     4: an imaginary ellipse (no real solutions);
        //         x^2 +   y^2       - 6x - 10y + 115 = 0;
        //     5: a pair of intersecting lines;
        //                        xy + 3x -   y -   3 = 0
        //     6: one point;
        //         x^2 +  2y^2       - 2x + 16y +  33 = 0;
        //     7: a pair of distinct parallel lines;
        //                 y^2            -  6y +   8 = 0
        //     8: a pair of imaginary parallel lines (no real solutions);
        //                 y^2            -  6y +  10 = 0
        //     9: a pair of coincident lines.
        //                 y^2            -  2y +   1 = 0
        //    10: a single line;
        //                             2x -   y +   1 = 0;
        //    11; all space;
        //                                          0 = 0;
        //    12; no solutions;
        //                                          1 = 0;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Daniel Zwillinger, editor,
        //    CRC Standard Mathematical Tables and Formulae,
        //    CRC Press, 30th Edition, 1996, pages 282-284.
        //
        //  Parameters:
        //
        //    Input, double A, B, C, D, E, F, the coefficients.
        //
        //    Output, int TYPE, indicates the type of the solution set.
        //
    {
        int type = 0;
        switch (a)
        {
            //
            //  Handle the degenerate case.
            //
            case 0.0 when b == 0.0 && c == 0.0:
            {
                switch (d)
                {
                    case 0.0 when e == 0.0:
                    {
                        type = f switch
                        {
                            0.0 => 11,
                            _ => 12
                        };

                        break;
                    }
                    default:
                        type = 10;
                        break;
                }

                return type;
            }
        }

        double delta = 8.0 * a * b * f
                       + 2.0 * c * e * d
                       - 2.0 * a * e * e
                       - 2.0 * b * d * d
                       - 2.0 * f * c * c;

        double j = 4.0 * a * b - c * c;

        if (delta != 0.0)
        {
            switch (j)
            {
                case < 0.0:
                    type = 1;
                    break;
                case 0.0:
                    type = 2;
                    break;
                case > 0.0 when Math.Abs(r8_sign(delta) - r8_sign(a + b)) > double.Epsilon:
                    type = 3;
                    break;
                case > 0.0:
                {
                    if (Math.Abs(r8_sign(delta) - r8_sign(a + b)) <= double.Epsilon)
                    {
                        type = 4;
                    }

                    break;
                }
            }
        }
        else
        {
            switch (delta)
            {
                case 0.0:
                    switch (j)
                    {
                        case < 0.0:
                            type = 5;
                            break;
                        case > 0.0:
                            type = 6;
                            break;
                        case 0.0:
                        {
                            double k = 4.0 * (a + b) * f - d * d - e * e;

                            type = k switch
                            {
                                < 0.0 => 7,
                                > 0.0 => 8,
                                0.0 => 9,
                                _ => type
                            };

                            break;
                        }
                    }

                    break;
            }
        }

        return type;
    }

    public static void r82poly2_type_print(int type)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82POLY2_TYPE_PRINT prints the meaning of the output from R82POLY2_TYPE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int TYPE, the type index returned by R82POLY2_TYPE.
        //
    {
        switch (type)
        {
            case 1:
                Console.WriteLine("  The set of solutions forms a hyperbola.");
                break;
            case 2:
                Console.WriteLine("  The set of solutions forms a parabola.");
                break;
            case 3:
                Console.WriteLine("  The set of solutions forms an ellipse.");
                break;
            case 4:
                Console.WriteLine("  The set of solutions forms an imaginary ellipse.");
                Console.WriteLine("  (There are no real solutions).");
                break;
            case 5:
                Console.WriteLine("  The set of solutions forms a pair of intersecting lines.");
                break;
            case 6:
                Console.WriteLine("  The set of solutions is a single point.");
                break;
            case 7:
                Console.WriteLine("  The set of solutions form a pair of distinct parallel lines.");
                break;
            case 8:
                Console.WriteLine("  The set of solutions forms a pair of imaginary parallel lines.");
                Console.WriteLine("  (There are no real solutions).");
                break;
            case 9:
                Console.WriteLine("  The set of solutions forms a pair of coincident lines.");
                break;
            case 10:
                Console.WriteLine("  The set of solutions forms a single line.");
                break;
            case 11:
                Console.WriteLine("  The set of solutions is all space.");
                break;
            case 12:
                Console.WriteLine("  The set of solutions is empty.");
                break;
            default:
                Console.WriteLine("  This type index is unknown.");
                break;
        }
    }
}