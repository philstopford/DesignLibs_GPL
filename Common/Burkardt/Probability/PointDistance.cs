using System;
using Burkardt.Types;

namespace Burkardt.Probability;

public static class PointDistance
{
    public static double point_distance_1d_pdf(double x, int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_DISTANCE_1D_PDF evaluates the point distance PDF in 1D.
        //
        //  Discussion:
        //
        //    It is assumed that a set of points has been generated in 1D
        //    according to a Poisson process.  The number of points in a region
        //    of size LENGTH is a Poisson variate with mean value B * LENGTH.
        //
        //    For a point chosen at random, we may now find the nearest
        //    Poisson point, the second nearest and so on.  We are interested
        //    in the PDF that governs the expected behavior of the distances
        //    of rank A = 1, 2, 3, ... with Poisson density B.
        //
        //    Note that this PDF is a form of the Gamma PDF.???
        //
        //    PDF(A,B;X) = B^A * X^( A - 1 ) * EXP ( - B * X ) / ( A - 1 )!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    0.0 <= X.
        //
        //    Input, int A, indicates the degree of nearness of the point.
        //    A = 1 means the nearest point, A = 2 the second nearest, and so on.
        //    0 < A.
        //
        //    Input, double B, the point density.  0.0 < B.
        //
        //    Output, double PDF, the value of the PDF.
        //
    {
        double pdf;

        switch (a)
        {
            case < 1:
                Console.WriteLine(" ");
                Console.WriteLine("POINT_DISTANCE_1D_PDF - Fatal error!");
                Console.WriteLine("  Input parameter A < 1.");
                return 1;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("POINT_DISTANCE_1D_PDF - Fatal error!");
                Console.WriteLine("  Input parameter B <= 0.0.");
                return 1;
        }

        pdf = x switch
        {
            < 0.0 => 0.0,
            _ => Math.Pow(b, a) * Math.Pow(x, a - 1) * Math.Exp(-b * x) / typeMethods.r8_factorial(a - 1)
        };

        return pdf;
    }

    public static double point_distance_2d_pdf(double x, int a, double b)

//****************************************************************************80
//
//  Purpose:
//
//    POINT_DISTANCE_2D_PDF evaluates the point distance PDF in 2D.
//
//  Discussion:
//
//    It is assumed that a set of points has been generated in 2D
//    according to a Poisson process.  The number of points in a region
//    of size AREA is a Poisson variate with mean value B * AREA.
//
//    For a point chosen at random, we may now find the nearest
//    Poisson point, the second nearest and so on.  We are interested
//    in the PDF that governs the expected behavior of the distances
//    of rank A = 1, 2, 3, ... with Poisson density B.
//
//    PDF(A,B;X) = 2 * ( B * PI )^A * X^( 2 * A - 1 )
//      * EXP ( - B * PI * X * X ) / ( A - 1 )!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 579.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X.
//
//    Input, int A, indicates the degree of nearness of the point.
//    A = 1 means the nearest point, A = 2 the second nearest, and so on.
//    0 < A.
//
//    Input, double B, the point density.  0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
    {
        double pdf;
            

        switch (a)
        {
            case < 1:
                Console.WriteLine(" ");
                Console.WriteLine("POINT_DISTANCE_2D_PDF - Fatal error!");
                Console.WriteLine("  Input parameter A < 1.");
                return 1;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("POINT_DISTANCE_2D_PDF - Fatal error!");
                Console.WriteLine("  Input parameter B <= 0.0.");
                return 1;
        }

        pdf = x switch
        {
            < 0.0 => 0.0,
            _ => 2.0 * Math.Pow(b * Math.PI, a) * Math.Pow(x, 2 * a - 1) * Math.Exp(-b * Math.PI * x * x) /
                 typeMethods.r8_factorial(a - 1)
        };

        return pdf;
    }

    public static double point_distance_3d_pdf(double x, int a, double b)

//****************************************************************************80
//
//  Purpose:
//
//    POINT_DISTANCE_3D_PDF evaluates the point distance PDF in the 3D.
//
//  Discussion:
//
//    It is assumed that a set of points has been generated in 3D
//    according to a Poisson process.  The number of points in a region
//    of size VOLUME is a Poisson variate with mean value B * VOLUME.
//
//    For a point chosen at random, we may now find the nearest
//    Poisson point, the second nearest and so on.  We are interested
//    in the PDF that governs the expected behavior of the distances
//    of rank A = 1, 2, 3, ... with Poisson density B.
//
//    PDF(A,B;X) = 3 * ( (4/3) * B * PI )^A * X^( 3 * A - 1 )
//      * EXP ( - (4/3) * B * PI * X * X * X ) / ( A - 1 )!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 580.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X.
//
//    Input, int A, indicates the degree of nearness of the point.
//    A = 1 means the nearest point, A = 2 the second nearest, and so on.
//    0 < A.
//
//    Input, double B, the Poisson point density.  0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
    {
        double pdf;
            

        switch (a)
        {
            case < 1:
                Console.WriteLine(" ");
                Console.WriteLine("POINT_DISTANCE_3D_PDF - Fatal error!");
                Console.WriteLine("  Input parameter A < 1.");
                return 1;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("POINT_DISTANCE_3D_PDF - Fatal error!");
                Console.WriteLine("  Input parameter B <= 0.0.");
                return 1;
        }

        pdf = x switch
        {
            < 0.0 => 0.0,
            _ => 3.0 * Math.Pow(4.0 / 3.0 * b * Math.PI, a) * Math.Pow(x, 3 * a - 1) *
                Math.Exp(-(4.0 / 3.0) * b * Math.PI * x * x * x) / typeMethods.r8_factorial(a - 1)
        };

        return pdf;
    }
}