using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
  public static double bcorr(double a0, double b0)

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    BCORR evaluates DEL(A0) + DEL(B0) - DEL(A0 + B0).
    //
    //  Discussion:
    //
    //    The function DEL(A) is a remainder term that is used in the expression:
    //
    //      ln ( Gamma ( A ) ) = ( A - 0.5 ) * ln ( A )
    //        - A + 0.5 * ln ( 2 * PI ) + DEL ( A ),
    //
    //    or, in other words, DEL ( A ) is defined as:
    //
    //      DEL ( A ) = ln ( Gamma ( A ) ) - ( A - 0.5 ) * ln ( A )
    //        + A + 0.5 * ln ( 2 * PI ).
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //
    //    14 February 2021
    //
    //  Author:
    //
    //    Barry Brown, James Lovato, Kathy Russell.
    //
    //  Input:
    //
    //    double *A0, *B0, the arguments.
    //    It is assumed that 8 <= A0 and 8 <= B0.
    //
    //  Output:
    //
    //    double BCORR, the value of the function.
    //
  {
    double a;
    double b;
    double bcorr;
    double c;
    double c0 = 0.833333333333333e-01;
    double c1 = -0.277777777760991e-02;
    double c2 = 0.793650666825390e-03;
    double c3 = -0.595202931351870e-03;
    double c4 = 0.837308034031215e-03;
    double c5 = -0.165322962780713e-02;
    double h;
    double s11;
    double s3;
    double s5;
    double s7;
    double s9;
    double t;
    double w;
    double x;
    double x2;

    a = Math.Min(a0, b0);
    b = Math.Max(a0, b0);
    h = a / b;
    c = h / (1.0e0 + h);
    x = 1.0e0 / (1.0e0 + h);
    x2 = x * x;
    //
    //  SET SN = (1 - X^N)/(1 - X)
    //
    s3 = 1.0e0 + (x + x2);
    s5 = 1.0e0 + (x + x2 * s3);
    s7 = 1.0e0 + (x + x2 * s5);
    s9 = 1.0e0 + (x + x2 * s7);
    s11 = 1.0e0 + (x + x2 * s9);
    //
    //  SET W = DEL(B) - DEL(A + B)
    //
    t = Math.Pow(1.0e0 / b, 2.0);

    w = ((((c5 * s11 * t + c4
            * s9) * t + c3
          * s7) * t + c2
        * s5) * t + c1
      * s3) * t + c0;
    w *= (c / b);
    //
    //  COMPUTE  DEL(A) + W
    //
    t = Math.Pow(1.0e0 / a, 2.0);

    bcorr = (((((c5 * t + c4)
            * t + c3)
          * t + c2)
        * t + c1)
      * t + c0) / a + w;
    return bcorr;
  }
}