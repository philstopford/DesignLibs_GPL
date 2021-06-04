using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double algdiv(double a, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ALGDIV computes ln ( Gamma ( B ) / Gamma ( A + B ) ) when 8 <= B.
            //
            //  Discussion:
            //
            //    In this algorithm, DEL(X) is the function defined by
            //
            //      ln ( Gamma(X) ) = ( X - 0.5 ) * ln ( X ) - X + 0.5 * ln ( 2 * PI )
            //                      + DEL(X).
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
            //    double *A, *B, define the arguments.
            //
            //  Output:
            //
            //    double ALGDIV, the value of ln(Gamma(B)/Gamma(A+B)).
            //
        {
            double algdiv;
            double c;
            double c0 = 0.833333333333333e-01;
            double c1 = -0.277777777760991e-02;
            double c2 = 0.793650666825390e-03;
            double c3 = -0.595202931351870e-03;
            double c4 = 0.837308034031215e-03;
            double c5 = -0.165322962780713e-02;
            double d;
            double h;
            double s11;
            double s3;
            double s5;
            double s7;
            double s9;
            double t;
            double T1;
            double u;
            double v;
            double w;
            double x;
            double x2;

            if (b <= a)
            {
                h = b / a;
                c = 1.0e0 / (1.0e0 + h);
                x = h / (1.0e0 + h);
                d = a + (b - 0.5e0);
            }
            else
            {
                h = a / b;
                c = h / (1.0e0 + h);
                x = 1.0e0 / (1.0e0 + h);
                d = b + (a - 0.5e0);
            }

            //
            //  SET SN = (1 - X**N)/(1 - X)
            //
            x2 = x * x;
            s3 = 1.0e0 + (x + x2);
            s5 = 1.0e0 + (x + x2 * s3);
            s7 = 1.0e0 + (x + x2 * s5);
            s9 = 1.0e0 + (x + x2 * s7);
            s11 = 1.0e0 + (x + x2 * s9);
            //
            //  SET W = DEL(B) - DEL(A + B)
            //
            t = Math.Pow(1.0e0 / b, 2.0);

            w = ((((c5 * s11 * t
                    + c4 * s9) * t
                   + c3 * s7) * t
                  + c2 * s5) * t
                 + c1 * s3) * t
                + c0;

            w = w * (c / b);
            //
            //  Combine the results.
            //
            T1 = a / b;
            u = d * alnrel(T1);
            v = a * (Math.Log(b) - 1.0e0);

            if (v < u)
            {
                algdiv = w - v - u;
            }
            else
            {
                algdiv = w - u - v;
            }

            return algdiv;

        }
    }
}