using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double alnrel(double a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ALNREL evaluates the function ln ( 1 + A ).
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
        //  Reference:
        //
        //    Armido DiDinato, Alfred Morris,
        //    Algorithm 708:
        //    Significant Digit Computation of the Incomplete Beta Function Ratios,
        //    ACM Transactions on Mathematical Software,
        //    Volume 18, 1993, pages 360-373.
        //
        //  Input:
        //
        //    double *A, the argument.
        //
        //  Output:
        //
        //    double ALNREL, the value of ln ( 1 + A ).
        //
    {
        double alnrel;
        const double p1 = -0.129418923021993e+01;
        const double p2 = 0.405303492862024e+00;
        const double p3 = -0.178874546012214e-01;
        const double q1 = -0.162752256355323e+01;
        const double q2 = 0.747811014037616e+00;
        const double q3 = -0.845104217945565e-01;

        switch (Math.Abs(a))
        {
            case <= 0.375e0:
                double t = a / (a + 2.0e0);
                double t2 = t * t;
                double w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.0e0)
                           / (((q3 * t2 + q2) * t2 + q1) * t2 + 1.0e0);
                alnrel = 2.0e0 * t * w;
                break;
            default:
                double x = 1.0e0 + a;
                alnrel = Math.Log(x);
                break;
        }

        return alnrel;
    }
}