using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double rexp ( double x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    REXP evaluates the function EXP(X) - 1.
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
        //    double *X, the argument of the function.
        //
        //  Output:
        //
        //    double REXP, the value of EXP(X)-1.
        //
    {
        const double p1 = .914041914819518e-09;
        const double p2 = .238082361044469e-01;
        const double q1 = -.499999999085958e+00;
        const double q2 = .107141568980644e+00;
        const double q3 = -.119041179760821e-01;
        const double q4 = .595130811860248e-03;
        double rexp;

        switch (Math.Abs(x))
        {
            case > 0.15e0:
                goto S10;
        }
        rexp = x*(((p2*x+p1)*x+1.0e0)/((((q4*x+q3)*x+q2)*x+q1)*x+1.0e0));
        return rexp;
        S10:
        double w = Math.Exp(x);
        switch (x)
        {
            case > 0.0e0:
                goto S20;
        }
        rexp = w-0.5e0-0.5e0;
        return rexp;
        S20:
        rexp = w*(0.5e0+(0.5e0-1.0e0/w));
        return rexp;
    }

}