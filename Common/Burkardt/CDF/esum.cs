using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double esum ( int mu, double x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ESUM evaluates exp ( MU + X ).
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
            //  Parameters:
            //
            //    Input, int *MU, part of the argument.
            //
            //    Input, double *X, part of the argument.
            //
            //    Output, double ESUM, the value of exp ( MU + X ).
            //
        {
            double esum,w;

            if(x > 0.0e0) goto S10;
            if(mu < 0) goto S20;
            w = (double)mu+x;
            if(w > 0.0e0) goto S20;
            esum = Math.Exp(w);
            return esum;
            S10:
            if(mu > 0) goto S20;
            w = (double)mu+x;
            if(w < 0.0e0) goto S20;
            esum = Math.Exp(w);
            return esum;
            S20:
            w = mu;
            esum = Math.Exp(w)*Math.Exp(x);
            return esum;
        }
    }
}