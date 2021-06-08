namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double fifdsign ( double mag, double sign )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FIFDSIGN transfers the sign of the variable "sign" to the variable "mag"
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
            //  mag     -     magnitude
            //  sign    -     sign to be transfered
            //
        {
            if (mag < 0) mag = -mag;
            if (sign < 0) mag = -mag;
            return mag;
        }
    }
}