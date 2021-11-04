using System;

namespace Burkardt.Elliptic
{
    public static class EK_inc
    {
        public static double evaluate(double phi, double k)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_INC_EK evaluates the incomplete elliptic integral E(PHI,K).
            //
            //  Discussion:
            //
            //    The value is computed using Carlson elliptic integrals:
            //
            //      E(phi,k) = 
            //                  sin ( phi )   RF ( cos^2 ( phi ), 1-k^2 sin^2 ( phi ), 1 ) 
            //        - 1/3 k^2 sin^3 ( phi ) RD ( cos^2 ( phi ), 1-k^2 sin^2 ( phi ), 1 ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double PHI, K, the arguments.
            //    0 <= PHI <= PI/2.
            //    0 <= K^2 * sin^2(PHI) <= 1.
            //
            //    Output, double ELLIPTIC_INC_EK, the function value.
            //
        {
            double cp;
            double errtol;
            int ierr = 0;
            
            double sp;
            double value;
            double value1;
            double value2;
            double x;
            double y;
            double z;

            cp = Math.Cos(phi);
            sp = Math.Sin(phi);
            x = cp * cp;
            y = (1.0 - k * sp) * (1.0 + k * sp);
            z = 1.0;
            errtol = 1.0E-03;

            value1 = Integral.rf(x, y, z, errtol, ref ierr);

            if (ierr != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("ELLIPTIC_INC_EK - Fatal error!");
                Console.WriteLine("  RF returned IERR = " + ierr + "");
                return(1);
            }

            value2 = Integral.rd(x, y, z, errtol, ref ierr);

            if (ierr != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("ELLIPTIC_INC_EK - Fatal error!");
                Console.WriteLine("  RD returned IERR = " + ierr + "");
                return(1);
            }

            value = sp * value1 - k * k * sp * sp * sp * value2 / 3.0;

            return value;
        }

        public static void values(ref int n_data, ref double phi, ref double k, ref double ek )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_INC_EK_VALUES: values of the incomplete elliptic integral E(PHI,K).
        //
        //  Discussion:
        //
        //    This is the incomplete elliptic integral of the second kind.
        //
        //      E(PHI,K) = integral ( 0 <= T <= PHI ) 
        //        sqrt ( 1 - K^2 * sin ( T )^2 ) dT
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    US Department of Commerce, 1964.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Wolfram Media / Cambridge University Press, 1999.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 
        //    before the first call.  On each call, the routine increments N_DATA by 1, 
        //    and returns the corresponding data when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &PHI, &K, the arguments.
        //
        //    Output, double &EK, the function value.
        //
        {
            int N_MAX = 20;

            double[] ek_vec =
            {
                0.2852345328295404,
                1.298690225567921,
                0.5508100202571943,
                0.3575401358115371,
                0.06801307805507453,
                0.09679584980231837,
                0.6003112504412838,
                0.8996717721794724,
                1.380715261453875,
                0.1191644625202453,
                1.196994838171557,
                0.1536260979667945,
                0.3546768920544152,
                0.1758756066650882,
                1.229819109410569,
                1.08381066114337,
                1.35023378157378,
                1.419775884709218,
                0.2824895528020034,
                0.5770427720982867
            }
            ;

            double[] k_vec =
            {
                2.712952582080266,
                0.1279518954120547,
                -1.429437513650137,
                -1.981659235625333,
                3.894801879555818,
                -1.042486024983672,
                0.8641142168759754,
                -1.049058412826877,
                -0.3024062128402472,
                -6.574288841527263,
                0.6987397421988888,
                -5.12558591600033,
                2.074947853793764,
                -1.670886158426681,
                -0.4843595000931672,
                0.1393061679635559,
                -0.0946527302537008,
                0.1977207111754007,
                1.788159919089993,
                -1.077780624681256
            }
            ;

            double[] phi_vec =
            {
                0.3430906586047127,
                1.302990057703935,
                0.6523628380743488,
                0.4046022501376546,
                0.06884642871852312,
                0.0969609046794745,
                0.630370432896175,
                1.252375418911598,
                1.409796082144801,
                0.1485105463502483,
                1.349466184634646,
                0.1933711786970301,
                0.4088829927466769,
                0.1785430666405224,
                1.292588374416351,
                1.087095515757691,
                1.352794600489329,
                1.432530166308616,
                0.2968093345769761,
                0.6235880396594726
            }
            ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            if (N_MAX <= n_data)
            {
                n_data = 0;
                ek = 0.0;
                k = 0.0;
                phi = 0.0;
            }
            else
            {
                ek = ek_vec[n_data];
                k = k_vec[n_data];
                phi = phi_vec[n_data];
                n_data = n_data + 1;
            }
        }
    }
}