using System;

namespace Burkardt.Elliptic;

public static class EM_inc
{
    public static double evaluate(double phi, double m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_INC_EM evaluates the incomplete elliptic integral E(PHI,M).
        //
        //  Discussion:
        //
        //    The value is computed using Carlson elliptic integrals:
        //
        //      E(phi,m) = 
        //                sin ( phi )   RF ( cos^2 ( phi ), 1-m sin^2 ( phi ), 1 ) 
        //        - 1/3 m sin^3 ( phi ) RD ( cos^2 ( phi ), 1-m sin^2 ( phi ), 1 ).
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
        //    Input, double PHI, M, the arguments.
        //    0 <= PHI <= PI/2.
        //    0 <= M * sin^2(PHI) <= 1.
        //
        //    Output, double ELLIPTIC_INC_EM, the function value.
        //
    {
        int ierr = 0;

        double cp = Math.Cos(phi);
        double sp = Math.Sin(phi);
        double x = cp * cp;
        double y = 1.0 - m * sp * sp;
        const double z = 1.0;
        const double errtol = 1.0E-03;

        double value1 = Integral.rf(x, y, z, errtol, ref ierr);

        if (ierr != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_EM - Fatal error!");
            Console.WriteLine("  RF returned IERR = " + ierr + "");
            return 1;
        }

        double value2 = Integral.rd(x, y, z, errtol, ref ierr);

        if (ierr != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_EM - Fatal error!");
            Console.WriteLine("  RD returned IERR = " + ierr + "");
            return 1;
        }

        double value = sp * value1 - m * sp * sp * sp * value2 / 3.0;

        return value;
    }

    public static void values(ref int n_data, ref double phi, ref double m, ref double em )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_INC_EM_VALUES: values of the incomplete elliptic integral E(PHI,M).
        //
        //  Discussion:
        //
        //    This is the incomplete elliptic integral of the second kind.
        //
        //      E(PHI,M) = integral ( 0 <= T <= PHI ) 
        //        sqrt ( 1 - M * sin ( T )^2 ) dT
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2018
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
        //    Output, double &PHI, &M, the arguments.
        //
        //    Output, double &EM, the function value.
        //
    {
        const int N_MAX = 20;

        double[] em_vec =
            {
                0.2732317284159052,
                1.124749725099781,
                0.6446601913679151,
                0.3968902354370061,
                0.06063960799944668,
                0.08909411577948728,
                0.532402014802015,
                1.251888640660265,
                1.28897116191626,
                0.1481718153599732,
                1.038090185639913,
                0.1931275771541276,
                0.3304419611986801,
                0.167394796063963,
                1.214501175324736,
                0.9516560179840655,
                1.203682959526176,
                1.206426326185419,
                0.2522791382096692,
                0.6026499038720986
            }
            ;

        double[] m_vec =
            {
                8.450689756874594,
                0.6039878267930615,
                0.1794126658351454,
                0.7095689301026752,
                133.9643389059188,
                47.96621393936416,
                2.172070586163255,
                0.002038130569431913,
                0.3600036705339421,
                0.6219544540067304,
                0.8834215943508453,
                0.2034290670379481,
                5.772526076430922,
                11.14853902343298,
                0.2889238477277305,
                0.7166617182589116,
                0.4760623731559658,
                0.6094948502068943,
                8.902276887883076,
                0.5434439226321253
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

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        if (N_MAX <= n_data)
        {
            n_data = 0;
            em = 0.0;
            m = 0.0;
            phi = 0.0;
        }
        else
        {
            em = em_vec[n_data];
            m = m_vec[n_data];
            phi = phi_vec[n_data];
            n_data += 1;
        }
    }
}