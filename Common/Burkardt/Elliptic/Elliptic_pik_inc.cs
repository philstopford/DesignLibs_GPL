using System;

namespace Burkardt.Elliptic;

public static class PIK_inc
{
    public static double evaluate(double phi, double n, double k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_INC_PIK evaluates the incomplete elliptic integral Pi(PHI,N,K).
        //
        //  Discussion:
        //
        //    The value is computed using Carlson elliptic integrals:
        //
        //      Pi(PHI,N,K) = integral ( 0 <= T <= PHI )
        //        dT / (1 - N sin^2(T) ) sqrt ( 1 - k^2 * sin ( T )^2 )
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
        //    Input, double PHI, N, K, the arguments.
        //
        //    Output, double ELLIPTIC_INC_PIK, the function value.
        //
    {
        double cp;
        double errtol;
        int ierr = 0;
        double p;
            
        double sp;
        double value = 0;
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
        p = 1.0 - n * sp * sp;
        errtol = 1.0E-03;

        value1 = Integral.rf(x, y, z, errtol, ref ierr);

        if (ierr != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_PIK - Fatal error!");
            Console.WriteLine("  RF returned IERR = " + ierr + "");
            return 1;
        }

        value2 = Integral.rj(x, y, z, p, errtol, ref ierr);

        if (ierr != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_PIK - Fatal error!");
            Console.WriteLine("  RJ returned IERR = " + ierr + "");
            return 1;
        }

        value = sp * value1 + n * sp * sp * sp * value2 / 3.0;

        return value;
    }

    public static void values(ref int n_data, ref double phi, ref double n, ref double k,
            ref double pik )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_INC_PIK_VALUES: values of incomplete elliptic integral Pi(PHI,N,K).
        //
        //  Discussion:
        //
        //    This is the incomplete elliptic integral of the third kind.
        //
        //      Pi(PHI,N,K) = integral ( 0 <= T <= PHI ) 
        //        dT / (1 - N sin^2(T) ) sqrt ( 1 - K^2 * sin ( T )^2 )
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
        //    Output, double &PHI, &N, &K, the arguments of the function.
        //
        //    Output, double &PIK, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] k_vec =
            {
                1.959036804709882,
                -1.123741823223131,
                -2.317629084640271,
                -0.1202582658444815,
                1.008702896970963,
                -103.3677494756118,
                4.853800240677973,
                -1.016577251056124,
                -1.94341484065839,
                -0.8876593284500023,
                0.8160487832898813,
                0.2994546721661018,
                -0.7044232294525243,
                -0.9266523277404759,
                -0.6962608926846425,
                -0.4453932031991797,
                -0.9104582513322106,
                0.6187501419936026,
                0.8672305032589989,
                -0.1996772638241632
            }
            ;

        double[] n_vec =
            {
                8.064681366127422,
                -0.2840588974558835,
                -5.034023488967104,
                -1.244606253942751,
                1.465981775919188,
                95338.12857321106,
                -44.43130633436311,
                -0.8029374966926196,
                5.218883222649502,
                2.345821782626782,
                0.157358332363011,
                1.926593468907062,
                6.113982855261652,
                1.805710621498681,
                -0.4072847419780592,
                -0.9416404038595624,
                0.7009655305226739,
                -1.019830985340273,
                -0.4510798219577842,
                0.6028821390092596
            }
            ;

        double[] phi_vec =
            {
                0.3430906586047127,
                0.8823091382756705,
                0.4046022501376546,
                0.9958310121985398,
                0.630370432896175,
                0.002887706662908567,
                0.1485105463502483,
                1.320800086884777,
                0.4088829927466769,
                0.552337007372852,
                1.087095515757691,
                0.7128175949111615,
                0.2968093345769761,
                0.2910907344062498,
                0.9695030752034163,
                1.122288759723523,
                1.295911610809573,
                1.116491437736542,
                1.170719322533712,
                1.199360682338851
            }
            ;

        double[] pik_vec =
            {
                0.7982975462595892,
                1.024022134726036,
                0.40158120852642,
                0.7772649487439858,
                0.8737159913132074,
                0.004733334297691273,
                0.1280656893638068,
                1.594376037512564,
                0.8521145133671923,
                0.8154325229803082,
                1.31594514075427,
                1.25394623148424,
                0.3796503567258643,
                0.3111034454739552,
                0.9442477901112342,
                0.9153111661980959,
                2.842080644328393,
                0.9263253777034376,
                1.212396018757624,
                1.628083572710471
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
            k = 0.0;
            n = 0.0;
            phi = 0.0;
            pik = 0.0;
        }
        else
        {
            k = k_vec[n_data];
            n = n_vec[n_data];
            phi = phi_vec[n_data];
            pik = pik_vec[n_data];
            n_data += 1;
        }
    }
}