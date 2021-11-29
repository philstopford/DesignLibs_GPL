using System;

namespace Burkardt.Elliptic;

public static class PIM_inc
{
    public static double evaluate(double phi, double n, double m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_INC_PIM evaluates the incomplete elliptic integral Pi(PHI,N,M).
        //
        //  Discussion:
        //
        //    The value is computed using Carlson elliptic integrals:
        //
        //      Pi(PHI,N,M) = integral ( 0 <= T <= PHI )
        //        dT / (1 - N sin^2(T) ) sqrt ( 1 - m * sin ( T )^2 )
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
        //    Input, double PHI, N, M, the arguments.
        //
        //    Output, double ELLIPTIC_INC_PIM, the function value.
        //
    {
        int ierr = 0;

        double cp = Math.Cos(phi);
        double sp = Math.Sin(phi);
        double x = cp * cp;
        double y = 1.0 - m * sp * sp;
        const double z = 1.0;
        double p = 1.0 - n * sp * sp;
        const double errtol = 1.0E-03;

        double value1 = Integral.rf(x, y, z, errtol, ref ierr);

        if (ierr != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_PIM - Fatal error!");
            Console.WriteLine("  RF returned IERR = " + ierr + "");
            return 1;
        }

        double value2 = Integral.rj(x, y, z, p, errtol, ref ierr);

        if (ierr != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_PIM - Fatal error!");
            Console.WriteLine("  RJ returned IERR = " + ierr + "");
            return 1;
        }

        double value = sp * value1 + n * sp * sp * sp * value2 / 3.0;

        return value;
    }

    public static void values(ref int n_data, ref double phi, ref double n, ref double m,
            ref double pim )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_INC_PIM_VALUES: values of incomplete elliptic integral Pi(PHI,N,M).
        //
        //  Discussion:
        //
        //    This is the incomplete elliptic integral of the third kind.
        //
        //      Pi(PHI,N,M) = integral ( 0 <= T <= PHI ) 
        //        dT / (1 - N sin^2(T) ) sqrt ( 1 - M * sin ( T )^2 )
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
        //    Output, double &PHI, &N, &M, the arguments of the function.
        //
        //    Output, double &PIM, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] m_vec =
            {
                7.330122710928245,
                0.1108806690614566,
                0.2828355944410993,
                0.6382999794812498,
                2.294718938593894,
                42062.55329826538,
                39.2394337789563,
                0.008002151065098688,
                0.7190579590867517,
                0.9703767630929055,
                1.098881295982823,
                1.398066725917478,
                4.641021931654496,
                4.455969064311461,
                0.3131448239736511,
                0.3686443684703166,
                0.06678210908100803,
                0.9635538974026796,
                1.060208762696207,
                0.4687160847955397
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

        double[] pim_vec =
            {
                1.0469349800785,
                0.842114448140669,
                0.3321642201520043,
                0.8483033529960849,
                1.055753817656772,
                0.005108896144265593,
                0.1426848042785896,
                1.031350958206424,
                0.7131013701418496,
                0.8268044665355507,
                1.57632867896015,
                1.542817120857211,
                0.4144629799126912,
                0.3313231611366746,
                0.9195822851915201,
                0.9422320754002217,
                2.036599002815859,
                1.076799231499882,
                1.416084462957852,
                1.824124922310891
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
            m = 0.0;
            n = 0.0;
            phi = 0.0;
            pim = 0.0;
        }
        else
        {
            m = m_vec[n_data];
            n = n_vec[n_data];
            phi = phi_vec[n_data];
            pim = pim_vec[n_data];
            n_data += 1;
        }
    }
}