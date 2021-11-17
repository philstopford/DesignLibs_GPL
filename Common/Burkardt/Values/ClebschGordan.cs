namespace Burkardt.Values;

public static class ClebschGordan
{
    public static void clebsch_gordan_values(ref int n_data, ref double j1, ref double j2, ref double j3,
            ref double m1, ref double m2, ref double m3, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLEBSCH_GORDAN_VALUES returns some values of the Clebsch-Gordan function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      ClebschGordan[{j1,m1},{j2,m2},{j3,m3}]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref double J1, &J2, &J3, &M1, &M2, &M3, the arguments
        //    of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 12;

        double[] fx_vec =
        {
            0.7071067811865475,
            1.000000000000000,
            0.5773502691896258,
            -0.2581988897471611,
            -0.6324555320336759,
            -0.7745966692414834,
            0.4082482904638630,
            0.8164965809277260,
            0.5345224838248488,
            0.2672612419124244,
            0.8944271909999159,
            0.3380617018914066
        };
        double[] j1_vec =
        {
            0.5,
            0.5,
            0.5,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            2.0,
            2.0,
            1.5,
            1.5
        };
        double[] j2_vec =
        {
            0.5,
            0.5,
            1.0,
            1.5,
            1.5,
            1.5,
            1.0,
            1.0,
            2.0,
            2.0,
            2.0,
            2.0
        };
        double[] j3_vec =
        {
            1.0,
            1.0,
            1.5,
            1.5,
            1.5,
            1.5,
            2.0,
            2.0,
            2.0,
            2.0,
            2.5,
            3.5
        };
        double[] m1_vec =
        {
            0.5,
            0.5,
            -0.5,
            0.0,
            -1.0,
            0.0,
            1.0,
            0.0,
            2.0,
            1.0,
            0.5,
            1.5
        };
        double[] m2_vec =
        {
            -0.5,
            0.5,
            1.0,
            0.5,
            1.5,
            1.5,
            -1.0,
            0.0,
            -2.0,
            -1.0,
            1.0,
            -1.0
        };
        double[] m3_vec =
        {
            -0.5,
            0.5,
            1.0,
            0.5,
            1.5,
            1.5,
            -1.0,
            0.0,
            -2.0,
            -1.0,
            1.0,
            -1.0
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            j1 = 0.0;
            j2 = 0.0;
            j3 = 0.0;
            m1 = 0.0;
            m2 = 0.0;
            m3 = 0.0;
            fx = 0.0;
        }
        else
        {
            j1 = j1_vec[n_data - 1];
            j2 = j2_vec[n_data - 1];
            j3 = j3_vec[n_data - 1];
            m1 = m1_vec[n_data - 1];
            m2 = m2_vec[n_data - 1];
            m3 = m3_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }
}