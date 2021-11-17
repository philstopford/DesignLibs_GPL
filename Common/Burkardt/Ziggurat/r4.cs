using System;

namespace Burkardt.Ziggurat;

public static class r4
{
    public static float r4_exp(ref int jsr, int[] ke, float[] fe, float[] we)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4_EXP returns an exponentially distributed single precision real value.
        //
        //  Discussion:
        //
        //    The underlying algorithm is the ziggurat method.
        //
        //    Before the first call to this function, the user must call R4_EXP_SETUP
        //    to determine the values of KE, FE and WE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Marsaglia, Wai Wan Tsang,
        //    The Ziggurat Method for Generating Random Variables,
        //    Journal of Statistical Software,
        //    Volume 5, Number 8, October 2000, seven pages.
        //
        //  Parameters:
        //
        //    Input/output, uint32_t &JSR, the seed.
        //
        //    Input, uint32_t KE[256], data computed by R4_EXP_SETUP.
        //
        //    Input, float FE[256], WE[256], data computed by R4_EXP_SETUP.
        //
        //    Output, float R4_EXP, an exponentially distributed random value.
        //
    {
        float value;

        int jz = SHR3.shr3_seeded(ref jsr);
        int iz = jz & 255;

        if (jz < ke[iz])
        {
            value = jz * we[iz];
        }
        else
        {
            for (;;)
            {
                if (iz == 0)
                {
                    value = (float) (7.69711 - Math.Log(r4_uni(ref jsr)));
                    break;
                }

                float x = jz * we[iz];

                if (fe[iz] + r4_uni(ref jsr) * (fe[iz - 1] - fe[iz]) < Math.Exp(-x))
                {
                    value = x;
                    break;
                }

                jz = SHR3.shr3_seeded(ref jsr);
                iz = jz & 255;

                if (jz >= ke[iz])
                {
                    continue;
                }

                value = jz * we[iz];
                break;
            }
        }

        return value;
    }

    public static void r4_exp_setup(ref int[] ke, ref float[] fe, ref float[] we)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4_EXP_SETUP sets data needed by R4_EXP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Marsaglia, Wai Wan Tsang,
        //    The Ziggurat Method for Generating Random Variables,
        //    Journal of Statistical Software,
        //    Volume 5, Number 8, October 2000, seven pages.
        //
        //  Parameters:
        //
        //    Output, uint32_t KE[256], data needed by R4_EXP.
        //
        //    Output, float FE[256], WE[256], data needed by R4_EXP.
        //
    {
        double de = 7.697117470131487;
        int i;
        const double m2 = 2147483648.0;
        double te = 7.697117470131487;
        const double ve = 3.949659822581572E-03;

        double q = ve / Math.Exp(-de);

        ke[0] = (int) (de / q * m2);
        ke[1] = 0;

        we[0] = (float) (q / m2);
        we[255] = (float) (de / m2);

        fe[0] = 1.0f;
        fe[255] = (float) Math.Exp(-de);

        for (i = 254; 1 <= i; i--)
        {
            de = -Math.Log(ve / de + Math.Exp(-de));
            ke[i + 1] = (int) (de / te * m2);
            te = de;
            fe[i] = (float) Math.Exp(-de);
            we[i] = (float) (de / m2);
        }
    }

    public static float r4_nor(ref int jsr, int[] kn, float[] fn, float[] wn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4_NOR returns a normally distributed single precision real value.
        //
        //  Discussion:
        //
        //    The value returned is generated from a distribution with mean 0 and 
        //    variance 1.
        //
        //    The underlying algorithm is the ziggurat method.
        //
        //    Before the first call to this function, the user must call R4_NOR_SETUP
        //    to determine the values of KN, FN and WN.
        //
        //    Thanks to Chad Wagner, 21 July 2014, for noticing a bug of the form
        //      if ( x * x <= y * y );   <-- Stray semicolon!
        //      {
        //        break;
        //      }
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Marsaglia, Wai Wan Tsang,
        //    The Ziggurat Method for Generating Random Variables,
        //    Journal of Statistical Software,
        //    Volume 5, Number 8, October 2000, seven pages.
        //
        //  Parameters:
        //
        //    Input/output, uint32_t &JSR, the seed.
        //
        //    Input, uint32_t KN[128], data computed by R4_NOR_SETUP.
        //
        //    Input, float FN[128], WN[128], data computed by R4_NOR_SETUP.
        //
        //    Output, float R4_NOR, a normally distributed random value.
        //
    {
        const float r = 3.442620f;
        float value;

        int hz = SHR3.shr3_seeded(ref jsr);
        int iz = hz & 127;

        if (Math.Abs(hz) < kn[iz])
        {
            value = hz * wn[iz];
        }
        else
        {
            for (;;)
            {
                float x;
                if (iz == 0)
                {
                    for (;;)
                    {
                        x = (float) (-0.2904764 * Math.Log(r4_uni(ref jsr)));
                        float y = (float) -Math.Log(r4_uni(ref jsr));
                        if (x * x <= y + y)
                        {
                            break;
                        }
                    }

                    value = hz switch
                    {
                        <= 0 => -r - x,
                        _ => +r + x
                    };

                    break;
                }

                x = hz * wn[iz];

                if (fn[iz] + r4_uni(ref jsr) * (fn[iz - 1] - fn[iz])
                    < Math.Exp(-0.5 * x * x))
                {
                    value = x;
                    break;
                }

                hz = SHR3.shr3_seeded(ref jsr);
                iz = hz & 127;

                if (Math.Abs(hz) >= kn[iz])
                {
                    continue;
                }

                value = hz * wn[iz];
                break;
            }
        }

        return value;
    }

    public static void r4_nor_setup(ref int[] kn, ref float[] fn, ref float[] wn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4_NOR_SETUP sets data needed by R4_NOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Marsaglia, Wai Wan Tsang,
        //    The Ziggurat Method for Generating Random Variables,
        //    Journal of Statistical Software,
        //    Volume 5, Number 8, October 2000, seven pages.
        //
        //  Parameters:
        //
        //    Output, uint32_t KN[128], data needed by R4_NOR.
        //
        //    Output, float FN[128], WN[128], data needed by R4_NOR.
        //
    {
        double dn = 3.442619855899;
        int i;
        const double m1 = 2147483648.0;
        double tn = 3.442619855899;
        const double vn = 9.91256303526217E-03;

        double q = vn / Math.Exp(-0.5 * dn * dn);

        kn[0] = (int) (dn / q * m1);
        kn[1] = 0;

        wn[0] = (float) (q / m1);
        wn[127] = (float) (dn / m1);

        fn[0] = 1.0f;
        fn[127] = (float) Math.Exp(-0.5 * dn * dn);

        for (i = 126; 1 <= i; i--)
        {
            dn = Math.Sqrt(-2.0 * Math.Log(vn / dn + Math.Exp(-0.5 * dn * dn)));
            kn[i + 1] = (int) (dn / tn * m1);
            tn = dn;
            fn[i] = (float) Math.Exp(-0.5 * dn * dn);
            wn[i] = (float) (dn / m1);
        }
    }

    public static float r4_uni(ref int jsr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4_UNI returns a uniformly distributed real value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Marsaglia, Wai Wan Tsang,
        //    The Ziggurat Method for Generating Random Variables,
        //    Journal of Statistical Software,
        //    Volume 5, Number 8, October 2000, seven pages.
        //
        //  Parameters:
        //
        //    Input/output, uint32_t &JSR, the seed.
        //
        //    Output, float R4_UNI, a uniformly distributed random value in
        //    the range [0,1].
        //
    {
        int jsr_input = jsr;

        jsr ^= jsr << 13;
        jsr ^= jsr >> 17;
        jsr ^= jsr << 5;

        float value = (float) (0.5
                               + (float) (jsr_input + jsr) / 65536.0 / 65536.0 % 1.0);

        return value;
    }
}