namespace Burkardt.TestValues
{
    public static class Sphere
    {
        public static void sphere_unit_area_values(ref int n_data, ref int n, ref double area)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_AREA_VALUES returns some areas of the unit sphere in ND.
            //
            //  Discussion:
            //
            //    The formula for the surface area of the unit sphere in N dimensions is:
            //
            //      Sphere_Unit_Area ( N ) = 2 * PI^(N/2) / Gamma ( N / 2 )
            //
            //    Some values of the function include:
            //
            //       N   Area
            //
            //       2    2        * PI
            //       3  ( 4 /    ) * PI
            //       4  ( 2 /   1) * PI^2
            //       5  ( 8 /   3) * PI^2
            //       6  ( 1 /   1) * PI^3
            //       7  (16 /  15) * PI^3
            //       8  ( 1 /   3) * PI^4
            //       9  (32 / 105) * PI^4
            //      10  ( 1 /  12) * PI^5
            //
            //    For the unit sphere, Area(N) = N * Volume(N)
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      2 * Pi^(n/2) / Gamma[n/2]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
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
            //    Input/output, ref int N_DATA.
            //    On input, if N_DATA is 0, the first test data is returned, and
            //    N_DATA is set to the index of the test data.  On each subsequent
            //    call, N_DATA is incremented and that test data is returned.  When
            //    there is no more test data, N_DATA is set to 0.
            //
            //    Output, ref int N, the spatial dimension.
            //
            //    Output, ref double AREA, the area of the unit sphere
            //    in that dimension.
            //
        {
            int N_MAX = 20;

            double[] area_vec =
            {
                0.2000000000000000E+01,
                0.6283185307179586E+01,
                0.1256637061435917E+02,
                0.1973920880217872E+02,
                0.2631894506957162E+02,
                0.3100627668029982E+02,
                0.3307336179231981E+02,
                0.3246969701133415E+02,
                0.2968658012464836E+02,
                0.2550164039877345E+02,
                0.2072514267328890E+02,
                0.1602315322625507E+02,
                0.1183817381218268E+02,
                0.8389703410491089E+01,
                0.5721649212349567E+01,
                0.3765290085742291E+01,
                0.2396678817591364E+01,
                0.1478625959000308E+01,
                0.8858104195716824E+00,
                0.5161378278002812E+00
            };

            int[] n_vec =
            {
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                area = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                area = area_vec[n_data - 1];
            }
        }

        public static void sphere_unit_volume_values(ref int n_data, ref int n, ref double volume)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_VOLUME_VALUES returns some volumes of the unit sphere in ND.
            //
            //  Discussion:
            //
            //    The formula for the volume of the unit sphere in N dimensions is
            //
            //      Volume(N) = 2 * PI^(N/2) / ( N * Gamma ( N / 2 ) )
            //
            //    This function satisfies the relationships:
            //
            //      Volume(N) = 2 * PI * Volume(N-2) / N
            //      Volume(N) = Area(N) / N
            //
            //    Some values of the function include:
            //
            //       N  Volume
            //
            //       1    1
            //       2    1        * PI
            //       3  ( 4 /   3) * PI
            //       4  ( 1 /   2) * PI^2
            //       5  ( 8 /  15) * PI^2
            //       6  ( 1 /   6) * PI^3
            //       7  (16 / 105) * PI^3
            //       8  ( 1 /  24) * PI^4
            //       9  (32 / 945) * PI^4
            //      10  ( 1 / 120) * PI^5
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      2 * Pi^(n/2) / ( n * Gamma[n/2] )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
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
            //    Input/output, ref int N_DATA.
            //    On input, if N_DATA is 0, the first test data is returned, and
            //    N_DATA is set to the index of the test data.  On each subsequent
            //    call, N_DATA is incremented and that test data is returned.  When
            //    there is no more test data, N_DATA is set to 0.
            //
            //    Output, ref int N, the spatial dimension.
            //
            //    Output, ref double VOLUME, the volume of the unit
            //    sphere in that dimension.
            //
        {
            int N_MAX = 20;

            int[] n_vec =
            {
                1, 2,
                3, 4,
                5, 6,
                7, 8,
                9, 10,
                11, 12,
                13, 14,
                15, 16,
                17, 18,
                19, 20
            };

            double[] volume_vec =
            {
                0.2000000000000000E+01,
                0.3141592653589793E+01,
                0.4188790204786391E+01,
                0.4934802200544679E+01,
                0.5263789013914325E+01,
                0.5167712780049970E+01,
                0.4724765970331401E+01,
                0.4058712126416768E+01,
                0.3298508902738707E+01,
                0.2550164039877345E+01,
                0.1884103879389900E+01,
                0.1335262768854589E+01,
                0.9106287547832831E+00,
                0.5992645293207921E+00,
                0.3814432808233045E+00,
                0.2353306303588932E+00,
                0.1409811069171390E+00,
                0.8214588661112823E-01,
                0.4662160103008855E-01,
                0.2580689139001406E-01
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                volume = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                volume = volume_vec[n_data - 1];
            }
        }

    }
}