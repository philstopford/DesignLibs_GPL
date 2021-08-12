namespace Burkardt.TestValues
{
    public static class Debye
    {
        public static void debye1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DEBYE1_VALUES returns some values of Debye's function of order 1.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      DEBYE1(x) = 1 / x * integral ( 0 <= t <= x ) t / ( exp ( t ) - 1 ) dt
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.99951182471380889183E+00,
                0.99221462647120597836E+00,
                0.96918395997895308324E+00,
                0.88192715679060552968E+00,
                0.77750463411224827642E+00,
                0.68614531078940204342E+00,
                0.60694728460981007205E+00,
                0.53878956907785587703E+00,
                0.48043521957304283829E+00,
                0.38814802129793784501E+00,
                0.36930802829242526815E+00,
                0.32087619770014612104E+00,
                0.29423996623154246701E+00,
                0.27126046678502189985E+00,
                0.20523930310221503723E+00,
                0.16444346567994602563E+00,
                0.10966194482735821276E+00,
                0.82246701178200016086E-01,
                0.54831135561510852445E-01,
                0.32898681336964528729E-01
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void debye2_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DEBYE2_VALUES returns some values of Debye's function of order 2.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      DEBYE2(x) = 2 / x^2 * integral ( 0 <= t <= x ) t^2 / ( exp ( t ) - 1 ) dt
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.99934911727904599738E+00,
                0.98962402299599181205E+00,
                0.95898426200345986743E+00,
                0.84372119334725358934E+00,
                0.70787847562782928288E+00,
                0.59149637225671282917E+00,
                0.49308264399053185014E+00,
                0.41079413579749669069E+00,
                0.34261396060786351671E+00,
                0.24055368752127897660E+00,
                0.22082770061202308232E+00,
                0.17232915939014138975E+00,
                0.14724346738730182894E+00,
                0.12666919046715789982E+00,
                0.74268805954862854626E-01,
                0.47971498020121871622E-01,
                0.21369201683658373846E-01,
                0.12020564476446432799E-01,
                0.53424751249537071952E-02,
                0.19232910450553508562E-02
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void debye3_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DEBYE3_VALUES returns some values of Debye's function of order 3.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      DEBYE3(x) = 3 / x^3 * integral ( 0 <= t <= x ) t^3 / ( exp ( t ) - 1 ) dt
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.99926776885985461940E+00,
                0.98833007755734698212E+00,
                0.95390610472023510237E+00,
                0.82496296897623372315E+00,
                0.67441556407781468010E+00,
                0.54710665141286285468E+00,
                0.44112847372762418113E+00,
                0.35413603481042394211E+00,
                0.28357982814342246206E+00,
                0.18173691382177474795E+00,
                0.16277924385112436877E+00,
                0.11759741179993396450E+00,
                0.95240802723158889887E-01,
                0.77581324733763020269E-01,
                0.36560295673194845002E-01,
                0.19295765690345489563E-01,
                0.57712632276188798621E-02,
                0.24352200674805479827E-02,
                0.72154882216335666096E-03,
                0.15585454565440389896E-03
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void debye4_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DEBYE4_VALUES returns some values of Debye's function of order 4.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      DEBYE4(x) = 4 / x^4 * integral ( 0 <= t <= x ) t^4 / ( exp ( t ) - 1 ) dt
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.99921896192761576256E+00,
                0.98755425280996071022E+00,
                0.95086788606389739976E+00,
                0.81384569172034042516E+00,
                0.65487406888673697092E+00,
                0.52162830964878715188E+00,
                0.41189273671788528876E+00,
                0.32295434858707304628E+00,
                0.25187863642883314410E+00,
                0.15185461258672022043E+00,
                0.13372661145921413299E+00,
                0.91471377664481164749E-01,
                0.71227828197462523663E-01,
                0.55676547822738862783E-01,
                0.21967566525574960096E-01,
                0.96736755602711590082E-02,
                0.19646978158351837850E-02,
                0.62214648623965450200E-03,
                0.12289514092077854510E-03,
                0.15927210319002161231E-04
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }
    }
}