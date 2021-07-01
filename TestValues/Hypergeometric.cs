namespace TestValues
{
    public static class Hypergeometric
    {
        public static void hyper_1f1_values(ref int n_data, ref double a, ref double b, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HYPER_1F1_VALUES returns some values of the hypergeometric function 1F1.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      fx = Hypergeometric1F1 [ a, b, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 March 2010
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
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition, CRC Press, 1996, pages 651-652.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double A, &B, &X, the parameters of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 24;

            double[] a_vec =
            {
                -2.500,
                -0.500,
                0.500,
                2.500,
                -2.500,
                -0.500,
                0.500,
                2.500,
                -2.500,
                -0.500,
                0.500,
                2.500,
                0.825,
                1.100,
                1.650,
                3.300,
                0.825,
                1.100,
                1.650,
                3.300,
                0.825,
                1.100,
                1.650,
                3.300
            };
            double[] b_vec =
            {
                3.3,
                1.1,
                1.1,
                3.3,
                3.3,
                1.1,
                1.1,
                3.3,
                3.3,
                1.1,
                1.1,
                3.3,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7
            };
            double[] fx_vec =
            {
                0.81879926689265186854,
                0.88283984828032972070,
                1.1245023764952626690,
                1.2101049301639599598,
                0.12723045536781567174,
                0.12326016871544045107,
                2.3297954665128293051,
                3.3890020264468009733,
                -0.18819510282516768874,
                -1.0764203806547022727,
                5.7521824680907968433,
                9.9998567403304086593,
                1.0317208964319891384,
                1.0424867029249952040,
                1.0643112000949092012,
                1.1321844369742336326,
                1.2328402688568452181,
                1.3200654482027340732,
                1.5104811522310825217,
                2.2307520785940524365,
                1.5197286298183137741,
                1.7364938170250847619,
                2.2492330307668135926,
                4.6377737119178965298
            };
            double[] x_vec =
            {
                0.25,
                0.25,
                0.25,
                0.25,
                1.55,
                1.55,
                1.55,
                1.55,
                2.85,
                2.85,
                2.85,
                2.85,
                0.25,
                0.25,
                0.25,
                0.25,
                1.55,
                1.55,
                1.55,
                1.55,
                2.85,
                2.85,
                2.85,
                2.85
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0.0;
                b = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void hyper_2f1_values(ref int n_data, ref double a, ref double b, ref double c,
                ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HYPER_2F1_VALUES returns some values of the hypergeometric function 2F1.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      fx = Hypergeometric2F1 [ a, b, c, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 September 2007
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
            //    Shanjie Zhang, Jianming Jin,
            //    Computation of Special Functions,
            //    Wiley, 1996,
            //    ISBN: 0-471-11963-6,
            //    LC: QA351.C45
            //
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition, CRC Press, 1996, pages 651-652.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double A, &B, &C, &X, the parameters of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 24;

            double[] a_vec =
            {
                -2.5,
                -0.5,
                0.5,
                2.5,
                -2.5,
                -0.5,
                0.5,
                2.5,
                -2.5,
                -0.5,
                0.5,
                2.5,
                3.3,
                1.1,
                1.1,
                3.3,
                3.3,
                1.1,
                1.1,
                3.3,
                3.3,
                1.1,
                1.1,
                3.3
            };
            double[] b_vec =
            {
                3.3,
                1.1,
                1.1,
                3.3,
                3.3,
                1.1,
                1.1,
                3.3,
                3.3,
                1.1,
                1.1,
                3.3,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7
            };
            double[] c_vec =
            {
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                -5.5,
                -0.5,
                0.5,
                4.5,
                -5.5,
                -0.5,
                0.5,
                4.5,
                -5.5,
                -0.5,
                0.5,
                4.5
            };
            double[] fx_vec =
            {
                0.72356129348997784913,
                0.97911109345277961340,
                1.0216578140088564160,
                1.4051563200112126405,
                0.46961431639821611095,
                0.95296194977446325454,
                1.0512814213947987916,
                2.3999062904777858999,
                0.29106095928414718320,
                0.92536967910373175753,
                1.0865504094806997287,
                5.7381565526189046578,
                15090.669748704606754,
                -104.31170067364349677,
                21.175050707768812938,
                4.1946915819031922850,
                1.0170777974048815592E+10,
                -24708.635322489155868,
                1372.2304548384989560,
                58.092728706394652211,
                5.8682087615124176162E+18,
                -4.4635010147295996680E+08,
                5.3835057561295731310E+06,
                20396.913776019659426
            };
            double[] x_vec =
            {
                0.25,
                0.25,
                0.25,
                0.25,
                0.55,
                0.55,
                0.55,
                0.55,
                0.85,
                0.85,
                0.85,
                0.85,
                0.25,
                0.25,
                0.25,
                0.25,
                0.55,
                0.55,
                0.55,
                0.55,
                0.85,
                0.85,
                0.85,
                0.85
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0.0;
                b = 0.0;
                c = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                c = c_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void hypergeometric_cdf_values(ref int n_data, ref int sam, ref int suc, ref int pop,
                ref int n, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
            //
            //  Discussion:
            //
            //    CDF(X)(A,B) is the probability of at most X successes in A trials,
            //    given that the probability of success on a single trial is B.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`DiscreteDistributions`]
            //      dist = HypergeometricDistribution [ sam, suc, pop ]
            //      CDF [ dist, n ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
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
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition, CRC Press, 1996, pages 651-652.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int SAM, ref int SUC, ref int POP, the sample size,
            //    success size, and population parameters of the function.
            //
            //    Output, ref int N, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 16;

            double[] fx_vec =
            {
                0.6001858177500578E-01,
                0.2615284665839845E+00,
                0.6695237889132748E+00,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.5332595856827856E+00,
                0.1819495964117640E+00,
                0.4448047017527730E-01,
                0.9999991751316731E+00,
                0.9926860896560750E+00,
                0.8410799901444538E+00,
                0.3459800113391901E+00,
                0.0000000000000000E+00,
                0.2088888139634505E-02,
                0.3876752992448843E+00,
                0.9135215248834896E+00
            };

            int[] n_vec =
            {
                7, 8, 9, 10,
                6, 6, 6, 6,
                6, 6, 6, 6,
                0, 0, 0, 0
            };

            int[] pop_vec =
            {
                100, 100, 100, 100,
                100, 100, 100, 100,
                100, 100, 100, 100,
                90, 200, 1000, 10000
            };

            int[] sam_vec =
            {
                10, 10, 10, 10,
                6, 7, 8, 9,
                10, 10, 10, 10,
                10, 10, 10, 10
            };

            int[] suc_vec =
            {
                90, 90, 90, 90,
                90, 90, 90, 90,
                10, 30, 50, 70,
                90, 90, 90, 90
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                sam = 0;
                suc = 0;
                pop = 0;
                n = 0;
                fx = 0.0;
            }
            else
            {
                sam = sam_vec[n_data - 1];
                suc = suc_vec[n_data - 1];
                pop = pop_vec[n_data - 1];
                n = n_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void hypergeometric_pdf_values(ref int n_data, ref int sam, ref int suc, ref int pop,
                ref int n, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HYPERGEOMETRIC_PDF_VALUES returns some values of the hypergeometric PDF.
            //
            //  Discussion:
            //
            //    CDF(X)(A,B) is the probability of X successes in A trials,
            //    given that the probability of success on a single trial is B.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      dist = HypergeometricDistribution [ sam, suc, pop ]
            //      PDF [ dist, n ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 January 2008
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
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition, CRC Press, 1996, pages 651-652.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int SAM, ref int SUC, ref int POP, the sample size,
            //    success size, and population parameters of the function.
            //
            //    Output, ref int N, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 16;

            double[] fx_vec =
            {
                0.05179370533242827E+00,
                0.2015098848089788E+00,
                0.4079953223292903E+00,
                0.3304762110867252E+00,
                0.5223047493549780E+00,
                0.3889503452643453E+00,
                0.1505614239732950E+00,
                0.03927689321042477E+00,
                0.00003099828465518108E+00,
                0.03145116093938197E+00,
                0.2114132170316862E+00,
                0.2075776621999210E+00,
                0.0000000000000000E+00,
                0.002088888139634505E+00,
                0.3876752992448843E+00,
                0.9135215248834896E+00
            };

            int[] n_vec =
            {
                7, 8, 9, 10,
                6, 6, 6, 6,
                6, 6, 6, 6,
                0, 0, 0, 0
            };

            int[] pop_vec =
            {
                100, 100, 100, 100,
                100, 100, 100, 100,
                100, 100, 100, 100,
                90, 200, 1000, 10000
            };

            int[] sam_vec =
            {
                10, 10, 10, 10,
                6, 7, 8, 9,
                10, 10, 10, 10,
                10, 10, 10, 10
            };

            int[] suc_vec =
            {
                90, 90, 90, 90,
                90, 90, 90, 90,
                10, 30, 50, 70,
                90, 90, 90, 90
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                sam = 0;
                suc = 0;
                pop = 0;
                n = 0;
                fx = 0.0;
            }
            else
            {
                sam = sam_vec[n_data - 1];
                suc = suc_vec[n_data - 1];
                pop = pop_vec[n_data - 1];
                n = n_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void hypergeometric_u_values(ref int n_data, ref double a, ref double b, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HYPERGEOMETRIC_U_VALUES: some values of the hypergeometric function U(a,b,x).
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      fx = HypergeometricU [ a, b, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 October 2011
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
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition, CRC Press, 1996, pages 651-652.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double A, &B, &X, the parameters of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 24;

            double[] a_vec =
            {
                -2.500,
                -0.500,
                0.500,
                2.500,
                -2.500,
                -0.500,
                0.500,
                2.500,
                -2.500,
                -0.500,
                0.500,
                2.500,
                0.825,
                1.100,
                1.650,
                3.300,
                0.825,
                1.100,
                1.650,
                3.300,
                0.825,
                1.100,
                1.650,
                3.300
            };
            double[] b_vec =
            {
                3.3,
                1.1,
                1.1,
                3.3,
                3.3,
                1.1,
                1.1,
                3.3,
                3.3,
                1.1,
                1.1,
                3.3,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7,
                6.7
            };
            double[] fx_vec =
            {
                -68.693628728078601389,
                -0.0029710551374761070801,
                1.5008631742177797301,
                20.614688244200596134,
                7.4563815469305551938,
                1.0155793767749293733,
                0.73446538936622668912,
                0.28046404941879399225,
                3.4508153741446547607,
                1.5156637368753063495,
                0.56042118587934993510,
                0.064897147735134223341,
                223432.02356977463356,
                263079.25980740811495,
                269802.90319351274132,
                82809.311335606553425,
                26.465684783131844524,
                28.093506172516056560,
                23.889164624518872504,
                4.5338847857070388229,
                3.0224469362694842535,
                2.8040650913713359934,
                1.9262578111480172682,
                0.23020518115860909098
            };
            double[] x_vec =
            {
                0.25,
                0.25,
                0.25,
                0.25,
                1.55,
                1.55,
                1.55,
                1.55,
                2.85,
                2.85,
                2.85,
                2.85,
                0.25,
                0.25,
                0.25,
                0.25,
                1.55,
                1.55,
                1.55,
                1.55,
                2.85,
                2.85,
                2.85,
                2.85
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0.0;
                b = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}