namespace TestValues
{
    public static class I
    {

        public static void i0ml0_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I0ML0_VALUES returns some values of the I0ML0 function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      I0ML0(x) = I0(x) - L0(x)
            //
            //    I0(x) is the modified Bessel function of the first kind of order 0,
            //    L0(x) is the modified Struve function of order 0.
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 August 2004
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
                0.99875755515461749793E+00,
                0.99011358230706643807E+00,
                0.92419435310023947018E+00,
                0.73624267134714273902E+00,
                0.55582269181411744686E+00,
                0.34215154434462160628E+00,
                0.17087174888774706539E+00,
                0.81081008709219208918E-01,
                0.53449421441089580702E-01,
                0.39950321008923244846E-01,
                0.39330637437584921392E-01,
                0.37582274342808670750E-01,
                0.31912486554480390343E-01,
                0.25506146883504738403E-01,
                0.21244480317825292412E-01,
                0.15925498348551684335E-01,
                0.12737506927242585015E-01,
                0.84897750814784916847E-02,
                0.63668349178454469153E-02,
                0.50932843163122551114E-02
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0156250000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                2.0000000000E+00,
                4.0000000000E+00,
                8.0000000000E+00,
                12.0000000000E+00,
                16.0000000000E+00,
                16.2500000000E+00,
                17.0000000000E+00,
                20.0000000000E+00,
                25.0000000000E+00,
                30.0000000000E+00,
                40.0000000000E+00,
                50.0000000000E+00,
                75.0000000000E+00,
                100.0000000000E+00,
                125.0000000000E+00
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

        public static void i1ml1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I1ML1_VALUES returns some values of the I1ML1 function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      I1ML1(x) = I1(x) - L1(x)
            //
            //    I1(x) is the modified Bessel function of the first kind of order 1,
            //    L1(x) is the modified Struve function of order 1.
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 August 2004
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
                0.97575346155386267134E-03,
                0.77609293280609272733E-02,
                0.59302966404545373770E-01,
                0.20395212276737365307E+00,
                0.33839472293667639038E+00,
                0.48787706726961324579E+00,
                0.59018734196576517506E+00,
                0.62604539530312149476E+00,
                0.63209315274909764698E+00,
                0.63410179313235359215E+00,
                0.63417966797578128188E+00,
                0.63439268632392089434E+00,
                0.63501579073257770690E+00,
                0.63559616677359459337E+00,
                0.63591001826697110312E+00,
                0.63622113181751073643E+00,
                0.63636481702133606597E+00,
                0.63650653499619902120E+00,
                0.63655609126300261851E+00,
                0.63657902087183929223E+00
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0156250000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                2.0000000000E+00,
                4.0000000000E+00,
                8.0000000000E+00,
                12.0000000000E+00,
                16.0000000000E+00,
                16.2500000000E+00,
                17.0000000000E+00,
                20.0000000000E+00,
                25.0000000000E+00,
                30.0000000000E+00,
                40.0000000000E+00,
                50.0000000000E+00,
                75.0000000000E+00,
                100.0000000000E+00,
                125.0000000000E+00
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