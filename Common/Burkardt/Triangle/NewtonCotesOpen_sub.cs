using System;

namespace Burkardt.TriangleNS;

public static partial class NewtonCotesOpen
{
    public static int[] triangle_nco_suborder(int rule, int suborder_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBORDER returns the suborders for an NCO rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, int TRIANGLE_NCO_SUBORDER[SUBORDER_NUM],
        //    the suborders of the rule.
        //
    {
        int[] suborder = new int[suborder_num];

        switch (rule)
        {
            case 1:
                suborder[0] = 1;
                break;
            case 2:
                suborder[0] = 3;
                break;
            case 3:
                suborder[0] = 3;
                suborder[1] = 3;
                break;
            case 4:
                suborder[0] = 3;
                suborder[1] = 6;
                suborder[2] = 1;
                break;
            case 5:
                suborder[0] = 3;
                suborder[1] = 6;
                suborder[2] = 3;
                suborder[3] = 3;
                break;
            case 6:
                suborder[0] = 3;
                suborder[1] = 6;
                suborder[2] = 6;
                suborder[3] = 3;
                suborder[4] = 3;
                break;
            case 7:
                suborder[0] = 3;
                suborder[1] = 6;
                suborder[2] = 6;
                suborder[3] = 3;
                suborder[4] = 3;
                suborder[5] = 6;
                suborder[6] = 1;
                break;
            case 8:
                suborder[0] = 3;
                suborder[1] = 6;
                suborder[2] = 6;
                suborder[3] = 3;
                suborder[4] = 6;
                suborder[5] = 6;
                suborder[6] = 3;
                suborder[7] = 3;
                break;
            case 9:
                suborder[0] = 3;
                suborder[1] = 6;
                suborder[2] = 6;
                suborder[3] = 3;
                suborder[4] = 6;
                suborder[5] = 6;
                suborder[6] = 3;
                suborder[7] = 6;
                suborder[8] = 3;
                suborder[9] = 3;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_NCO_SUBORDER - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                break;
        }

        return suborder;
    }

    public static int triangle_nco_suborder_num(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBORDER_NUM returns the number of suborders for an NCO rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Output, int TRIANGLE_NCO_SUBORDER_NUM, the number of suborders
        //    of the rule.
        //
    {
        int suborder_num;

        switch (rule)
        {
            case 1:
            case 2:
                suborder_num = 1;
                break;
            case 3:
                suborder_num = 2;
                break;
            case 4:
                suborder_num = 3;
                break;
            case 5:
                suborder_num = 4;
                break;
            case 6:
                suborder_num = 5;
                break;
            case 7:
                suborder_num = 7;
                break;
            case 8:
                suborder_num = 8;
                break;
            case 9:
                suborder_num = 10;
                break;
            default:
                suborder_num = -1;
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_NCO_SUBORDER_NUM - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                break;
        }

        return suborder_num;
    }

    public static void triangle_nco_subrule(int rule, int suborder_num, ref double[] suborder_xyz,
            ref double[] suborder_w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBRULE returns a compressed NCO rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
        //    the barycentric coordinates of the abscissas.
        //
        //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
        //
    {
        int s;
        int suborder_w_d = 0;
        int suborder_xyz_d = 0;

        int[] suborder_xyz_n = new int[3 * suborder_num];
        int[] suborder_w_n = new int[suborder_num];

        switch (rule)
        {
            case 1:
                triangle_nco_subrule_01(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
                break;
            case 2:
                triangle_nco_subrule_02(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
                break;
            case 3:
                triangle_nco_subrule_03(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
                break;
            case 4:
                triangle_nco_subrule_04(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
                break;
            case 5:
                triangle_nco_subrule_05(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
                break;
            case 6:
                triangle_nco_subrule_06(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
                break;
            case 7:
                triangle_nco_subrule_07(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
                break;
            case 8:
                triangle_nco_subrule_08(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
                break;
            case 9:
                triangle_nco_subrule_09(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_NCO_SUBRULE - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return;
        }

        for (s = 0; s < suborder_num; s++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                suborder_xyz[i + s * 3] =
                    (1 + suborder_xyz_n[i + s * 3])
                    / (double) (3 + suborder_xyz_d);
            }
        }

        for (s = 0; s < suborder_num; s++)
        {
            suborder_w[s] = suborder_w_n[s] / (double) suborder_w_d;
        }
    }

    public static void triangle_nco_subrule_01(int suborder_num, ref int[] suborder_xyz_n,
            ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBRULE_01 returns a compressed NCO rule 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
        //    the numerators of the barycentric coordinates of the abscissas.
        //
        //    Output, int[] SUBORDER_XYZ_D,
        //    the denominator of the barycentric coordinates of the abscissas.
        //
        //    Output, int SUBORDER_W_N[SUBORDER_NUM],
        //    the numerator of the suborder weights.
        //
        //    Output, int SUBORDER_W_D,
        //    the denominator of the suborder weights.
        //
    {
        int s;
        int[] suborder_xyz_n_01 =
        {
            0, 0, 0
        };
        const int suborder_xyz_d_01 = 0;
        int[] suborder_w_n_01 = {1};
        const int suborder_w_d_01 = 1;

        for (s = 0; s < suborder_num; s++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                suborder_xyz_n[i + s * 3] = suborder_xyz_n_01[i + s * 3];
            }
        }

        suborder_xyz_d = suborder_xyz_d_01;

        for (s = 0; s < suborder_num; s++)
        {
            suborder_w_n[s] = suborder_w_n_01[s];
        }

        suborder_w_d = suborder_w_d_01;

    }

    public static void triangle_nco_subrule_02(int suborder_num, ref int[] suborder_xyz_n,
            ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBRULE_02 returns a compressed NCO rule 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
        //    the numerators of the barycentric coordinates of the abscissas.
        //
        //    Output, int[] SUBORDER_XYZ_D,
        //    the denominator of the barycentric coordinates of the abscissas.
        //
        //    Output, int SUBORDER_W_N[SUBORDER_NUM],
        //    the numerator of the suborder weights.
        //
        //    Output, int SUBORDER_W_D,
        //    the denominator of the suborder weights.
        //
    {
        int s;
        int[] suborder_xyz_n_02 =
        {
            1, 0, 0
        };
        const int suborder_xyz_d_02 = 1;
        int[] suborder_w_n_02 = {1};
        const int suborder_w_d_02 = 3;

        for (s = 0; s < suborder_num; s++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                suborder_xyz_n[i + s * 3] = suborder_xyz_n_02[i + s * 3];
            }
        }

        suborder_xyz_d = suborder_xyz_d_02;

        for (s = 0; s < suborder_num; s++)
        {
            suborder_w_n[s] = suborder_w_n_02[s];
        }

        suborder_w_d = suborder_w_d_02;

    }

    public static void triangle_nco_subrule_03(int suborder_num, ref int[] suborder_xyz_n,
            ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBRULE_03 returns a compressed NCO rule 3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
        //    the numerators of the barycentric coordinates of the abscissas.
        //
        //    Output, int[] SUBORDER_XYZ_D,
        //    the denominator of the barycentric coordinates of the abscissas.
        //
        //    Output, int SUBORDER_W_N[SUBORDER_NUM],
        //    the numerator of the suborder weights.
        //
        //    Output, int SUBORDER_W_D,
        //    the denominator of the suborder weights.
        //
    {
        int s;
        int[] suborder_xyz_n_03 =
        {
            2, 0, 0,
            1, 1, 0
        };
        const int suborder_xyz_d_03 = 2;
        int[] suborder_w_n_03 = {7, -3};
        const int suborder_w_d_03 = 12;

        for (s = 0; s < suborder_num; s++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                suborder_xyz_n[i + s * 3] = suborder_xyz_n_03[i + s * 3];
            }
        }

        suborder_xyz_d = suborder_xyz_d_03;

        for (s = 0; s < suborder_num; s++)
        {
            suborder_w_n[s] = suborder_w_n_03[s];
        }

        suborder_w_d = suborder_w_d_03;

    }

    public static void triangle_nco_subrule_04(int suborder_num, ref int[] suborder_xyz_n,
            ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBRULE_04 returns a compressed NCO rule 4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
        //    the numerators of the barycentric coordinates of the abscissas.
        //
        //    Output, int[] SUBORDER_XYZ_D,
        //    the denominator of the barycentric coordinates of the abscissas.
        //
        //    Output, int SUBORDER_W_N[SUBORDER_NUM],
        //    the numerator of the suborder weights.
        //
        //    Output, int SUBORDER_W_D,
        //    the denominator of the suborder weights.
        //
    {
        int s;
        int[] suborder_xyz_n_04 =
        {
            3, 0, 0,
            2, 1, 0,
            1, 1, 1
        };
        const int suborder_xyz_d_04 = 3;
        int[] suborder_w_n_04 = {8, 3, -12};
        const int suborder_w_d_04 = 30;

        for (s = 0; s < suborder_num; s++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                suborder_xyz_n[i + s * 3] = suborder_xyz_n_04[i + s * 3];
            }
        }

        suborder_xyz_d = suborder_xyz_d_04;

        for (s = 0; s < suborder_num; s++)
        {
            suborder_w_n[s] = suborder_w_n_04[s];
        }

        suborder_w_d = suborder_w_d_04;

    }

    public static void triangle_nco_subrule_05(int suborder_num, ref int[] suborder_xyz_n,
            ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBRULE_05 returns a compressed NCO rule 5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
        //    the numerators of the barycentric coordinates of the abscissas.
        //
        //    Output, int[] SUBORDER_XYZ_D,
        //    the denominator of the barycentric coordinates of the abscissas.
        //
        //    Output, int SUBORDER_W_N[SUBORDER_NUM],
        //    the numerator of the suborder weights.
        //
        //    Output, int SUBORDER_W_D,
        //    the denominator of the suborder weights.
        //
    {
        int s;
        int[] suborder_xyz_n_05 =
        {
            4, 0, 0,
            3, 1, 0,
            2, 2, 0,
            2, 1, 1
        };
        const int suborder_xyz_d_05 = 4;
        int[] suborder_w_n_05 = {307, -316, 629, -64};
        const int suborder_w_d_05 = 720;

        for (s = 0; s < suborder_num; s++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                suborder_xyz_n[i + s * 3] = suborder_xyz_n_05[i + s * 3];
            }
        }

        suborder_xyz_d = suborder_xyz_d_05;

        for (s = 0; s < suborder_num; s++)
        {
            suborder_w_n[s] = suborder_w_n_05[s];
        }

        suborder_w_d = suborder_w_d_05;

    }

    public static void triangle_nco_subrule_06(int suborder_num, ref int[] suborder_xyz_n,
            ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBRULE_06 returns a compressed NCO rule 6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
        //    the numerators of the barycentric coordinates of the abscissas.
        //
        //    Output, int[] SUBORDER_XYZ_D,
        //    the denominator of the barycentric coordinates of the abscissas.
        //
        //    Output, int SUBORDER_W_N[SUBORDER_NUM],
        //    the numerator of the suborder weights.
        //
        //    Output, int SUBORDER_W_D,
        //    the denominator of the suborder weights.
        //
    {
        int s;
        int[] suborder_xyz_n_06 =
        {
            5, 0, 0,
            4, 1, 0,
            3, 2, 0,
            3, 1, 1,
            2, 2, 1
        };
        const int suborder_xyz_d_06 = 5;
        int[] suborder_w_n_06 = {71, -13, 57, -167, 113};
        const int suborder_w_d_06 = 315;

        for (s = 0; s < suborder_num; s++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                suborder_xyz_n[i + s * 3] = suborder_xyz_n_06[i + s * 3];
            }
        }

        suborder_xyz_d = suborder_xyz_d_06;

        for (s = 0; s < suborder_num; s++)
        {
            suborder_w_n[s] = suborder_w_n_06[s];
        }

        suborder_w_d = suborder_w_d_06;

    }

    public static void triangle_nco_subrule_07(int suborder_num, ref int[] suborder_xyz_n,
            ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBRULE_07 returns a compressed NCO rule 7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
        //    the numerators of the barycentric coordinates of the abscissas.
        //
        //    Output, int[] SUBORDER_XYZ_D,
        //    the denominator of the barycentric coordinates of the abscissas.
        //
        //    Output, int SUBORDER_W_N[SUBORDER_NUM],
        //    the numerator of the suborder weights.
        //
        //    Output, int SUBORDER_W_D,
        //    the denominator of the suborder weights.
        //
    {
        int s;
        int[] suborder_xyz_n_07 =
        {
            6, 0, 0,
            5, 1, 0,
            4, 2, 0,
            4, 1, 1,
            3, 3, 0,
            3, 2, 1,
            2, 2, 2
        };
        const int suborder_xyz_d_07 = 6;
        int[] suborder_w_n_07 = {767, -1257, 2901, 387, -3035, -915, 3509};
        const int suborder_w_d_07 = 2240;

        for (s = 0; s < suborder_num; s++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                suborder_xyz_n[i + s * 3] = suborder_xyz_n_07[i + s * 3];
            }
        }

        suborder_xyz_d = suborder_xyz_d_07;

        for (s = 0; s < suborder_num; s++)
        {
            suborder_w_n[s] = suborder_w_n_07[s];
        }

        suborder_w_d = suborder_w_d_07;

    }

    public static void triangle_nco_subrule_08(int suborder_num, ref int[] suborder_xyz_n,
            ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBRULE_08 returns a compressed NCO rule 8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
        //    the numerators of the barycentric coordinates of the abscissas.
        //
        //    Output, int[] SUBORDER_XYZ_D,
        //    the denominator of the barycentric coordinates of the abscissas.
        //
        //    Output, int SUBORDER_W_N[SUBORDER_NUM],
        //    the numerator of the suborder weights.
        //
        //    Output, int SUBORDER_W_D,
        //    the denominator of the suborder weights.
        //
    {
        int s;
        int[] suborder_xyz_n_08 =
        {
            7, 0, 0,
            6, 1, 0,
            5, 2, 0,
            5, 1, 1,
            4, 3, 0,
            4, 2, 1,
            3, 3, 1,
            3, 2, 2
        };
        const int suborder_xyz_d_08 = 7;
        int[] suborder_w_n_08 = {898, -662, 1573, -2522, -191, 2989, -5726, 1444};
        const int suborder_w_d_08 = 4536;

        for (s = 0; s < suborder_num; s++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                suborder_xyz_n[i + s * 3] = suborder_xyz_n_08[i + s * 3];
            }
        }

        suborder_xyz_d = suborder_xyz_d_08;

        for (s = 0; s < suborder_num; s++)
        {
            suborder_w_n[s] = suborder_w_n_08[s];
        }

        suborder_w_d = suborder_w_d_08;

    }

    public static void triangle_nco_subrule_09(int suborder_num, ref int[] suborder_xyz_n,
            ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCO_SUBRULE_09 returns a compressed NCO rule 9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int SUBORDER_NUM, the number of suborders of the rule.
        //
        //    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
        //    the numerators of the barycentric coordinates of the abscissas.
        //
        //    Output, int[] SUBORDER_XYZ_D,
        //    the denominator of the barycentric coordinates of the abscissas.
        //
        //    Output, int SUBORDER_W_N[SUBORDER_NUM],
        //    the numerator of the suborder weights.
        //
        //    Output, int SUBORDER_W_D,
        //    the denominator of the suborder weights.
        //
    {
        int s;
        int[] suborder_xyz_n_09 =
        {
            8, 0, 0,
            7, 1, 0,
            6, 2, 0,
            6, 1, 1,
            5, 3, 0,
            5, 2, 1,
            4, 4, 0,
            4, 3, 1,
            4, 2, 2,
            3, 3, 2
        };
        const int suborder_xyz_d_09 = 8;
        int[] suborder_w_n_09 =
        {
            1051445, -2366706, 6493915, 1818134, -9986439, -3757007, 12368047,
            478257, 10685542, -6437608
        };
        const int suborder_w_d_09 = 3628800;

        for (s = 0; s < suborder_num; s++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                suborder_xyz_n[i + s * 3] = suborder_xyz_n_09[i + s * 3];
            }
        }

        suborder_xyz_d = suborder_xyz_d_09;

        for (s = 0; s < suborder_num; s++)
        {
            suborder_w_n[s] = suborder_w_n_09[s];
        }

        suborder_w_d = suborder_w_d_09;

    }

}