using System;
using Burkardt.Types;

namespace Burkardt.TriangleNS;

public static partial class LynessRule
{
    public static int[] lyness_suborder(int rule, int suborder_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LYNESS_SUBORDER returns the suborders for a Lyness rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Lyness, Dennis Jespersen,
        //    Moderate Degree Symmetric Quadrature Rules for the Triangle,
        //    Journal of the Institute of Mathematics and its Applications,
        //    Volume 15, Number 1, February 1975, pages 19-32.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Input, int SUBORDER_NUM, the number of suborders 
        //    of the rule.
        //
        //    Output, int LYNESS_SUBORDER[SUBORDER_NUM], the suborders 
        //    of the rule.
        //
    {
        int[] suborder;
        int[] suborder_00 = {1};
        int[] suborder_01 = {3};
        int[] suborder_02 = {1, 3};
        int[] suborder_03 = {1, 3};
        int[] suborder_04 = {1, 3, 3};
        int[] suborder_05 = {3, 3};
        int[] suborder_06 = {1, 3, 6};
        int[] suborder_07 = {3, 3, 3};
        int[] suborder_08 = {1, 3, 3};
        int[] suborder_09 = {1, 3, 3, 3};
        int[] suborder_10 = {3, 3, 6};
        int[] suborder_11 = {1, 3, 3, 3, 6};
        int[] suborder_12 = {1, 3, 3, 6};
        int[] suborder_13 = {1, 3, 3, 6};
        int[] suborder_14 = {1, 3, 3, 3, 6};
        int[] suborder_15 = {1, 3, 3, 3, 6};
        int[] suborder_16 = {3, 3, 3, 3, 3, 6};
        int[] suborder_17 = {1, 3, 3, 3, 6};
        int[] suborder_18 = {1, 3, 3, 3, 3, 6};
        int[] suborder_19 = {1, 3, 3, 3, 3, 3, 6};
        int[] suborder_20 = {3, 3, 3, 3, 3, 6, 6};
        int[] suborder_21 = {1, 3, 3, 3, 3, 3, 6, 6};

        suborder = new int[suborder_num];

        switch (rule)
        {
            case 0:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_00);
                break;
            case 1:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_01);
                break;
            case 2:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_02);
                break;
            case 3:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_03);
                break;
            case 4:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_04);
                break;
            case 5:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_05);
                break;
            case 6:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_06);
                break;
            case 7:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_07);
                break;
            case 8:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_08);
                break;
            case 9:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_09);
                break;
            case 10:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_10);
                break;
            case 11:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_11);
                break;
            case 12:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_12);
                break;
            case 13:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_13);
                break;
            case 14:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_14);
                break;
            case 15:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_15);
                break;
            case 16:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_16);
                break;
            case 17:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_17);
                break;
            case 18:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_18);
                break;
            case 19:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_19);
                break;
            case 20:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_20);
                break;
            case 21:
                suborder = typeMethods.i4vec_copy_new(suborder_num, suborder_21);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LYNESS_SUBORDER - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return null;
        }

        return suborder;
    }

    public static int lyness_suborder_num(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LYNESS_SUBORDER_NUM returns the number of suborders for a Lyness rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Lyness, Dennis Jespersen,
        //    Moderate Degree Symmetric Quadrature Rules for the Triangle,
        //    Journal of the Institute of Mathematics and its Applications,
        //    Volume 15, Number 1, February 1975, pages 19-32.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Output, int LYNESS_SUBORDER_NUM, the number of suborders 
        //    of the rule.
        //
    {
        int suborder_num;

        switch (rule)
        {
            case 0:
            case 1:
                suborder_num = 1;
                break;
            case 2:
            case 3:
                suborder_num = 2;
                break;
            case 4:
                suborder_num = 3;
                break;
            case 5:
                suborder_num = 2;
                break;
            case 6:
            case 7:
            case 8:
                suborder_num = 3;
                break;
            case 9:
                suborder_num = 4;
                break;
            case 10:
                suborder_num = 3;
                break;
            case 11:
                suborder_num = 5;
                break;
            case 12:
            case 13:
                suborder_num = 4;
                break;
            case 14:
            case 15:
                suborder_num = 5;
                break;
            case 16:
                suborder_num = 6;
                break;
            case 17:
                suborder_num = 5;
                break;
            case 18:
                suborder_num = 6;
                break;
            case 19:
            case 20:
                suborder_num = 7;
                break;
            case 21:
                suborder_num = 8;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LYNESS_SUBORDER_NUM - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return 1;
        }

        return suborder_num;
    }

    public static void lyness_subrule(int rule, int suborder_num, ref double[] sub_xyz,
            ref double[] sub_w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LYNESS_SUBRULE returns a compressed Lyness rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Lyness, Dennis Jespersen,
        //    Moderate Degree Symmetric Quadrature Rules for the Triangle,
        //    Journal of the Institute of Mathematics and its Applications,
        //    Volume 15, Number 1, February 1975, pages 19-32.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Input, int SUB_ORDER_NUM, the number of suborders 
        //    of the rule.
        //
        //    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
        //    the barycentric coordinates of the abscissas.
        //
        //    Output, double SUBORDER_W[SUBORDER_NUM], the
        //    suborder weights.
        //
    {
        double[] sub_w_00 =
        {
            1.0000000000E+00
        };
        double[] sub_w_01 =
        {
            1.0000000000E+00
        };
        double[] sub_w_02 =
        {
            0.7500000000E+00,
            0.2500000000E+00
        };
        double[] sub_w_03 =
        {
            -0.5625000000E+00,
            1.5625000000E+00
        };
        double[] sub_w_04 =
        {
            0.45E+00,
            0.15E+00,
            0.40E+00
        };
        double[] sub_w_05 =
        {
            3.298552309659655E-01,
            6.701447690340345E-01
        };
        double[] sub_w_06 =
        {
            0.45E+00,
            -0.05E+00,
            0.60E+00
        };
        double[] sub_w_07 =
        {
            6.16204060378000851E-02,
            0.18592649660480146E+00,
            0.75245309735739840E+00
        };
        double[] sub_w_08 =
        {
            0.22500000000000000E+00,
            0.37781754163448150E+00,
            0.39718245836551852E+00
        };
        double[] sub_w_09 =
        {
            0.25312500000000000E+00,
            0.03333333333333333E+00,
            0.21333333333333333E+00,
            0.50020833333333333E+00
        };
        double[] sub_w_10 =
        {
            3.503588271790222E-01,
            1.525347191106164E-01,
            4.971064537103575E-01
        };
        double[] sub_w_11 =
        {
            -0.57857142857142863E+00,
            -5.95238095238095205E-02,
            0.16190476190476191E+00,
            1.2190476190476192E+00,
            0.25714285714285712E+00
        };
        double[] sub_w_12 =
        {
            1.527089667883523E-01,
            2.944076042366762E-01,
            3.887052878418766E-01,
            1.641781411330949E-01
        };
        double[] sub_w_13 =
        {
            -1.495700444677495E-01,
            5.268457722996828E-01,
            1.600417068265167E-01,
            4.626825653415500E-01
        };
        double[] sub_w_14 =
        {
            1.763126156005252E-01,
            1.210901532768310E-02,
            3.499561757697094E-01,
            3.195119754425220E-01,
            1.421102178595603E-01
        };
        double[] sub_w_15 =
        {
            1.443156076777862E-01,
            2.852749028018549E-01,
            9.737549286959440E-02,
            3.096521116041552E-01,
            1.633818850466092E-01
        };
        double[] sub_w_16 =
        {
            1.207273935292775E-02,
            -8.491579879151455E-01,
            1.042367468891334E+00,
            1.947229791412260E-01,
            4.511852767201322E-01,
            1.488095238095238E-01
        };
        double[] sub_w_17 =
        {
            -2.834183851113958E-01,
            2.097208857979572E-01,
            5.127273801480265E-02,
            6.564896469913506E-01,
            3.659351143072855E-01
        };
        double[] sub_w_18 =
        {
            9.713579628279610E-02,
            9.400410068141950E-02,
            2.334826230143263E-01,
            2.389432167816273E-01,
            7.673302697609430E-02,
            2.597012362637364E-01
        };
        double[] sub_w_19 =
        {
            1.133624844599192E-01,
            1.062573789846380E-03,
            4.803411513859279E-02,
            2.524243006337300E-01,
            7.819254371487040E-02,
            2.472227459993048E-01,
            2.597012362637364E-01
        };
        double[] sub_w_20 =
        {
            4.097919300803106E-02,
            1.085536215102866E-01,
            2.781018986881812E-03,
            1.779689321422668E-01,
            2.314486047444677E-01,
            3.140226717732234E-01,
            1.242459578348437E-01
        };
        double[] sub_w_21 =
        {
            8.797730116222190E-02,
            2.623293466120857E-02,
            1.142447159818060E-01,
            5.656634416839376E-02,
            2.164790926342230E-01,
            2.079874161166116E-01,
            4.417430269980344E-02,
            2.463378925757316E-01
        };
        double[] sub_xyz_00 =
        {
            0.3333333333E+00, 0.3333333333E+00, 0.3333333334E+00
        };
        double[] sub_xyz_01 =
        {
            0.0000000000E+00, 0.5000000000E+00, 0.5000000000E+00
        };
        double[] sub_xyz_02 =
        {
            0.3333333333E+00, 0.3333333333E+00, 0.3333333334E+00,
            1.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00
        };
        double[] sub_xyz_03 =
        {
            0.3333333333E+00, 0.3333333333E+00, 0.3333333334E+00,
            0.6000000000E+00, 0.2000000000E+00, 0.2000000000E+00
        };
        double[] sub_xyz_04 =
        {
            0.33333333333333333E+00, 0.33333333333333333E+00, 0.33333333333333333E+00,
            1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00,
            0.00000000000000000E+00, 0.50000000000000000E+00, 0.50000000000000000E+00
        };
        double[] sub_xyz_05 =
        {
            8.168475729804585E-01, 9.157621350977073E-02, 9.15762135097707569E-02,
            1.081030181680702E-01, 4.459484909159649E-01, 0.44594849091596489E+00
        };
        double[] sub_xyz_06 =
        {
            0.33333333333333333E+00, 0.33333333333333333E+00, 0.33333333333333333E+00,
            1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00,
            0.00000000000000000E+00, 0.78867513459481281E+00, 0.21132486540518719E+00
        };
        double[] sub_xyz_07 =
        {
            1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00,
            0.00000000000000000E+00, 0.50000000000000000E+00, 0.50000000000000000E+00,
            0.62283903060710999E+00, 0.18858048469644506E+00, 0.18858048469644506E+00
        };
        double[] sub_xyz_08 =
        {
            0.33333333333333333E+00, 0.33333333333333333E+00, 0.33333333333333333E+00,
            0.79742698535308720E+00, 0.10128650732345633E+00, 0.10128650732345633E+00,
            5.97158717897698088E-02, 0.47014206410511505E+00, 0.47014206410511505E+00
        };
        double[] sub_xyz_09 =
        {
            0.33333333333333333E+00, 0.33333333333333333E+00, 0.33333333333333333E+00,
            1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00,
            0.00000000000000000E+00, 0.50000000000000000E+00, 0.50000000000000000E+00,
            0.71428571428571430E-00, 0.14285714285714285E+00, 0.14285714285714285E+00
        };
        double[] sub_xyz_10 =
        {
            5.014265096581342E-01, 2.492867451709329E-01, 0.24928674517093291E+00,
            8.738219710169965E-01, 6.308901449150177E-02, 6.30890144915016854E-02,
            6.365024991213939E-01, 5.314504984483216E-02, 0.31035245103377396E+00
        };
        double[] sub_xyz_11 =
        {
            0.33333333333333333E+00, 0.33333333333333333E+00, 0.33333333333333333E+00,
            1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00,
            0.00000000000000000E+00, 0.50000000000000000E+00, 0.50000000000000000E+00,
            0.50000000000000000E+00, 0.25000000000000000E+00, 0.25000000000000000E+00,
            0.00000000000000000E+00, 0.90824829046386302E+00, 9.17517095361370244E-02
        };
        double[] sub_xyz_12 =
        {
            3.333333333333333E-01, 3.333333333333333E-01, 0.33333333333333333E+00,
            5.233837209269747E-02, 4.738308139536513E-01, 0.47383081395365129E+00,
            6.557646607383649E-01, 1.721176696308175E-01, 0.17211766963081757E+00,
            0.000000000000000E+00, 8.653073540834571E-01, 0.13469264591654295E+00
        };
        double[] sub_xyz_13 =
        {
            3.333333333333333E-01, 3.333333333333333E-01, 0.33333333333333333E+00,
            4.793080678419067E-01, 2.603459660790466E-01, 0.26034596607904670E+00,
            8.697397941955675E-01, 6.513010290221623E-02, 6.51301029022163108E-02,
            6.384441885698096E-01, 4.869031542531756E-02, 0.31286549600487290E+00
        };
        double[] sub_xyz_14 =
        {
            3.333333333333333E-01, 3.333333333333333E-01, 0.33333333333333333E+00,
            1.000000000000000E+00, 0.000000000000000E+00, 0.00000000000000000E+00,
            6.901278795524791E-01, 1.549360602237604E-01, 0.15493606022376047E+00,
            6.169850771237593E-02, 4.691507461438120E-01, 0.46915074614381203E+00,
            0.000000000000000E+00, 8.392991722729236E-01, 0.16070082772707639E+00
        };
        double[] sub_xyz_15 =
        {
            3.333333333333333E-01, 3.333333333333333E-01, 0.33333333333333333E+00,
            8.141482341455413E-02, 4.592925882927229E-01, 0.45929258829272301E+00,
            8.989055433659379E-01, 5.054722831703103E-02, 5.05472283170311024E-02,
            6.588613844964797E-01, 1.705693077517601E-01, 0.17056930775176021E+00,
            8.394777409957211E-03, 7.284923929554041E-01, 0.26311282963463867E+00
        };
        double[] sub_xyz_16 =
        {
            1.000000000000000E+00, 0.000000000000000E+00, 0.00000000000000000E+00,
            0.000000000000000E+00, 0.500000000000000E+00, 0.50000000000000000E+00,
            8.637211648883667E-03, 4.956813941755582E-01, 0.49568139417555818E+00,
            8.193444849714693E-01, 9.032775751426533E-02, 9.03277575142653888E-02,
            5.316905005853895E-01, 2.341547497073052E-01, 0.23415474970730532E+00,
            0.000000000000000E-01, 7.236067977499790E-01, 0.27639320225002095E+00
        };
        double[] sub_xyz_17 =
        {
            3.333333333333333E-01, 3.333333333333333E-01, 0.33333333333333333E+00,
            4.666912123569507E-02, 4.766654393821525E-01, 0.47666543938215239E+00,
            9.324563118910393E-01, 3.377184405448033E-02, 3.37718440544803877E-02,
            4.593042216691921E-01, 2.703478891654040E-01, 0.27034788916540398E+00,
            5.146433548666149E-02, 7.458294907672514E-01, 0.20270617374608713E+00
        };
        double[] sub_xyz_18 =
        {
            3.333333333333333E-01, 3.333333333333333E-01, 0.33333333333333333E+00,
            2.063496160252593E-02, 4.896825191987370E-01, 0.48968251919873701E+00,
            1.258208170141290E-01, 4.370895914929355E-01, 0.43708959149293541E+00,
            6.235929287619356E-01, 1.882035356190322E-01, 0.18820353561903219E+00,
            9.105409732110941E-01, 4.472951339445297E-02, 4.47295133944529688E-02,
            3.683841205473626E-02, 7.411985987844980E-01, 0.22196298916076573E+00
        };
        double[] sub_xyz_19 =
        {
            3.333333333333333E-01, 3.333333333333333E-01, 0.33333333333333333E+00,
            1.000000000000000E+00, 0.000000000000000E+00, 0.00000000000000000E+00,
            0.000000000000000E+00, 0.500000000000000E+00, 0.50000000000000000E+00,
            1.004413236259677E-01, 4.497793381870162E-01, 0.44977933818701610E+00,
            9.061051136018193E-01, 4.694744319909033E-02, 4.69474431990903329E-02,
            6.162561745251021E-01, 1.918719127374489E-01, 0.19187191273744902E+00,
            3.683841205473626E-02, 7.411985987844980E-01, 0.22196298916076573E+00
        };
        double[] sub_xyz_20 =
        {
            9.352701037774565E-01, 3.236494811127173E-02, 3.23649481112718157E-02,
            7.612981754348137E-01, 1.193509122825931E-01, 0.11935091228259319E+00,
            -6.922209654151433E-02, 5.346110482707572E-01, 0.53461104827075701E+00,
            5.933801991374367E-01, 2.033099004312816E-01, 0.20330990043128172E+00,
            2.020613940682885E-01, 3.989693029658558E-01, 0.39896930296585570E+00,
            5.017813831049474E-02, 5.932012134282132E-01, 0.35662064826129203E+00,
            2.102201653616613E-02, 8.074890031597923E-01, 0.17148898030404158E+00
        };
        double[] sub_xyz_21 =
        {
            3.333333333333333E-01, 3.333333333333333E-01, 0.33333333333333333E+00,
            9.480217181434233E-01, 2.598914092828833E-02, 2.59891409282883845E-02,
            8.114249947041546E-01, 9.428750264792270E-02, 9.42875026479226691E-02,
            1.072644996557060E-02, 4.946367750172147E-01, 0.49463677501721470E+00,
            5.853132347709715E-01, 2.073433826145142E-01, 0.20734338261451427E+00,
            1.221843885990187E-01, 4.389078057004907E-01, 0.43890780570049059E+00,
            0.000000000000000E+00, 8.588702812826364E-01, 0.14112971871736357E+00,
            4.484167758913055E-02, 6.779376548825902E-01, 0.27722066752827923E+00
        };

        switch (rule)
        {
            case 0:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_00, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_00, ref sub_w);
                break;
            case 1:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_01, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_01, ref sub_w);
                break;
            case 2:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_02, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_02, ref sub_w);
                break;
            case 3:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_03, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_03, ref sub_w);
                break;
            case 4:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_04, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_04, ref sub_w);
                break;
            case 5:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_05, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_05, ref sub_w);
                break;
            case 6:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_06, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_06, ref sub_w);
                break;
            case 7:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_07, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_07, ref sub_w);
                break;
            case 8:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_08, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_08, ref sub_w);
                break;
            case 9:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_09, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_09, ref sub_w);
                break;
            case 10:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_10, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_10, ref sub_w);
                break;
            case 11:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_11, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_11, ref sub_w);
                break;
            case 12:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_12, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_12, ref sub_w);
                break;
            case 13:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_13, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_13, ref sub_w);
                break;
            case 14:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_14, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_14, ref sub_w);
                break;
            case 15:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_15, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_15, ref sub_w);
                break;
            case 16:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_16, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_16, ref sub_w);
                break;
            case 17:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_17, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_17, ref sub_w);
                break;
            case 18:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_18, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_18, ref sub_w);
                break;
            case 19:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_19, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_19, ref sub_w);
                break;
            case 20:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_20, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_20, ref sub_w);
                break;
            case 21:
                typeMethods.r8mat_copy(3, suborder_num, sub_xyz_21, ref sub_xyz);
                typeMethods.r8vec_copy(suborder_num, sub_w_21, ref sub_w);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LYNESS_SUBRULE - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                break;
        }

    }
}