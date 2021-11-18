namespace Burkardt.TetrahedronNS;

public static class Integrand
{
    public static double[] tetrahedron_integrand_01(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_INTEGRAND_01 evaluates 1 integrand function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P[3*P_NUM], the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP[F_NUM*P_NUM], the integrand values.
        //
    {
        int j;

        double[] fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = 1.0;
        }

        return fp;
    }

    public static double[] tetrahedron_integrand_02(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_INTEGRAND_02 evaluates 3 integrand functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P[3*P_NUM], the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP[F_NUM*P_NUM], the integrand values.
        //
    {
        int j;

        double[] fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = p[0 + j * 3];
            fp[1 + j * f_num] = p[1 + j * 3];
            fp[2 + j * f_num] = p[2 + j * 3];
        }

        return fp;
    }

    public static double[] tetrahedron_integrand_03(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_INTEGRAND_03 evaluates 6 integrand functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P[3*P_NUM], the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP[F_NUM*P_NUM], the integrand values.
        //
    {
        int j;

        double[] fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = p[0 + j * 3] * p[0 + j * 3];
            fp[1 + j * f_num] = p[0 + j * 3] * p[1 + j * 3];
            fp[2 + j * f_num] = p[0 + j * 3] * p[2 + j * 3];
            fp[3 + j * f_num] = p[1 + j * 3] * p[1 + j * 3];
            fp[4 + j * f_num] = p[1 + j * 3] * p[2 + j * 3];
            fp[5 + j * f_num] = p[2 + j * 3] * p[2 + j * 3];
        }

        return fp;
    }

    public static double[] tetrahedron_integrand_04(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_INTEGRAND_04 evaluates 10 integrand functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P[3*P_NUM], the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP[F_NUM*P_NUM], the integrand values.
        //
    {
        int j;

        double[] fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = p[0 + j * 3] * p[0 + j * 3] * p[0 + j * 3];
            fp[1 + j * f_num] = p[0 + j * 3] * p[0 + j * 3] * p[1 + j * 3];
            fp[2 + j * f_num] = p[0 + j * 3] * p[0 + j * 3] * p[2 + j * 3];
            fp[3 + j * f_num] = p[0 + j * 3] * p[1 + j * 3] * p[1 + j * 3];
            fp[4 + j * f_num] = p[0 + j * 3] * p[1 + j * 3] * p[2 + j * 3];
            fp[5 + j * f_num] = p[0 + j * 3] * p[2 + j * 3] * p[2 + j * 3];
            fp[6 + j * f_num] = p[1 + j * 3] * p[1 + j * 3] * p[1 + j * 3];
            fp[7 + j * f_num] = p[1 + j * 3] * p[1 + j * 3] * p[2 + j * 3];
            fp[8 + j * f_num] = p[1 + j * 3] * p[2 + j * 3] * p[2 + j * 3];
            fp[9 + j * f_num] = p[2 + j * 3] * p[2 + j * 3] * p[2 + j * 3];
        }

        return fp;
    }

    public static double[] tetrahedron_integrand_05(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_INTEGRAND_05 evaluates 15 integrand functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P[3*P_NUM], the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP[F_NUM*P_NUM], the integrand values.
        //
    {
        int j;

        double[] fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = p[0 + j * 3] * p[0 + j * 3] * p[0 + j * 3] * p[0 + j * 3];
            fp[1 + j * f_num] = p[0 + j * 3] * p[0 + j * 3] * p[0 + j * 3] * p[1 + j * 3];
            fp[2 + j * f_num] = p[0 + j * 3] * p[0 + j * 3] * p[0 + j * 3] * p[2 + j * 3];
            fp[3 + j * f_num] = p[0 + j * 3] * p[0 + j * 3] * p[1 + j * 3] * p[1 + j * 3];
            fp[4 + j * f_num] = p[0 + j * 3] * p[0 + j * 3] * p[1 + j * 3] * p[2 + j * 3];
            fp[5 + j * f_num] = p[0 + j * 3] * p[0 + j * 3] * p[2 + j * 3] * p[2 + j * 3];
            fp[6 + j * f_num] = p[0 + j * 3] * p[1 + j * 3] * p[1 + j * 3] * p[1 + j * 3];
            fp[7 + j * f_num] = p[0 + j * 3] * p[1 + j * 3] * p[1 + j * 3] * p[2 + j * 3];
            fp[8 + j * f_num] = p[0 + j * 3] * p[1 + j * 3] * p[2 + j * 3] * p[2 + j * 3];
            fp[9 + j * f_num] = p[0 + j * 3] * p[2 + j * 3] * p[2 + j * 3] * p[2 + j * 3];
            fp[10 + j * f_num] = p[1 + j * 3] * p[1 + j * 3] * p[1 + j * 3] * p[1 + j * 3];
            fp[11 + j * f_num] = p[1 + j * 3] * p[1 + j * 3] * p[1 + j * 3] * p[2 + j * 3];
            fp[12 + j * f_num] = p[1 + j * 3] * p[1 + j * 3] * p[2 + j * 3] * p[2 + j * 3];
            fp[13 + j * f_num] = p[1 + j * 3] * p[2 + j * 3] * p[2 + j * 3] * p[2 + j * 3];
            fp[14 + j * f_num] = p[2 + j * 3] * p[2 + j * 3] * p[2 + j * 3] * p[2 + j * 3];
        }

        return fp;
    }
}