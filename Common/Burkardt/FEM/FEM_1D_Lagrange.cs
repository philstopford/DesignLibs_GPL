using System;
using Burkardt.Quadrature;

namespace Burkardt.FEM;

public static class FEM_1D_Lagrange
{
    public static void fem1d_lagrange_stiffness(int x_num, double[] x, int q_num,
            Func<double, double> f, ref double[] a, ref double[] m, ref double[] b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM1D_LAGRANGE_STIFFNESS evaluates the Lagrange polynomial stiffness matrix.
        //
        //  Discussion:
        //
        //    The finite element method is to be applied over a given interval that
        //    has been meshed with X_NUM points X.
        //
        //    The finite element basis functions are to be the X_NUM Lagrange
        //    basis polynomials L(i)(X), such that
        //      L(i)(X(j)) = delta(i,j).
        //
        //    The following items are computed:
        //    * A, the stiffness matrix, with A(I,J) = integral L'(i)(x) L'(j)(x)
        //    * M, the mass matrix, with M(I,J) = integral L(i)(x) L(j)(x)
        //    * B, the load matrix, with B(I) = integral L(i)(x) F(x)
        //
        //    The integrals are approximated by quadrature.
        //
        //    Boundary conditions are not handled by this routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of nodes.
        //
        //    Input, double X[X_NUM], the coordinates of the nodes.
        //
        //    Input, int Q_NUM, the number of quadrature points to use.
        //
        //    Input, double F ( double X ), the right hand side function.
        //
        //    Output, double A[X_NUM*X_NUM], the stiffness matrix.
        //
        //    Output, double M[X_NUM*X_NUM], the mass matrix.
        //
        //    Output, double B[X_NUM], the right hand side vector.
        //
    {
        //
        //  Get the quadrature rule for [-1,+1].
        //
        double[] q_x = new double[q_num];
        double[] q_w = new double[q_num];

        LegendreQuadrature.legendre_set(q_num, ref q_x, ref q_w);
        //
        //  Adjust the quadrature rule to the interval [ x(1), x(x_num) }
        //
        for (int q_i = 0; q_i < q_num; q_i++)
        {
            q_x[q_i] = ((1.0 - q_x[q_i]) * x[0]
                        + (1.0 + q_x[q_i]) * x[x_num - 1])
                       / 2.0;

            q_w[q_i] = q_w[q_i] * (x[x_num - 1] - x[0]) / 2.0;
        }

        //
        //  Evaluate all the Lagrange basis polynomials at all the quadrature points.
        //
        double[] l = lagrange_value(x_num, x, q_num, q_x);
        //
        //  Evaluate all the Lagrange basis polynomial derivatives at all 
        //  the quadrature points.
        //
        double[] lp = lagrange_derivative(x_num, x, q_num, q_x);
        //
        //  Assemble the matrix and right hand side.
        //
        for (int x_j = 0; x_j < x_num; x_j++)
        {
            for (int x_i = 0; x_i < x_num; x_i++)
            {
                a[x_i + x_j * x_num] = 0.0;
                m[x_i + x_j * x_num] = 0.0;
            }

            b[x_j] = 0.0;
        }

        for (int x_i = 0; x_i < x_num; x_i++)
        {
            for (int q_i = 0; q_i < q_num; q_i++)
            {
                double li = l[q_i + x_i * q_num];
                double lpi = lp[q_i + x_i * q_num];
                for (int x_j = 0; x_j < x_num; x_j++)
                {
                    double lj = l[q_i + x_j * q_num];
                    double lpj = lp[q_i + x_j * q_num];
                    a[x_i + x_j * x_num] += q_w[q_i] * lpi * lpj;
                    m[x_i + x_j * x_num] += q_w[q_i] * li * lj;
                }

                b[x_i] += q_w[q_i] * li * f(q_x[q_i]);
            }
        }
    }

    public static double[] lagrange_derivative(int nd, double[] xd, int ni, double[] xi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_DERIVATIVE evaluates the Lagrange basis derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //    ND must be at least 1.
        //
        //    Input, double XD[ND], the data points.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points.
        //
        //    Output, double LAGRANGE_DERIVATIVE[NI*ND], the values
        //    of the Lagrange basis derivatives at the interpolation points.
        //
    {
        double[] lpi = new double[ni * nd];

        for (int j = 0; j < nd; j++)
        {
            for (int i = 0; i < ni; i++)
            {
                lpi[i + j * ni] = 0.0;
            }
        }

        for (int i = 0; i < ni; i++)
        {
            for (int j = 0; j < nd; j++)
            {
                int j1;
                for (j1 = 0; j1 < nd; j1++)
                {
                    if (j1 != j)
                    {
                        double p = 1.0;
                        int j2;
                        for (j2 = 0; j2 < nd; j2++)
                        {
                            if (j2 == j1)
                            {
                                p /= (xd[j] - xd[j2]);
                            }
                            else if (j2 != j)
                            {
                                p = p * (xi[i] - xd[j2]) / (xd[j] - xd[j2]);
                            }
                        }

                        lpi[i + j * ni] += p;
                    }
                }
            }
        }

        return lpi;
    }

    public static double[] lagrange_value_OLD(int data_num, double[] t_data, int interp_num,
            double[] t_interp)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_VALUE evaluates the Lagrange polynomials.
        //
        //  Discussion:
        //
        //    Given DATA_NUM distinct abscissas, T_DATA(1:DATA_NUM),
        //    the I-th Lagrange polynomial L(I)(T) is defined as the polynomial of
        //    degree DATA_NUM - 1 which is 1 at T_DATA(I) and 0 at the DATA_NUM - 1
        //    other abscissas.
        //
        //    A formal representation is:
        //
        //      L(I)(T) = Product ( 1 <= J <= DATA_NUM, I /= J )
        //       ( T - T(J) ) / ( T(I) - T(J) )
        //
        //    This routine accepts a set of INTERP_NUM values at which all the Lagrange
        //    polynomials should be evaluated.
        //
        //    Given data values P_DATA at each of the abscissas, the value of the
        //    Lagrange interpolating polynomial at each of the interpolation points
        //    is then simple to compute by matrix multiplication:
        //
        //      P_INTERP(1:INTERP_NUM) =
        //        P_DATA(1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
        //
        //    or, in the case where P is multidimensional:
        //
        //      P_INTERP(1:M,1:INTERP_NUM) =
        //        P_DATA(1:M,1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DATA_NUM, the number of data points.
        //    DATA_NUM must be at least 1.
        //
        //    Input, double T_DATA[DATA_NUM], the data points.
        //
        //    Input, int INTERP_NUM, the number of
        //    interpolation points.
        //
        //    Input, double T_INTERP[INTERP_NUM], the
        //    interpolation points.
        //
        //    Output, double LAGRANGE_VALUE[DATA_NUM*INTERP_NUM], the values
        //    of the Lagrange polynomials at the interpolation points.
        //
    {
        int i;
        int i1;
        int i2;
        int j;
        double[] l_interp;

        l_interp = new double[data_num * interp_num];
        //
        //  Evaluate the polynomial.
        //
        for (j = 0; j < interp_num; j++)
        {
            for (i = 0; i < data_num; i++)
            {
                l_interp[i + j * data_num] = 1.0;
            }
        }

        for (i1 = 0; i1 < data_num; i1++)
        {
            for (i2 = 0; i2 < data_num; i2++)
            {
                if (i1 != i2)
                {
                    for (j = 0; j < interp_num; j++)
                    {
                        l_interp[i1 + j * data_num] = l_interp[i1 + j * data_num]
                            * (t_interp[j] - t_data[i2]) / (t_data[i1] - t_data[i2]);
                    }
                }
            }
        }

        return l_interp;
    }


    public static double[] lagrange_value(int nd, double[] xd, int ni, double[] xi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_VALUE evaluates the Lagrange basis polynomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //    ND must be at least 1.
        //
        //    Input, double XD[ND], the data points.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points.
        //
        //    Output, double LAGRANGE_BASIS[NI*ND], the values
        //    of the Lagrange basis polynomials at the interpolation points.
        //
    {
        //  Evaluate the polynomial.
        //
        double[] li = new double[ni * nd];

        for (int j = 0; j < nd; j++)
        {
            for (int i = 0; i < ni; i++)
            {
                li[i + j * ni] = 1.0;
            }
        }

        for (int i = 0; i < nd; i++)
        {
            for (int j = 0; j < nd; j++)
            {
                if (j != i)
                {
                    for (int k = 0; k < ni; k++)
                    {
                        li[k + i * ni] = li[k + i * ni] * (xi[k] - xd[j]) / (xd[i] - xd[j]);
                    }
                }
            }
        }

        return li;
    }
}