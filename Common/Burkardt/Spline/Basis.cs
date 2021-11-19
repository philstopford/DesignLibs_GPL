using System;
using Burkardt.Types;

namespace Burkardt.Spline;

public static class Basis
{
    public static double basis_function_b_val(double[] tdata, double tval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_FUNCTION_B_VAL evaluates the B spline basis function.
        //
        //  Discussion:
        //
        //    The B spline basis function is a piecewise cubic which
        //    has the properties that:
        //
        //    * it equals 2/3 at TDATA(3), 1/6 at TDATA(2) and TDATA(4);
        //    * it is 0 for TVAL <= TDATA(1) or TDATA(5) <= TVAL;
        //    * it is strictly increasing from TDATA(1) to TDATA(3),
        //      and strictly decreasing from TDATA(3) to TDATA(5);
        //    * the function and its first two derivatives are continuous
        //      at each node TDATA(I).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Davies, Philip Samuels,
        //    An Introduction to Computational Geometry for Curves and Surfaces,
        //    Clarendon Press, 1996,
        //    ISBN: 0-19-851478-6,
        //    LC: QA448.D38.
        //
        //  Input:
        //
        //    double TDATA(5), the nodes associated with the basis function.
        //    The entries of TDATA are assumed to be distinct and increasing.
        //
        //    double TVAL, a point at which the B spline basis function is
        //    to be evaluated.
        //
        //  Output:
        //
        //    double BASIS_FUNCTION_B_VAL, the value of the function at TVAL.
        //
    {
        const int NDATA = 5;

        int left = 0;
        int right = 0;
        double yval = 0;

        if (tval <= tdata[0] || tdata[NDATA - 1] <= tval)
        {
            yval = 0.0;
            return yval;
        }

        //
        //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
        //
        typeMethods.r8vec_bracket(NDATA, tdata, tval, ref left, ref right);
        //
        //  U is the normalized coordinate of TVAL in this interval.
        //
        double u = (tval - tdata[left - 1]) / (tdata[right - 1] - tdata[left - 1]);
        //
        //  Now evaluate the function.
        //
        if (tval < tdata[1])
        {
            yval = Math.Pow(u, 3) / 6.0;
        }
        else if (tval < tdata[2])
        {
            yval = (((-3.0
                        * u + 3.0)
                    * u + 3.0)
                * u + 1.0) / 6.0;
        }
        else if (tval < tdata[3])
        {
            yval = (((+3.0
                        * u - 6.0)
                    * u + 0.0)
                * u + 4.0) / 6.0;
        }
        else if (tval < tdata[4])
        {
            yval = Math.Pow(1.0 - u, 3) / 6.0;
        }

        return yval;

    }

    public static double basis_function_beta_val(double beta1, double beta2, double[] tdata,
            double tval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_FUNCTION_BETA_VAL evaluates the beta spline basis function.
        //
        //  Discussion:
        //
        //    With BETA1 = 1 and BETA2 = 0, the beta spline basis function
        //    equals the B spline basis function.
        //
        //    With BETA1 large, and BETA2 = 0, the beta spline basis function
        //    skews to the right, that is, its maximum increases, and occurs
        //    to the right of the center.
        //
        //    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
        //    a linear basis function; that is, its value in the outer two intervals
        //    goes to zero, and its behavior in the inner two intervals approaches
        //    a piecewise linear function that is 0 at the second node, 1 at the
        //    third, and 0 at the fourth.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Davies, Philip Samuels,
        //    An Introduction to Computational Geometry for Curves and Surfaces,
        //    Clarendon Press, 1996,
        //    ISBN: 0-19-851478-6,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, double BETA1, the skew or bias parameter.
        //    BETA1 = 1 for no skew or bias.
        //
        //    Input, double BETA2, the tension parameter.
        //    BETA2 = 0 for no tension.
        //
        //    Input, double TDATA[5], the nodes associated with the basis function.
        //    The entries of TDATA are assumed to be distinct and increasing.
        //
        //    Input, double TVAL, a point at which the B spline basis function is
        //    to be evaluated.
        //
        //    Output, double BASIS_FUNCTION_BETA_VAL, the value of the function at TVAL.
        //
    {
        const int NDATA = 5;

        double a;
        double b;
        double c;
        double d;
        int left = 0;
        int right = 0;
        double yval = 0;

        if (tval <= tdata[0] || tdata[NDATA - 1] <= tval)
        {
            yval = 0.0;
            return yval;
        }

        //
        //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
        //
        typeMethods.r8vec_bracket(NDATA, tdata, tval, ref left, ref right);
        //
        //  U is the normalized coordinate of TVAL in this interval.
        //
        double u = (tval - tdata[left - 1]) / (tdata[right - 1] - tdata[left - 1]);
        //
        //  Now evaluate the function.
        //
        if (tval < tdata[1])
        {
            yval = 2.0 * u * u * u;
        }
        else if (tval < tdata[2])
        {
            a = beta2 + 4.0 * beta1 + 4.0 * beta1 * beta1
                + 6.0 * (1.0 - beta1 * beta1)
                - 3.0 * (2.0 + beta2 + 2.0 * beta1)
                + 2.0 * (1.0 + beta2 + beta1 + beta1 * beta1);

            b = -6.0 * (1.0 - beta1 * beta1)
                + 6.0 * (2.0 + beta2 + 2.0 * beta1)
                - 6.0 * (1.0 + beta2 + beta1 + beta1 * beta1);

            c = -3.0 * (2.0 + beta2 + 2.0 * beta1)
                + 6.0 * (1.0 + beta2 + beta1 + beta1 * beta1);

            d = -2.0 * (1.0 + beta2 + beta1 + beta1 * beta1);

            yval = a + b * u + c * u * u + d * u * u * u;
        }
        else if (tval < tdata[3])
        {
            a = beta2 + 4.0 * beta1 + 4.0 * beta1 * beta1;

            b = -6.0 * beta1 * (1.0 - beta1 * beta1);

            c = -3.0 * (beta2 + 2.0 * beta1 * beta1
                              + 2.0 * beta1 * beta1 * beta1);

            d = 2.0 * (beta2 + beta1 + beta1 * beta1 + beta1 * beta1 * beta1);

            yval = a + b * u + c * u * u + d * u * u * u;
        }
        else if (tval < tdata[4])
        {
            yval = 2.0 * Math.Pow(beta1 * (1.0 - u), 3);
        }

        yval /= 2.0 + beta2 + 4.0 * beta1 + 4.0 * beta1 * beta1
                + 2.0 * beta1 * beta1 * beta1;

        return yval;
    }

    public static double[] basis_matrix_b_uni()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_B_UNI sets up the uniform B spline basis matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Foley, Andries vanDam, Steven Feiner, John Hughes,
        //    Computer Graphics, Principles and Practice,
        //    Second Edition,
        //    Addison Wesley, 1995,
        //    ISBN: 0201848406,
        //    LC: T385.C5735.
        //
        //  Parameters:
        //
        //    Output, double BASIS_MATRIX_B_UNI[4*4], the basis matrix.
        //
    {
        int j;
        double[] mbasis_save =
        {
            -1.0 / 6.0,
            3.0 / 6.0,
            -3.0 / 6.0,
            1.0 / 6.0,
            3.0 / 6.0,
            -6.0 / 6.0,
            0.0,
            4.0 / 6.0,
            -3.0 / 6.0,
            3.0 / 6.0,
            3.0 / 6.0,
            1.0 / 6.0,
            1.0 / 6.0,
            0.0,
            0.0,
            0.0
        };

        double[] mbasis = new double[4 * 4];

        for (j = 0; j < 4; j++)
        {
            int i;
            for (i = 0; i < 4; i++)
            {
                mbasis[i + j * 4] = mbasis_save[i + j * 4];
            }
        }

        return mbasis;
    }

    public static double[] basis_matrix_beta_uni(double beta1, double beta2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_BETA_UNI sets up the uniform beta spline basis matrix.
        //
        //  Discussion:
        //
        //    If BETA1 = 1 and BETA2 = 0, then the beta spline reduces to
        //    the B spline.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Foley, Andries vanDam, Steven Feiner, John Hughes,
        //    Computer Graphics, Principles and Practice,
        //    Second Edition,
        //    Addison Wesley, 1995,
        //    ISBN: 0201848406,
        //    LC: T385.C5735.
        //
        //  Parameters:
        //
        //    Input, double BETA1, the skew or bias parameter.
        //    BETA1 = 1 for no skew or bias.
        //
        //    Input, double BETA2, the tension parameter.
        //    BETA2 = 0 for no tension.
        //
        //    Output, double BASIS_MATRIX_BETA_UNI[4*4], the basis matrix.
        //
    {
        int j;

        double[] mbasis = new double[4 * 4];

        mbasis[0 + 0 * 4] = -2.0 * beta1 * beta1 * beta1;
        mbasis[0 + 1 * 4] = 2.0 * beta2
                            + 2.0 * beta1 * (beta1 * beta1 + beta1 + 1.0);
        mbasis[0 + 2 * 4] = -2.0 * (beta2 + beta1 * beta1 + beta1 + 1.0);
        mbasis[0 + 3 * 4] = 2.0;

        mbasis[1 + 0 * 4] = 6.0 * beta1 * beta1 * beta1;
        mbasis[1 + 1 * 4] = -3.0 * beta2
                            - 6.0 * beta1 * beta1 * (beta1 + 1.0);
        mbasis[1 + 2 * 4] = 3.0 * beta2 + 6.0 * beta1 * beta1;
        mbasis[1 + 3 * 4] = 0.0;

        mbasis[2 + 0 * 4] = -6.0 * beta1 * beta1 * beta1;
        mbasis[2 + 1 * 4] = 6.0 * beta1 * (beta1 - 1.0) * (beta1 + 1.0);
        mbasis[2 + 2 * 4] = 6.0 * beta1;
        mbasis[2 + 3 * 4] = 0.0;

        mbasis[3 + 0 * 4] = 2.0 * beta1 * beta1 * beta1;
        mbasis[3 + 1 * 4] = 4.0 * beta1 * (beta1 + 1.0) + beta2;
        mbasis[3 + 2 * 4] = 2.0;
        mbasis[3 + 3 * 4] = 0.0;

        double delta = ((2.0
                    * beta1 + 4.0)
                * beta1 + 4.0)
            * beta1 + 2.0 + beta2;

        for (j = 0; j < 4; j++)
        {
            int i;
            for (i = 0; i < 4; i++)
            {
                mbasis[i + j * 4] /= delta;
            }
        }

        return mbasis;
    }

    public static double[] basis_matrix_bezier()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_BEZIER_UNI sets up the cubic Bezier spline basis matrix.
        //
        //  Discussion:
        //
        //    This basis matrix assumes that the data points are stored as
        //    ( P1, P2, P3, P4 ).  P1 is the function value at T = 0, while
        //    P2 is used to approximate the derivative at T = 0 by
        //    dP/dt = 3 * ( P2 - P1 ).  Similarly, P4 is the function value
        //    at T = 1, and P3 is used to approximate the derivative at T = 1
        //    by dP/dT = 3 * ( P4 - P3 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Foley, Andries vanDam, Steven Feiner, John Hughes,
        //    Computer Graphics, Principles and Practice,
        //    Second Edition,
        //    Addison Wesley, 1995,
        //    ISBN: 0201848406,
        //    LC: T385.C5735.
        //
        //  Parameters:
        //
        //    Output, double BASIS_MATRIX_BEZIER[4*4], the basis matrix.
        //
    {
        double[] mbasis = new double[4 * 4];

        mbasis[0 + 0 * 4] = -1.0;
        mbasis[0 + 1 * 4] = 3.0;
        mbasis[0 + 2 * 4] = -3.0;
        mbasis[0 + 3 * 4] = 1.0;

        mbasis[1 + 0 * 4] = 3.0;
        mbasis[1 + 1 * 4] = -6.0;
        mbasis[1 + 2 * 4] = 3.0;
        mbasis[1 + 3 * 4] = 0.0;

        mbasis[2 + 0 * 4] = -3.0;
        mbasis[2 + 1 * 4] = 3.0;
        mbasis[2 + 2 * 4] = 0.0;
        mbasis[2 + 3 * 4] = 0.0;

        mbasis[3 + 0 * 4] = 1.0;
        mbasis[3 + 1 * 4] = 0.0;
        mbasis[3 + 2 * 4] = 0.0;
        mbasis[3 + 3 * 4] = 0.0;

        return mbasis;
    }

    public static double[] basis_matrix_hermite()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_HERMITE sets up the Hermite spline basis matrix.
        //
        //  Discussion:
        //
        //    This basis matrix assumes that the data points are stored as
        //    ( P1, P2, P1', P2' ), with P1 and P1' being the data value and 
        //    the derivative dP/dT at T = 0, while P2 and P2' apply at T = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Foley, Andries vanDam, Steven Feiner, John Hughes,
        //    Computer Graphics, Principles and Practice,
        //    Second Edition,
        //    Addison Wesley, 1995,
        //    ISBN: 0201848406,
        //    LC: T385.C5735.
        //
        //  Parameters:
        //
        //    Output, double BASIS_MATRIX_HERMITE[4*4], the basis matrix.
        //
    {
        double[] mbasis = new double[4 * 4];

        mbasis[0 + 0 * 4] = 2.0;
        mbasis[0 + 1 * 4] = -2.0;
        mbasis[0 + 2 * 4] = 1.0;
        mbasis[0 + 3 * 4] = 1.0;

        mbasis[1 + 0 * 4] = -3.0;
        mbasis[1 + 1 * 4] = 3.0;
        mbasis[1 + 2 * 4] = -2.0;
        mbasis[1 + 3 * 4] = -1.0;

        mbasis[2 + 0 * 4] = 0.0;
        mbasis[2 + 1 * 4] = 0.0;
        mbasis[2 + 2 * 4] = 1.0;
        mbasis[2 + 3 * 4] = 0.0;

        mbasis[3 + 0 * 4] = 1.0;
        mbasis[3 + 1 * 4] = 0.0;
        mbasis[3 + 2 * 4] = 0.0;
        mbasis[3 + 3 * 4] = 0.0;

        return mbasis;
    }

    public static double[] basis_matrix_overhauser_nonuni(double alpha, double beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_OVERHAUSER_NONUNI sets the nonuniform Overhauser spline basis matrix.
        //
        //  Discussion:
        //
        //    This basis matrix assumes that the data points P1, P2, P3 and
        //    P4 are not uniformly spaced in T, and that P2 corresponds to T = 0,
        //    and P3 to T = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ALPHA, BETA.
        //    ALPHA = || P2 - P1 || / ( || P3 - P2 || + || P2 - P1 || )
        //    BETA  = || P3 - P2 || / ( || P4 - P3 || + || P3 - P2 || ).
        //
        //    Output, double BASIS_MATRIX_OVERHAUSER_NONUNI[4*4], the basis matrix.
        //
    {
        double[] mbasis = new double[4 * 4];

        mbasis[0 + 0 * 4] = -(1.0 - alpha) * (1.0 - alpha) / alpha;
        mbasis[0 + 1 * 4] = beta + (1.0 - alpha) / alpha;
        mbasis[0 + 2 * 4] = alpha - 1.0 / (1.0 - beta);
        mbasis[0 + 3 * 4] = beta * beta / (1.0 - beta);

        mbasis[1 + 0 * 4] = 2.0 * (1.0 - alpha) * (1.0 - alpha) / alpha;
        mbasis[1 + 1 * 4] = (-2.0 * (1.0 - alpha) - alpha * beta) / alpha;
        mbasis[1 + 2 * 4] = (2.0 * (1.0 - alpha)
                             - beta * (1.0 - 2.0 * alpha)) / (1.0 - beta);
        mbasis[1 + 3 * 4] = -beta * beta / (1.0 - beta);

        mbasis[2 + 0 * 4] = -(1.0 - alpha) * (1.0 - alpha) / alpha;
        mbasis[2 + 1 * 4] = (1.0 - 2.0 * alpha) / alpha;
        mbasis[2 + 2 * 4] = alpha;
        mbasis[2 + 3 * 4] = 0.0;

        mbasis[3 + 0 * 4] = 0.0;
        mbasis[3 + 1 * 4] = 1.0;
        mbasis[3 + 2 * 4] = 0.0;
        mbasis[3 + 3 * 4] = 0.0;

        return mbasis;
    }

    public static double[] basis_matrix_overhauser_nul(double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_OVERHAUSER_NUL sets the nonuniform left Overhauser spline basis matrix.
        //
        //  Discussion:
        //
        //    This basis matrix assumes that the data points P1, P2, and
        //    P3 are not uniformly spaced in T, and that P1 corresponds to T = 0,
        //    and P2 to T = 1. (???)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ALPHA.
        //    ALPHA = || P2 - P1 || / ( || P3 - P2 || + || P2 - P1 || )
        //
        //    Output, double BASIS_MATRIX_OVERHAUSER_NUL[3*3], the basis matrix.
        //
    {
        double[] mbasis = new double[3 * 3];

        mbasis[0 + 0 * 3] = 1.0 / alpha;
        mbasis[0 + 1 * 3] = -1.0 / (alpha * (1.0 - alpha));
        mbasis[0 + 2 * 3] = 1.0 / (1.0 - alpha);

        mbasis[1 + 0 * 3] = -(1.0 + alpha) / alpha;
        mbasis[1 + 1 * 3] = 1.0 / (alpha * (1.0 - alpha));
        mbasis[1 + 2 * 3] = -alpha / (1.0 - alpha);

        mbasis[2 + 0 * 3] = 1.0;
        mbasis[2 + 1 * 3] = 0.0;
        mbasis[2 + 2 * 3] = 0.0;

        return mbasis;
    }

    public static double[] basis_matrix_overhauser_nur(double beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_OVERHAUSER_NUR sets the nonuniform right Overhauser spline basis matrix.
        //
        //  Discussion:
        //
        //    This basis matrix assumes that the data points PN-2, PN-1, and
        //    PN are not uniformly spaced in T, and that PN-1 corresponds to T = 0,
        //    and PN to T = 1. (???)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double BETA.
        //    BETA = || P(N) - P(N-1) || / ( || P(N) - P(N-1) || + || P(N-1) - P(N-2) || )
        //
        //    Output, double BASIS_MATRIX_OVERHAUSER_NUR[3*3], the basis matrix.
        //
    {
        double[] mbasis = new double[3 * 3];

        mbasis[0 + 0 * 3] = 1.0 / beta;
        mbasis[0 + 1 * 3] = -1.0 / (beta * (1.0 - beta));
        mbasis[0 + 2 * 3] = 1.0 / (1.0 - beta);

        mbasis[1 + 0 * 3] = -(1.0 + beta) / beta;
        mbasis[1 + 1 * 3] = 1.0 / (beta * (1.0 - beta));
        mbasis[1 + 2 * 3] = -beta / (1.0 - beta);

        mbasis[2 + 0 * 3] = 1.0;
        mbasis[2 + 1 * 3] = 0.0;
        mbasis[2 + 2 * 3] = 0.0;

        return mbasis;
    }

    public static double[] basis_matrix_overhauser_uni()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_OVERHAUSER_UNI sets the uniform Overhauser spline basis matrix.
        //
        //  Discussion:
        //
        //    This basis matrix assumes that the data points P1, P2, P3 and
        //    P4 are uniformly spaced in T, and that P2 corresponds to T = 0,
        //    and P3 to T = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Foley, Andries vanDam, Steven Feiner, John Hughes,
        //    Computer Graphics, Principles and Practice,
        //    Second Edition,
        //    Addison Wesley, 1995,
        //    ISBN: 0201848406,
        //    LC: T385.C5735.
        //
        //  Parameters:
        //
        //    Output, double BASIS_MATRIX_OVERHASUER_UNI[4*4], the basis matrix.
        //
    {
        double[] mbasis = new double[4 * 4];

        mbasis[0 + 0 * 4] = -1.0 / 2.0;
        mbasis[0 + 1 * 4] = 3.0 / 2.0;
        mbasis[0 + 2 * 4] = -3.0 / 2.0;
        mbasis[0 + 3 * 4] = 1.0 / 2.0;

        mbasis[1 + 0 * 4] = 2.0 / 2.0;
        mbasis[1 + 1 * 4] = -5.0 / 2.0;
        mbasis[1 + 2 * 4] = 4.0 / 2.0;
        mbasis[1 + 3 * 4] = -1.0 / 2.0;

        mbasis[2 + 0 * 4] = -1.0 / 2.0;
        mbasis[2 + 1 * 4] = 0.0;
        mbasis[2 + 2 * 4] = 1.0 / 2.0;
        mbasis[2 + 3 * 4] = 0.0;

        mbasis[3 + 0 * 4] = 0.0;
        mbasis[3 + 1 * 4] = 2.0 / 2.0;
        mbasis[3 + 2 * 4] = 0.0;
        mbasis[3 + 3 * 4] = 0.0;

        return mbasis;
    }

    public static double[] basis_matrix_overhauser_uni_l()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_OVERHAUSER_UNI_L sets the left uniform Overhauser spline basis matrix.
        //
        //  Discussion:
        //
        //    This basis matrix assumes that the data points P1, P2, and P3
        //    are not uniformly spaced in T, and that P1 corresponds to T = 0,
        //    and P2 to T = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double BASIS_MATRIX_OVERHASUER_UNI_L[3*3], the basis matrix.
        //
    {
        double[] mbasis = new double[3 * 3];

        mbasis[0 + 0 * 3] = 2.0;
        mbasis[0 + 1 * 3] = -4.0;
        mbasis[0 + 2 * 3] = 2.0;

        mbasis[1 + 0 * 3] = -3.0;
        mbasis[1 + 1 * 3] = 4.0;
        mbasis[1 + 2 * 3] = -1.0;

        mbasis[2 + 0 * 3] = 1.0;
        mbasis[2 + 1 * 3] = 0.0;
        mbasis[2 + 2 * 3] = 0.0;

        return mbasis;
    }

    public static double[] basis_matrix_overhauser_uni_r()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_OVERHAUSER_UNI_R sets the right uniform Overhauser spline basis matrix.
        //
        //  Discussion:
        //
        //    This basis matrix assumes that the data points P(N-2), P(N-1),
        //    and P(N) are uniformly spaced in T, and that P(N-1) corresponds to
        //    T = 0, and P(N) to T = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double BASIS_MATRIX_OVERHASUER_UNI_R[3*3], the basis matrix.
        //
    {
        double[] mbasis = new double[3 * 3];

        mbasis[0 + 0 * 3] = 2.0;
        mbasis[0 + 1 * 3] = -4.0;
        mbasis[0 + 2 * 3] = 2.0;

        mbasis[1 + 0 * 3] = -3.0;
        mbasis[1 + 1 * 3] = 4.0;
        mbasis[1 + 2 * 3] = -1.0;

        mbasis[2 + 0 * 3] = 1.0;
        mbasis[2 + 1 * 3] = 0.0;
        mbasis[2 + 2 * 3] = 0.0;

        return mbasis;
    }

    public static double basis_matrix_tmp(int left, int n, double[] mbasis, int ndata,
            double[] tdata, double[] ydata, double tval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MATRIX_TMP computes Q = T * MBASIS * P
        //
        //  Discussion:
        //
        //    YDATA is a vector of data values, most frequently the values of some
        //    function sampled at uniformly spaced points.  MBASIS is the basis
        //    matrix for a particular kind of spline.  T is a vector of the
        //    powers of the normalized difference between TVAL and the left
        //    endpoint of the interval.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEFT, indicats that TVAL is in the interval
        //    [ TDATA(LEFT), TDATA(LEFT+1) ], or that this is the "nearest"
        //    interval to TVAL.
        //    For TVAL < TDATA(1), use LEFT = 1.
        //    For TDATA(NDATA) < TVAL, use LEFT = NDATA - 1.
        //
        //    Input, int N, the order of the basis matrix.
        //
        //    Input, double MBASIS[N*N], the basis matrix.
        //
        //    Input, int NDATA, the dimension of the vectors TDATA and YDATA.
        //
        //    Input, double TDATA[NDATA], the abscissa values.  This routine
        //    assumes that the TDATA values are uniformly spaced, with an
        //    increment of 1.0.
        //
        //    Input, double YDATA[NDATA], the data values to be interpolated or
        //    approximated.
        //
        //    Input, double TVAL, the value of T at which the spline is to be
        //    evaluated.
        //
        //    Output, double BASIS_MATRIX_TMP, the value of the spline at TVAL.
        //
    {
        double arg = 0;
        int first = 0;
        int i;
        int j;

        double[] tvec = new double[n];

        switch (left)
        {
            case 1:
                arg = 0.5 * (tval - tdata[left - 1]);
                first = left;
                break;
            default:
            {
                if (left < ndata - 1)
                {
                    arg = tval - tdata[left - 1];
                    first = left - 1;
                }
                else if (left == ndata - 1)
                {
                    arg = 0.5 * (1.0 + tval - tdata[left - 1]);
                    first = left - 1;
                }

                break;
            }
        }

        //
        //  TVEC(I) = ARG**(N-I).
        //
        tvec[n - 1] = 1.0;
        for (i = n - 2; 0 <= i; i--)
        {
            tvec[i] = arg * tvec[i + 1];
        }

        double yval = 0.0;
        for (j = 0; j < n; j++)
        {
            double tm = 0.0;
            for (i = 0; i < n; i++)
            {
                tm += tvec[i] * mbasis[i + j * n];
            }

            yval += tm * ydata[first - 1 + j];
        }

        return yval;
    }
}