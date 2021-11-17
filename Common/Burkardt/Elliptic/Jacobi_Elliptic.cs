using System;

namespace Burkardt.Elliptic;

public static class Jacobi
{
    public static void sncndn(double u, double m, ref double sn, ref double cn, ref double dn )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SNCNDN evaluates Jacobi elliptic functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 June 2018
        //
        //  Author:
        //
        //    Original ALGOL version by Roland Bulirsch.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Roland Bulirsch,
        //    Numerical calculation of elliptic integrals and elliptic functions,
        //    Numerische Mathematik,
        //    Volume 7, Number 1, 1965, pages 78-90.
        //
        //  Parameters:
        //
        //    Input, double U, M, the arguments.
        //
        //    Output, double &SN, &CN, &DN, the value of the Jacobi
        //    elliptic functions sn(u,m), cn(u,m), and dn(u,m).
        //
    {
        double a;
        double b;
        double c = 0;
        double ca;
        double d = 0;
        int i;
        int l;
        double[] m_array;
        double m_comp;
        double[] n_array;
        const double r8_epsilon = 2.220446049250313E-16;
        double u_copy;

        m_comp = 1.0 - m;
        u_copy = u;

        switch (m_comp)
        {
            case 0.0:
                cn = 1.0 / Math.Cosh(u_copy);
                dn = cn;
                sn = Math.Tanh(u_copy);
                return;
        }

        switch (m)
        {
            case > 1.0:
                d = 1.0 - m_comp;
                m_comp = -m_comp / d;
                d = Math.Sqrt(d);
                u_copy = d * u_copy;
                break;
        }

        ca = Math.Sqrt(r8_epsilon);

        a = 1.0;
        dn = 1.0;
        l = 24;

        m_array = new double[25];
        n_array = new double[25];

        for (i = 0; i < 25; i++)
        {
            m_array[i] = a;
            m_comp = Math.Sqrt(m_comp);
            n_array[i] = m_comp;
            c = 0.5 * (a + m_comp);
            if (Math.Abs(a - m_comp) <= ca * a)
            {
                l = i;
                break;
            }

            m_comp = a * m_comp;
            a = c;
        }

        u_copy = c * u_copy;
        sn = Math.Sin(u_copy);
        cn = Math.Cos(u_copy);

        if (sn != 0.0)
        {
            a = cn / sn;
            c = a * c;

            for (i = l; 0 <= i; i--)
            {
                b = m_array[i];
                a = c * a;
                c = dn * c;
                dn = (n_array[i] + a) / (b + a);
                a = c / b;
            }

            a = 1.0 / Math.Sqrt(c * c + 1.0);

            sn = sn switch
            {
                < 0.0 => -a,
                _ => a
            };

            cn = c * sn;
        }

        switch (m)
        {
            case > 1.0:
                a = dn;
                dn = cn;
                cn = a;
                sn /= d;
                break;
        }
    }
}