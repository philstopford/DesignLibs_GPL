using System;
using Burkardt.MatrixNS;

namespace Burkardt.Quadrature
{
    public static class HermiteQuadrature
    {
        public static void h_quadrature_rule ( int nt, ref double[] t, ref double[] wts )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    H_QUADRATURE_RULE: quadrature for H(i,x).
        //
        //  Discussion:
        //
        //    H(i,x) is the physicist's Hermite polynomial of degree I.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NT, the order of the rule.
        //
        //    Output, double T[NT], WTS[NT], the points and weights 
        //    of the rule.
        //
        {
            double[] bj;
            int i;
            const double r8_pi = 3.141592653589793;

            for ( i = 0; i < nt; i++ )
            {
                t[i] = 0.0;
            }

            bj = new double[nt];

            for ( i = 0; i < nt; i++ )
            {
                bj[i] = Math.Sqrt ( ( double ) ( i + 1 ) / 2.0 );
            }

            for ( i = 0; i < nt; i++ )
            {
                wts[i] = 0.0;
            }

            wts[0] = Math.Sqrt ( Math.Sqrt ( r8_pi ) );

            IMTQLX.imtqlx ( nt, ref t, ref bj, ref wts );

            for ( i = 0; i < nt; i++ )
            {
                wts[i] = wts[i] * wts[i];
            }
        }
        
        public static void he_quadrature_rule ( int nt, ref double[] t, ref double[] wts )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HE_QUADRATURE_RULE: quadrature for He(i,x).
        //
        //  Discussion:
        //
        //    He(i,x) represents the probabilist's Hermite polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NT, the order of the rule.
        //
        //    Output, double T[NT], WTS[NT], the points and weights 
        //    of the rule.
        //
        {
            double[] bj;
            int i;
            const double r8_pi = 3.141592653589793;

            for ( i = 0; i < nt; i++ )
            {
                t[i] = 0.0;
            }

            bj = new double[nt];

            for ( i = 0; i < nt; i++ )
            {
                bj[i] = Math.Sqrt ( ( double ) ( i + 1 ) / 2.0 );
            }

            for ( i = 0; i < nt; i++ )
            {
                wts[i] = 0.0;
            }
            wts[0] = Math.Sqrt ( Math.Sqrt ( r8_pi ) );

            IMTQLX.imtqlx ( nt, ref t, ref bj, ref wts );

            for ( i = 0; i < nt; i++ )
            {
                t[i] = t[i] * Math.Sqrt ( 2.0 );
            }
            for ( i = 0; i < nt; i++ )
            {
                wts[i] = wts[i] * wts[i] * Math.Sqrt ( 2.0 );
            }
        }
        
        public static void hf_quadrature_rule ( int nt, ref double[] t, ref double[] wts )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HF_QUADRATURE_RULE: quadrature for Hf(i,x).
        //
        //  Discussion:
        //
        //    Hf(i,x) represents the Hermite function of "degree" I.   
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NT, the order of the rule.
        //
        //    Output, double T[NT], WTS[NT], the points and weights 
        //    of the rule.
        //
        {
            double[] bj;
            int i;
            const double r8_pi = 3.141592653589793;

            for ( i = 0; i < nt; i++ )
            {
                t[i] = 0.0;
            }

            bj = new double[nt];

            for ( i = 0; i < nt; i++ )
            {
                bj[i] = Math.Sqrt ( ( double ) ( i + 1 ) / 2.0 );
            }

            for ( i = 0; i < nt; i++ )
            {
                wts[i] = 0.0;
            }

            wts[0] = Math.Sqrt ( Math.Sqrt ( r8_pi ) );

            IMTQLX.imtqlx ( nt, ref t, ref bj, ref wts );

            for ( i = 0; i < nt; i++ )
            {
                wts[i] = wts[i] * wts[i] * Math.Exp ( t[i] * t[i] );
            }
        }
    }
}