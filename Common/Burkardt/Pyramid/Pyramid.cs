using Burkardt.Types;

namespace Burkardt.PyramidNS
{
    public static class Pyramid
    {
        public static double pyra_unit_monomial(int[] expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRA_UNIT_MONOMIAL: monomial integral in a unit pyramid.
            //
            //  Discussion:
            //
            //    This function returns the value of the integral of X^ALPHA Y^BETA Z^GAMMA
            //    over the unit pyramid.
            //
            //    The integration region is:
            //
            //    - ( 1 - Z ) <= X <= 1 - Z
            //    - ( 1 - Z ) <= Y <= 1 - Z
            //              0 <= Z <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, int EXPON[3], the exponents.
            //
            //    Output, double PYRA_UNIT_MONOMIAL, the integral of the monomial
            //    over the pyramid.
            //
        {
            int i;
            int i_hi;
            double value;

            value = 0.0;

            if ((expon[0] % 2) == 0 && (expon[1] % 2) == 0)
            {
                i_hi = 2 + expon[0] + expon[1];

                for (i = 0; i <= i_hi; i++)
                {
                    value = value + typeMethods.r8_mop(i) * typeMethods.r8_choose(i_hi, i)
                        / (double)(i + expon[2] + 1);
                }

                value = value
                    * 2.0 / (double)(expon[0] + 1)
                    * 2.0 / (double)(expon[1] + 1);
            }

            return value;
        }

        public static double pyra_unit_volume()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRA_UNIT_VOLUME: volume of a unit pyramid with square base.
            //
            //  Discussion:
            //
            //    The volume of this unit pyramid is 4/3.
            //
            //    The integration region is:
            //
            //      - ( 1 - Z ) <= X <= 1 - Z
            //      - ( 1 - Z ) <= Y <= 1 - Z
            //                0 <= Z <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double PYRA_UNIT_VOLUME, the volume of the pyramid.
            //
        {
            double volume;

            volume = 4.0 / 3.0;

            return volume;
        }
    }
}