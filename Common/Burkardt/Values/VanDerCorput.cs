﻿namespace Burkardt.Values
{
    public static class VanDerCorput
    {
        public static void van_der_corput_values(ref int n_data, ref int base_, ref int seed, ref double value)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VAN_DER_CORPUT_VALUES returns some values of the van der Corput sequence.
            //
            //  Discussion:
            //
            //    The van der Corput sequence is often used to generate a "subrandom"
            //    sequence of points which have a better covering property
            //    than pseudorandom points.
            //
            //    The van der Corput sequence generates a sequence of points in [0,1]
            //    which (theoretically) never repeats.  Except for SEED = 0, the
            //    elements of the van der Corput sequence are strictly between 0 and 1.
            //
            //    The van der Corput sequence writes an int *in a given base B,
            //    and then its digits are "reflected" about the decimal point.
            //    This maps the numbers from 1 to N into a set of numbers in [0,1],
            //    which are especially nicely distributed if N is one less
            //    than a power of the base.
            //
            //    Hammersley suggested generating a set of N nicely distributed
            //    points in two dimensions by setting the first component of the
            //    Ith point to I/N, and the second to the van der Corput
            //    value of I in base 2.
            //
            //    Halton suggested that in many cases, you might not know the number
            //    of points you were generating, so Hammersley's formulation was
            //    not ideal.  Instead, he suggested that to generate a nicely
            //    distributed sequence of points in M dimensions, you simply
            //    choose the first M primes, P(1:M), and then for the J-th component of
            //    the I-th point in the sequence, you compute the van der Corput
            //    value of I in base P(J).
            //
            //    Thus, to generate a Halton sequence in a 2 dimensional space,
            //    it is typical practice to generate a pair of van der Corput sequences,
            //    the first with prime base 2, the second with prime base 3.
            //    Similarly, by using the first K primes, a suitable sequence
            //    in K-dimensional space can be generated.
            //
            //    The generation is quite simple.  Given an int *SEED, the expansion
            //    of SEED in base BASE is generated.  Then, essentially, the result R
            //    is generated by writing a decimal point followed by the digits of
            //    the expansion of SEED, in reverse order.  This decimal value is actually
            //    still in base BASE, so it must be properly interpreted to generate
            //    a usable value.
            //
            //  Example:
            //
            //    BASE = 2
            //
            //    SEED     SEED      van der Corput
            //    decimal  binary    binary   decimal
            //    -------  ------    ------   -------
            //        0  =     0  =>  .0     = 0.0;
            //        1  =     1  =>  .1     = 0.5
            //        2  =    10  =>  .01    = 0.25
            //        3  =    11  =>  .11    = 0.75
            //        4  =   100  =>  .001   = 0.125
            //        5  =   101  =>  .101   = 0.625
            //        6  =   110  =>  .011   = 0.375
            //        7  =   111  =>  .111   = 0.875
            //        8  =  1000  =>  .0001  = 0.0625
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 June 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    John Halton,
            //    On the efficiency of certain quasi-random sequences of points
            //    in evaluating multi-dimensional integrals,
            //    Numerische Mathematik,
            //    Volume 2, pages 84-90, 1960.
            //
            //    John Hammersley,
            //    Monte Carlo methods for solving multivariable problems,
            //    Proceedings of the New York Academy of Science,
            //    Volume 86, pages 844-874, 1960.
            //
            //    J G van der Corput,
            //    Verteilungsfunktionen,
            //    Proc Akad Amsterdam,
            //    Volume 38, 1935,
            //    Volume 39, 1936.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int BASE, the base of the sequence.
            //
            //    Output, ref int SEED, the index of the element of the sequence.
            //
            //    Output, ref double VALUE, the value of the SEED-th element of the
            //    van der Corput sequence in base BASE.
            //
        {
            int N_MAX = 75;

            int[] base_vec =
            {
                2, 2, 2, 2, 2,
                2, 2, 2, 2, 3,
                3, 3, 3, 3, 3,
                3, 3, 3, 4, 4,
                4, 4, 4, 4, 4,
                4, 4, 2, 3, 4,
                5, 7, 11, 13, 2,
                3, 4, 5, 7, 11,
                13, 2, 3, 4, 5,
                7, 11, 13, 2, 3,
                4, 5, 7, 11, 13,
                29, 29, 29, 29, 29,
                71, 71, 71, 71, 71,
                173, 173, 173, 173, 173,
                409, 409, 409, 409, 409
            };

            int[] seed_vec =
            {
                0, 1, 2, 3, 4,
                5, 6, 7, 8, 0,
                1, 2, 3, 4, 5,
                6, 7, 8, 0, 1,
                2, 3, 4, 5, 6,
                7, 8, 10, 10, 10,
                10, 10, 10, 10, 100,
                100, 100, 100, 100, 100,
                100, 1000, 1000, 1000, 1000,
                1000, 1000, 1000, 10000, 10000,
                10000, 10000, 10000, 10000, 10000,
                1000, 1001, 1002, 1003, 1004,
                1000, 1001, 1002, 1003, 1004,
                1000, 1001, 1002, 1003, 1004,
                1000, 1001, 1002, 1003, 1004
            };

            double[] value_vec =
            {
                0.0000000000000000E+00,
                0.5000000000000000E+00,
                0.2500000000000000E+00,
                0.7500000000000000E+00,
                0.1250000000000000E+00,
                0.6250000000000000E+00,
                0.3750000000000000E+00,
                0.8750000000000000E+00,
                0.0625000000000000E+00,
                0.0000000000000000E+00,
                0.3333333333333333E+00,
                0.6666666666666666E+00,
                0.1111111111111111E+00,
                0.4444444444444444E+00,
                0.7777777777777777E+00,
                0.2222222222222222E+00,
                0.5555555555555556E+00,
                0.8888888888888888E+00,
                0.0000000000000000E+00,
                0.2500000000000000E+00,
                0.5000000000000000E+00,
                0.7500000000000000E+00,
                0.0625000000000000E+00,
                0.3125000000000000E+00,
                0.5625000000000000E+00,
                0.8125000000000000E+00,
                0.1250000000000000E+00,
                0.3125000000000000E+00,
                0.3703703703703703E+00,
                0.6250000000000000E+00,
                0.0800000000000000E+00,
                0.4489795918367347E+00,
                0.9090909090909092E+00,
                0.7692307692307693E+00,
                0.1484375000000000E+00,
                0.4115226337448559E+00,
                0.0976562500000000E+00,
                0.0320000000000000E+00,
                0.2915451895043731E+00,
                0.1652892561983471E+00,
                0.7337278106508875E+00,
                0.0927734375000000E+00,
                0.3475080018289895E+00,
                0.1708984375000000E+00,
                0.0051200000000000E+00,
                0.9162848812994586E+00,
                0.9316303531179565E+00,
                0.9904415111515704E+00,
                0.0347290039062500E+00,
                0.3861200020322105E+00,
                0.0189208984375000E+00,
                0.0005120000000000E+00,
                0.5749985125245433E+00,
                0.1529950140017758E+00,
                0.2459297643639929E+00,
                0.4887449259912255E+00,
                0.5232276846119153E+00,
                0.5577104432326049E+00,
                0.5921932018532945E+00,
                0.6266759604739842E+00,
                0.0872842689942472E+00,
                0.1013687760365007E+00,
                0.1154532830787542E+00,
                0.1295377901210077E+00,
                0.1436222971632613E+00,
                0.7805138828560928E+00,
                0.7862942296769020E+00,
                0.7920745764977113E+00,
                0.7978549233185205E+00,
                0.8036352701393298E+00,
                0.4449997309915651E+00,
                0.4474447187666262E+00,
                0.4498897065416874E+00,
                0.4523346943167484E+00,
                0.4547796820918096E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                base_ = 0;
                seed = 0;
                value = 0.0;
            }
            else
            {
                base_ = base_vec[n_data - 1];
                seed = seed_vec[n_data - 1];
                value = value_vec[n_data - 1];
            }
        }

    }
}