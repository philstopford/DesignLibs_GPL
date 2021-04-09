using System;

namespace utility
{
    // https://stackoverflow.com/questions/1552738/is-there-a-java-equivalent-of-frexp
    // This one seems to work better than DoubleConverter for some users :
    // mantissa is a double, exponent is int.
    public static class FrExp
    {
        public class FRexpResult
        {
            public int exponent = 0;
            public double mantissa = 0.0;
        }

        public static FRexpResult calculate(double value)
        {
            FRexpResult result = new FRexpResult();
            long bits = BitConverter.DoubleToInt64Bits(value);
            double realMant = 1.0;

            // Test for NaN, infinity, and zero.
            if (Double.IsNaN(value) ||
                Math.Abs(value + value - value) < Double.Epsilon ||
                Double.IsInfinity(value))
            {
                result.exponent = 0;
                result.mantissa = value;
            }
            else
            {

                bool neg = (bits < 0);
                int exponent = (int)((bits >> 52) & 0x7ffL);
                long mantissa = bits & 0xfffffffffffffL;

                if (exponent == 0)
                {
                    exponent++;
                }
                else
                {
                    mantissa = mantissa | (1L << 52);
                }

                // bias the exponent - actually biased by 1023.
                // we are treating the mantissa as m.0 instead of 0.m
                //  so subtract another 52.
                exponent -= 1075;
                realMant = mantissa;

                // normalize
                while (realMant > 1.0)
                {
                    mantissa >>= 1;
                    realMant /= 2.0;
                    exponent++;
                }

                if (neg)
                {
                    realMant = realMant * -1;
                }

                result.exponent = exponent;
                result.mantissa = realMant;
            }

            return result;
        }
    }
}
