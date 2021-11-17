using System;

namespace utility;

// https://stackoverflow.com/questions/1552738/is-there-a-java-equivalent-of-frexp
// This one seems to work better than DoubleConverter for some users :
// mantissa is a double, exponent is int.
public static class FrExp
{
    public class FRexpResult
    {
        public int exponent;
        public double mantissa;
    }

    public static FRexpResult calculate(double value)
    {
        FRexpResult result = new();
        long bits = BitConverter.DoubleToInt64Bits(value);

        // Test for NaN, infinity, and zero.
        if (double.IsNaN(value) ||
            Math.Abs(value + value - value) <= double.Epsilon ||
            double.IsInfinity(value))
        {
            result.exponent = 0;
            result.mantissa = value;
        }
        else
        {

            bool neg = bits < 0;
            int exponent = (int)((bits >> 52) & 0x7ffL);
            long mantissa = bits & 0xfffffffffffffL;

            switch (exponent)
            {
                case 0:
                    exponent++;
                    break;
                default:
                    mantissa |= (1L << 52);
                    break;
            }

            // bias the exponent - actually biased by 1023.
            // we are treating the mantissa as m.0 instead of 0.m
            //  so subtract another 52.
            exponent -= 1075;
            double realMant = mantissa;

            // normalize
            while (realMant > 1.0)
            {
                mantissa >>= 1;
                realMant /= 2.0;
                exponent++;
            }

            switch (neg)
            {
                case true:
                    realMant *= -1;
                    break;
            }

            result.exponent = exponent;
            result.mantissa = realMant;
        }

        return result;
    }
}