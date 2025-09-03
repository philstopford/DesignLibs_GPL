using System;
using System.Runtime.CompilerServices;

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

    // Performance optimization: Use object pooling for results to reduce allocations
    [ThreadStatic] private static FRexpResult? _cachedResult;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static FRexpResult calculate(double value)
    {
        // Performance optimization: Reuse result object to reduce allocations
        FRexpResult result = _cachedResult ??= new FRexpResult();
        
        // Performance optimization: Fast path for special values
        if (value == 0.0)
        {
            result.exponent = 0;
            result.mantissa = 0.0;
            return result;
        }
        
        if (double.IsNaN(value))
        {
            result.exponent = 0;
            result.mantissa = double.NaN;
            return result;
        }
        
        if (double.IsInfinity(value))
        {
            result.exponent = 0;
            result.mantissa = value;
            return result;
        }

        long bits = BitConverter.DoubleToInt64Bits(value);
        bool neg = bits < 0;
        int exponent = (int)((bits >> 52) & 0x7ffL);
        long mantissa = bits & 0xfffffffffffffL;

        switch (exponent)
        {
            case 0:
                exponent++;
                break;
            default:
                mantissa |= 1L << 52;
                break;
        }

        // bias the exponent - actually biased by 1023.
        // we are treating the mantissa as m.0 instead of 0.m
        //  so subtract another 52.
        exponent -= 1075;
        double realMant = mantissa;

        // Performance optimization: Use bit shifting for normalization when possible
        while (realMant > 1.0)
        {
            mantissa >>= 1;
            realMant *= 0.5; // Multiplication is typically faster than division
            exponent++;
        }

        if (neg)
        {
            realMant = -realMant;
        }

        result.exponent = exponent;
        result.mantissa = realMant;

        return result;
    }
}