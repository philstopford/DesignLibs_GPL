using System;
using System.Globalization;

namespace MiscUtil.Conversion;

/// <summary>
/// A class to allow the conversion of doubles to string representations of
/// their exact decimal values. The implementation aims for readability over
/// efficiency.
/// </summary>
public class DoubleConverter
{
    /// <doc> This provides the results from the convertor - the exponent, mantissa and whether the result is negative or not.</doc>
    public class Result
    {
        /// <summary>
        /// exponent
        /// </summary>
        public int exponent;
        /// <summary>
        /// mantissa
        /// </summary>
        public long mantissa;
        /// <summary>
        /// negative
        /// </summary>
        public bool negative;
    }

    /// <summary>
    /// Takes a double and returns a result from the conversion
    /// </summary>
    /// <param name="d"></param>
    /// <returns></returns>
    public static Result calculate(double d)
    {
        Result r = new();
        // Translate the double into sign, exponent and mantissa.
        long bits = BitConverter.DoubleToInt64Bits(d);
        // Note that the shift is sign-extended, hence the test against -1 not 1
        r.negative = bits < 0;
        int exponent = (int)((bits >> 52) & 0x7ffL);
        long mantissa = bits & 0xfffffffffffffL;

        switch (exponent)
        {
            // Subnormal numbers; exponent is effectively one higher,
            // but there's no extra normalisation bit in the mantissa
            case 0:
                exponent++;
                break;
            // Normal numbers; leave exponent as it is but add extra
            default:
                mantissa |= 1L << 52;
                break;
        }

        // Bias the exponent. It's actually biased by 1023, but we're
        // treating the mantissa as m.0 rather than 0.m, so we need
        // to subtract another 52 from it.
        exponent -= 1075;

        switch (mantissa)
        {
            case 0:
                r.mantissa = 0;
                break;
            default:
            {
                /* Normalize */
                while ((mantissa & 1) == 0)
                {    /*  i.e., Mantissa is even */
                    mantissa >>= 1;
                    exponent++;
                }

                break;
            }
        }
        r.exponent = exponent;
        r.mantissa = mantissa;

        return r;
    }

    /// <summary>
    /// Converts the given double to a string representation of its
    /// exact decimal value.
    /// </summary>
    /// <param name="d">The double to convert.</param>
    /// <returns>A string representation of the double's exact decimal value.</returns>
    public static string ToExactString(double d)
    {
        if (double.IsPositiveInfinity(d))
        {
            return "+Infinity";
        }

        if (double.IsNegativeInfinity(d))
        {
            return "-Infinity";
        }

        switch (d)
        {
            case double.NaN:
                return "NaN";
        }

        Result r = calculate(d);

        // Construct a new decimal expansion with the mantissa
        ArbitraryDecimal ad = new(r.mantissa);

        switch (r.exponent)
        {
            // If the exponent is less than 0, we need to repeatedly
            // divide by 2 - which is the equivalent of multiplying
            // by 5 and dividing by 10.
            case < 0:
            {
                for (int i = 0; i < -r.exponent; i++)
                    ad.MultiplyBy(5);
                ad.Shift(-r.exponent);
                break;
            }
            // Otherwise, we need to repeatedly multiply by 2
            default:
            {
                for (int i = 0; i < r.exponent; i++)
                    ad.MultiplyBy(2);
                break;
            }
        }

        return r.negative switch
        {
            // Finally, return the string with an appropriate sign
            true => "-" + ad,
            _ => ad.ToString()
        };
    }

    /// <summary>
    /// Converts the given double to a string representation of its
    /// exact decimal value.
    /// </summary>
    /// <param>The double to convert.</param>
    /// <returns>A string representation of the double's exact decimal value.</returns>

    /// <summary>
    /// Private class used for manipulating sequences of decimal digits.
    /// </summary>
    private class ArbitraryDecimal
    {
        /// <summary>Digits in the decimal expansion, one byte per digit</summary>
        private byte[] digits;
        /// <summary> 
        /// How many digits are *after* the decimal point
        /// </summary>
        private int decimalPoint;

        /// <summary> 
        /// Constructs an arbitrary decimal expansion from the given long.
        /// The long must not be negative.
        /// </summary>
        internal ArbitraryDecimal(long x)
        {
            string tmp = x.ToString(CultureInfo.InvariantCulture);
            digits = new byte[tmp.Length];
            for (int i = 0; i < tmp.Length; i++)
                digits[i] = (byte)(tmp[i] - '0');
            Normalize();
        }

        /// <summary>
        /// Multiplies the current expansion by the given amount, which should
        /// only be 2 or 5.
        /// </summary>
        internal void MultiplyBy(int amount)
        {
            byte[] result = new byte[digits.Length + 1];
            for (int i = digits.Length - 1; i >= 0; i--)
            {
                int resultDigit = digits[i] * amount + result[i + 1];
                result[i] = (byte)(resultDigit / 10);
                result[i + 1] = (byte)(resultDigit % 10);
            }
            if (result[0] != 0)
            {
                digits = result;
            }
            else
            {
                Array.Copy(result, 1, digits, 0, digits.Length);
            }
            Normalize();
        }

        /// <summary>
        /// Shifts the decimal point; a negative value makes
        /// the decimal expansion bigger (as fewer digits come after the
        /// decimal place) and a positive value makes the decimal
        /// expansion smaller.
        /// </summary>
        internal void Shift(int amount)
        {
            decimalPoint += amount;
        }

        /// <summary>
        /// Removes leading/trailing zeroes from the expansion.
        /// </summary>
        internal void Normalize()
        {
            int first;
            for (first = 0; first < digits.Length; first++)
                if (digits[first] != 0)
                {
                    break;
                }

            int last;
            for (last = digits.Length - 1; last >= 0; last--)
                if (digits[last] != 0)
                {
                    break;
                }

            switch (first)
            {
                case 0 when last == digits.Length - 1:
                    return;
            }

            byte[] tmp = new byte[last - first + 1];
            for (int i = 0; i < tmp.Length; i++)
                tmp[i] = digits[i + first];

            decimalPoint -= digits.Length - (last + 1);
            digits = tmp;
        }

        /// <summary>
        /// Converts the value to a proper decimal string representation.
        /// </summary>
        public override string ToString()
        {
            char[] digitString = new char[digits.Length];
            for (int i = 0; i < digits.Length; i++)
                digitString[i] = (char)(digits[i] + '0');

            switch (decimalPoint)
            {
                // Simplest case - nothing after the decimal point,
                // and last real digit is non-zero, eg value=35
                case 0:
                    return new string(digitString);
                // Fairly simple case - nothing after the decimal
                // point, but some 0s to add, eg value=350
                case < 0:
                    return new string(digitString) +
                           new string('0', -decimalPoint);
            }

            // Nothing before the decimal point, eg 0.035
            if (decimalPoint >= digitString.Length)
            {
                return "0." +
                       new string('0', decimalPoint - digitString.Length) +
                       new string(digitString);
            }

            // Most complicated case - part of the string comes
            // before the decimal point, part comes after it,
            // eg 3.5
            return new string(digitString, 0,
                       digitString.Length - decimalPoint) +
                   "." +
                   new string(digitString,
                       digitString.Length - decimalPoint,
                       decimalPoint);
        }
    }
}