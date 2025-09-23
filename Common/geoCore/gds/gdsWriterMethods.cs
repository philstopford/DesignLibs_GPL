using System;

namespace gds;

public partial class gdsWriter
{

    public void write8ByteReal(double d)
    {
        pWrite8ByteReal(d);
    }

    private void pWrite8ByteReal(double d)
    {
        byte[] b = new byte[8];

        // Handle zero as special case
        if (d == 0.0)
        {
            bw.Write(b); // All zeros
            return;
        }

        // Extract sign and work with absolute value
        bool negative = d < 0;
        if (negative)
        {
            d = -d;
        }

        // Find the exponent such that mantissa is in range [1/16, 1)
        // We want: 1/16 <= d * 16^(-exponent) < 1
        // This means: log16(d) - 1 < exponent <= log16(d) + 4
        double log16d = Math.Log(d) / Math.Log(16.0);
        int exponent = (int)Math.Floor(log16d + 4);
        
        // Ensure mantissa is in proper range [1/16, 1)
        double mantissa = d / Math.Pow(16.0, exponent);
        while (mantissa >= 1.0)
        {
            exponent++;
            mantissa = d / Math.Pow(16.0, exponent);
        }
        while (mantissa < 1.0/16.0 && exponent > -64)
        {
            exponent--;
            mantissa = d / Math.Pow(16.0, exponent);
        }

        // Convert to excess-64 format
        int excessExponent = exponent + 64;
        if (excessExponent < 0) excessExponent = 0;
        if (excessExponent > 127) excessExponent = 127;

        // Set sign and exponent in first byte
        b[0] = (byte)(excessExponent & 0x7f);
        if (negative)
        {
            b[0] |= 0x80;
        }

        // Convert mantissa to 56-bit integer (7 bytes)
        // Binary point is to the left of bit 8, so multiply by 2^56
        ulong mantissaBits = (ulong)(mantissa * (1UL << 56) + 0.5);
        
        // Store mantissa in bytes 1-7 (big-endian)
        for (int i = 7; i > 0; i--)
        {
            b[i] = (byte)(mantissaBits & 0xff);
            mantissaBits >>= 8;
        }

        bw.Write(b);
    }

    public void writeString(string s, int type)
    {
        int len = s.Length;
        bool add = len % 2 == 1;
        if (add)
        {
            bw.Write((ushort)(len + 5));
        }
        else
        {
            bw.Write((ushort)(len + 4));
        }

        bw.Write((byte)type);
        bw.Write((byte)6);
        bw.Write(s.ToCharArray());
        if (add)
        {
            bw.Write((byte)0);
        }
    }
}