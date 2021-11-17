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

        b[0] = 0;
        switch (d)
        {
            case < 0:
                b[0] = 0x80;
                d = -d;
                break;
        }

        //  compute the next power of 16 that that value will fit in
        int e = 0;
        switch (d)
        {
            /*~16^-64*/
            case < 1e-77:
                d = 0;
                break;
            default:
            {
                double lg16 = Math.Log(d) / Math.Log(16.0);
                e = (int)Math.Ceiling(Math.Log(d) / Math.Log(16.0));
                switch (Math.Abs(e - lg16))
                {
                    case <= double.Epsilon:
                        ++e;
                        break;
                }

                break;
            }
        }

        d /= Math.Pow(16.0, e - 14);

        b[0] |= (byte)((e + 64) & 0x7f);

        ulong m = (ulong)(d + 0.5);
        for (int i = 7; i > 0; --i)
        {
            b[i] = (byte)(m & 0xff);
            m >>= 8;
        }

        bw.Write(b);
    }

    public void writeString(string s, int type)
    {
        int len = s.Length;
        bool add = len % 2 == 1;
        switch (add)
        {
            case true:
                bw.Write((ushort)(len + 5));
                break;
            default:
                bw.Write((ushort)(len + 4));
                break;
        }
        bw.Write((byte)type);
        bw.Write((byte)6);
        bw.Write(s.ToCharArray());
        switch (add)
        {
            case true:
                bw.Write((byte)0);
                break;
        }
    }
}