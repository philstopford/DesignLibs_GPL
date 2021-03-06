﻿using System;

namespace gds
{
    public partial class gdsWriter
    {

        public void write8ByteReal(double d)
        {
            pWrite8ByteReal(d);
        }

        void pWrite8ByteReal(double d)
        {
            byte[] b = new byte[8];

            b[0] = 0;
            if (d < 0)
            {
                b[0] = 0x80;
                d = -d;
            }

            //  compute the next power of 16 that that value will fit in
            int e = 0;
            if (d < 1e-77 /*~16^-64*/)
            {
                d = 0;
            }
            else
            {
                double lg16 = Math.Log(d) / Math.Log(16.0);
                e = (int)(Math.Ceiling(Math.Log(d) / Math.Log(16.0)));
                if (e == lg16)
                {
                    ++e;
                }
            }

            d /= Math.Pow(16.0, e - 14);

            b[0] |= (byte)((e + 64) & 0x7f);

            UInt64 m = (UInt64)(d + 0.5);
            for (int i = 7; i > 0; --i)
            {
                b[i] = (byte)(m & 0xff);
                m >>= 8;
            }

            bw.Write(b);

            /*
            return;
            utility.FrExp.FRexpResult r = utility.FrExp.calculate(d);


            double mantissa = r.mantissa;
            int exponent = r.exponent;

            int i = Math.Abs(exponent) % 4;
            int exptwo = exponent;
            exponent = exponent / 4;
            if (i != 0)
            {
                exponent++;
                if (exptwo > 0)
                {
                    mantissa /= Math.Pow(2, 4 - i);
                }
                else
                {
                    mantissa /= Math.Pow(2, i);
                    exponent--;
                }
            }

            byte help = (byte)(exponent + 64);
            if (0 > mantissa)
            {
                help |= 0x80;
                mantissa *= -1;
            }
            bw.Write(help);
            for (int k = 0; k < 7; k++)
            {
                mantissa *= 256;
                help = (byte)mantissa;
                bw.Write(help);
                mantissa -= help;
            }
            */
        }

        public void writeString(string s, int type)
        {
            int len = s.Length;
            bool add = false;
            if (len % 2 == 1)
            {
                add = true;
                //len++;
            }
            if (add)
            {
                bw.Write((UInt16)(len + 5));
            }
            else
            {
                bw.Write((UInt16)(len + 4));
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
}
