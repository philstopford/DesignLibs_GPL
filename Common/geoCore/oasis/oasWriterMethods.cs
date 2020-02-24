using geoLib;
using geoWrangler;
using System;

namespace oasis
{
    public partial class oasWriter
    {
        public void writeUnsignedInteger(uint i)
        {
            while (i > 127)
            {
                uint h = i & (127);
                i = i >> 7;
                h += 128;
                bw.Write((byte)(h));
            }
            bw.Write((byte)(i));
        }

        public void writeString(string text)
        {
            writeUnsignedInteger((uint)text.Length);
            bw.Write(text.ToCharArray());
        }

        public void writeSignedInteger(int i)
        {
            bool sig = false;
            if (i < 0)
            {
                sig = true;
            }
            i = Math.Abs(i);
            int h = (i & (63)) << 1;
            i = i >> 6;
            if (sig)
            {
                h++;
            }
            while (i != 0)
            {
                h += 128;
                bw.Write((byte)(h));
                h = i & (127);
                i = i >> 7;
            }
            bw.Write((byte)(h));
        }

        public void writeRaw(byte byte_)
        {
            bw.Write(byte_);
        }

        public void writeRaw(byte[] bytes)
        {
            bw.Write(bytes);
        }

        public void writeReal(double d)
        {
            if (Math.Abs(d) >= 0.5 && Math.Abs(Math.Floor(d + 0.5) - d) < 1e-6 && Math.Abs(d) < Double.MaxValue)
            {
                //  whole number (negative or positive)
                if (d < 0.0)
                {
                    writeUnsignedInteger(1);
                    writeUnsignedInteger((uint)Math.Floor(-d + 0.5));
                }
                else
                {
                    writeUnsignedInteger(0);
                    writeUnsignedInteger((uint)Math.Floor(d + 0.5));
                }
            }
            else
            {
                writeUnsignedInteger(7);
                byte[] a = new byte[8];
                a = BitConverter.GetBytes(d);
                bw.Write(a);
            }
        }

        public void writePointArray(GeoLibPoint[] p, bool excludeImplicid)
        {
            bool type0 = true;
            bool type1 = true;
            bool type2 = true;
            bool type3 = true;
            if (p.Length % 2 == 0)
            {
                type0 = false;
                type1 = false;
            }
            if (p.Length < 4)
            {
                type0 = false;
                type1 = false;
            }
            for (int i = 0; i < p.Length - 1; i++)
            {
                GeoLibPointF pd = GeoWrangler.distanceBetweenPoints_point(p[i + 1], p[i]);
                if ((pd.Y != 0) && (i % 2 != 0))
                {
                    type1 = false;
                }

                if ((pd.X != 0) && (i % 2 != 0))
                {
                    type0 = false;
                }

                if ((pd.X != 0) && (i % 2 != 1))
                {
                    type1 = false;
                }

                if ((pd.Y != 0) && (i % 2 != 1))
                {
                    type0 = false;
                }

                if ((pd.X != 0) && (pd.Y != 0))
                {
                    type2 = false;
                }

                if ((pd.X != 0) && (pd.Y != 0) && (pd.X != pd.Y) && (pd.X != -pd.Y))
                {
                    type3 = false;
                }
            }
            if (type0)
            {
                writeUnsignedInteger(0);
                int count = 0;
                if (excludeImplicid)
                {
                    count = p.Length - 3;
                }
                else
                {
                    count = p.Length - 1;
                }
                writeUnsignedInteger((uint)count);
                GeoLibPoint last = new GeoLibPoint(p[0]);
                for (int i = 1; i <= count; i++)
                {
                    GeoLibPointF h = GeoWrangler.distanceBetweenPoints_point(p[i], last);
                    write1Delta(h, i % 2 == 0);
                    last = new GeoLibPoint(p[i]);
                }
            }
            else if (type1)
            {
                writeUnsignedInteger(1);
                int count = 0;
                if (excludeImplicid)
                {
                    count = p.Length - 3;
                }
                else
                {
                    count = p.Length - 1;
                }
                writeUnsignedInteger((uint)count);
                GeoLibPoint last = new GeoLibPoint(p[0]);
                for (int i = 1; i <= count; i++)
                {
                    GeoLibPointF h = GeoWrangler.distanceBetweenPoints_point(p[i], last);
                    write1Delta(h, i % 2 == 1);
                    last = new GeoLibPoint(p[i]);
                }
            }
            else if (type2)
            {
                writeUnsignedInteger(2);
                int count = 0;
                if (excludeImplicid)
                {
                    count = p.Length - 2;
                }
                else
                {
                    count = p.Length - 1;
                }
                writeUnsignedInteger((uint)count);
                GeoLibPoint last = new GeoLibPoint(p[0]);
                for (int i = 1; i <= count; i++)
                {
                    GeoLibPointF h = GeoWrangler.distanceBetweenPoints_point(p[i], last);
                    write2Delta(h);
                    last = new GeoLibPoint(p[i]);
                }
            }
            else if (type3)
            {
                writeUnsignedInteger(3);
                int count = 0;
                if (excludeImplicid)
                {
                    count = p.Length - 2;
                }
                else
                {
                    count = p.Length - 1;
                }
                writeUnsignedInteger((uint)count);
                GeoLibPoint last = new GeoLibPoint(p[0]);
                for (int i = 1; i <= count; i++)
                {
                    GeoLibPointF h = GeoWrangler.distanceBetweenPoints_point(p[i], last);
                    write3Delta(h);
                    last = new GeoLibPoint(p[i]);
                }
            }
            else
            {
                writeUnsignedInteger(4);
                int count = 0;
                if (excludeImplicid)
                {
                    count = p.Length - 2;
                }
                else
                {
                    count = p.Length - 1;
                }
                writeUnsignedInteger((uint)count);
                GeoLibPoint last = new GeoLibPoint(p[0]);
                for (int i = 1; i <= count; i++)
                {
                    GeoLibPointF h = GeoWrangler.distanceBetweenPoints_point(p[i], last);
                    writeGDelta(h);
                    last = new GeoLibPoint(p[i]);
                }
            }
        }

        void write1Delta(GeoLibPointF p, bool dir)
        {
            int w;
            if (dir)
            {
                w = (int)p.Y;
            }
            else
            {
                w = (int)p.X;
            }
            if (w > 0)
            {
                w = (w << 1);
            }
            else
            {
                w = ((-w) << 1) + 1;
            }
            writeUnsignedInteger((uint)w);
        }

        void write2Delta(GeoLibPointF p)
        {
            int w;
            if (p.X == 0)
            {
                if (p.Y > 0)
                {
                    w = ((int)p.Y << 2) + 1;
                }
                else
                {
                    w = ((int)(-p.Y) << 2) + 3;
                }
            }
            else
            {
                if (p.X > 0)
                {
                    w = ((int)p.X << 2) + 0;
                }
                else
                {
                    w = ((int)(-p.X) << 2) + 2;
                }
            }
            writeUnsignedInteger((uint)w);
        }

        void write3Delta(GeoLibPointF p)
        {
            int w;
            if (p.X == 0)
            {
                if (p.Y > 0)
                {
                    w = ((int)p.Y << 3) + 1;
                }
                else
                {
                    w = ((int)-p.Y << 3) + 3;
                }
            }
            else if (p.Y == 0)
            {
                if (p.X > 0)
                {
                    w = ((int)p.X << 3) + 0;
                }
                else
                {
                    w = ((int)-p.X << 3) + 2;
                }
            }
            else if (p.Y == p.X)
            {
                if (p.Y > 0)
                {
                    w = ((int)p.Y << 3) + 4;
                }
                else
                {
                    w = ((int)-p.Y << 3) + 6;
                }
            }
            else
            {
                if (p.X > 0)
                {
                    w = ((int)p.X << 3) + 7;
                }
                else
                {
                    w = ((int)-p.X << 3) + 5;
                }
            }
            writeUnsignedInteger((uint)w);
        }

        void writeGDelta(GeoLibPointF p)
        {
            if ((p.Y == 0))
            {
                int i = 0;
                if (p.X >= 0)
                {
                    i = 0;
                    i = ((int)p.X << 4) + 2 * i;
                }
                else
                {
                    i = 2;
                    i = ((int)-p.X << 4) + 2 * i;
                }
                writeUnsignedInteger((uint)i);
            }
            else if (p.X == 0)
            {
                int i = 0;
                if (p.Y >= 0)
                {
                    i = 1;
                    i = ((int)p.Y << 4) + 2 * i;
                }
                else
                {
                    i = 3;
                    i = ((int)-p.Y << 4) + 2 * i;
                }
                writeUnsignedInteger((uint)i);
            }
            else if (p.Y == p.X)
            {
                int i = 0;
                if (p.X >= 0)
                {
                    i = 4;
                    i = ((int)p.X << 4) + 2 * i;
                }
                else
                {
                    i = 6;
                    i = ((int)-p.X << 4) + 2 * i;
                }
                writeUnsignedInteger((uint)i);
            }
            else if (p.Y == -p.X)
            {
                int i = 0;
                if (p.X >= 0)
                {
                    i = 7;
                    i = ((int)p.X << 4) + 2 * i;
                }
                else
                {
                    i = 5;
                    i = ((int)-p.X << 4) + 2 * i;
                }
                writeUnsignedInteger((uint)i);
            }
            else
            {
                int i = 0, k = 0, j = 0;
                if (p.X <= 0)
                {
                    i = 2;
                    k = ((int)-p.X << 2) + i + 1;
                }
                else
                {
                    k = ((int)p.X << 2) + i + 1;
                }
                writeUnsignedInteger((uint)k);
                if (p.Y <= 0)
                {
                    j = 1;
                    k = ((int)-p.Y << 1) + j;
                }
                else
                {
                    k = ((int)p.Y << 1) + j;
                }
                writeUnsignedInteger((uint)k);
            }
        }

        public void writeCtrapezoid(int layerNum, int type, int x, int y, int w, int h, int d)
        {
            if (!modal.absoluteMode)
            {
                modal.absoluteMode = true;
            }
            byte info_byte = 0;  //write pointlist;
            if (layerNum != modal.layer)
            {
                info_byte += 1;
            }
            if (d != modal.datatype)
            {
                info_byte += 2;
            }
            if (x != modal.geometry_x)
            {
                info_byte += 16;
            }
            if (y != modal.geometry_y)
            {
                info_byte += 8;
            }
            if (type != modal.ctrapezoid_type)
            {
                info_byte += 128;
            }
            if ((type != 20) && (type != 21))
            {
                if (w != modal.geometry_w)
                {
                    info_byte += 64;
                }
            }
            if ((type != 16) && (type != 17) && (type != 18) && (type != 19) && (type != 22) && (type != 23) && (type != 25))
            {
                if (h != modal.geometry_h)
                {
                    info_byte += 32;
                }
            }
            writeUnsignedInteger(26);
            writeRaw(info_byte);
            if ((info_byte & 1) > 0)
            {
                modal.layer = layerNum;
                writeUnsignedInteger((uint)modal.layer);
            }
            if ((info_byte & 2) > 0)
            {
                modal.datatype = d;
                writeUnsignedInteger((uint)d);
            }
            if ((info_byte & 128) > 0)
            {
                modal.ctrapezoid_type = type;
                writeUnsignedInteger((uint)modal.ctrapezoid_type);
            }
            if ((info_byte & 64) > 0)
            {
                modal.geometry_w = w;
                writeUnsignedInteger((uint)modal.geometry_w);
            }
            if ((info_byte & 32) > 0)
            {
                modal.geometry_h = h;
                writeUnsignedInteger((uint)modal.geometry_h);
            }
            if ((info_byte & 16) > 0)
            {
                modal.geometry_x = x;
                writeSignedInteger(modal.geometry_x);
            }
            if ((info_byte & 8) > 0)
            {
                modal.geometry_y = y;
                writeSignedInteger(modal.geometry_y);
            }
        }

        public void writeTrapezoid(int layerNum, int type, int x, int y, int w, int h, int da, int db, int d)
        {
            if (!modal.absoluteMode)
            {
                modal.absoluteMode = true;
            }
            byte info_byte = 0;  //write pointlist;
            if (layerNum != modal.layer)
            {
                info_byte += 1;
            }
            if (d != modal.datatype)
            {
                info_byte += 2;
            }
            if (x != modal.geometry_x)
            {
                info_byte += 16;
            }
            if (y != modal.geometry_y)
            {
                info_byte += 8;
            }
            if (type == 1)
            {
                info_byte += 128;
            }
            if (w != modal.geometry_w)
            {
                info_byte += 64;
            }
            if (h != modal.geometry_h)
            {
                info_byte += 32;
            }
            if (da == 0)
            {
                writeUnsignedInteger(25);
            }
            else if (db == 0)
            {
                writeUnsignedInteger(24);
            }
            else
            {
                writeUnsignedInteger(23);
            }
            writeRaw(info_byte);
            if ((info_byte & 1) > 0)
            {
                modal.layer = layerNum;
                writeUnsignedInteger((uint)modal.layer);
            }
            if ((info_byte & 2) > 0)
            {
                modal.datatype = d;
                writeUnsignedInteger((uint)d);
            }
            if ((info_byte & 64) > 0)
            {
                modal.geometry_w = w;
                writeUnsignedInteger((uint)modal.geometry_w);
            }
            if ((info_byte & 32) > 0)
            {
                modal.geometry_h = h;
                writeUnsignedInteger((uint)modal.geometry_h);
            }
            if (da == 0)
            {
                write1Delta(new GeoLibPointF(db, 0), false);
            }
            else if (db == 0)
            {
                write1Delta(new GeoLibPointF(da, 0), false);
            }
            else
            {
                write1Delta(new GeoLibPointF(da, 0), false);
                write1Delta(new GeoLibPointF(db, 0), false);
            }
            if ((info_byte & 16) > 0)
            {
                modal.geometry_x = x;
                writeSignedInteger(modal.geometry_x);
            }
            if ((info_byte & 8) > 0)
            {
                modal.geometry_y = y;
                writeSignedInteger(modal.geometry_y);
            }
        }
    }
}
