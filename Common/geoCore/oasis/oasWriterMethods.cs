using geoWrangler;
using System;
using Clipper2Lib;

namespace oasis;

public partial class oasWriter
{
    public void writeUnsignedInteger(uint i)
    {
        while (i > 127)
        {
            uint h = i & 127;
            i >>= 7;
            h += 128;
            bw.Write((byte)h);
        }
        bw.Write((byte)i);
    }

    public void writeString(string text)
    {
        writeUnsignedInteger((uint)text.Length);
        bw.Write(text.ToCharArray());
    }

    public void writeSignedInteger(int i)
    {
        bool sig = i < 0;
        i = Math.Abs(i);
        int h = (i & 63) << 1;
        i >>= 6;
        if (sig)
        {
            h++;
        }

        while (i != 0)
        {
            h += 128;
            bw.Write((byte)h);
            h = i & 127;
            i >>= 7;
        }
        bw.Write((byte)h);
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
        if (Math.Abs(d) >= 0.5 && (Math.Abs(Math.Floor(d + 0.5) - d) < 1e-6 && Math.Abs(d) < double.MaxValue))
        {
            if (d <
                //  whole number (negative or positive)
                0.0)
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
            byte[] a = BitConverter.GetBytes(d);
            bw.Write(a);
        }
    }

    public void writePointArray(Path64 p, bool excludeImplicid)
    {
        bool type0 = true;
        bool type1 = true;
        bool type2 = true;
        bool type3 = true;
        if (p.Count % 2 == 0)
        {
            type0 = false;
            type1 = false;
        }

        if (p.Count < 4)
        {
            type0 = false;
            type1 = false;
        }

        for (int i = 0; i < p.Count - 1; i++)
        {
            PointD pd = GeoWrangler.distanceBetweenPoints_point(p[i + 1], p[i]);
            if (pd.y != 0 && i % 2 != 0)
            {
                type1 = false;
            }

            if (pd.x != 0 && i % 2 != 0)
            {
                type0 = false;
            }

            if (pd.x != 0 && i % 2 != 1)
            {
                type1 = false;
            }

            if (pd.y != 0 && i % 2 != 1)
            {
                type0 = false;
            }

            if (pd.x != 0 && pd.y != 0)
            {
                type2 = false;
            }

            if (pd.x != 0 && pd.y != 0 && Math.Abs(pd.x - pd.y) > double.Epsilon && Math.Abs(pd.x - -pd.y) > double.Epsilon)
            {
                type3 = false;
            }
        }

        if (type0)
        {
            writeUnsignedInteger(0);
            int count = excludeImplicid switch
            {
                true => p.Count - 3,
                _ => p.Count - 1
            };
            writeUnsignedInteger((uint)count);
            Point64 last = new(p[0]);
            for (int i = 1; i <= count; i++)
            {
                PointD h = GeoWrangler.distanceBetweenPoints_point(p[i], last);
                write1Delta(h, i % 2 == 0);
                last = new(p[i]);
            }
        }
        else
        {
            if (type1)
            {
                writeUnsignedInteger(1);
                int count = excludeImplicid switch
                {
                    true => p.Count - 3,
                    _ => p.Count - 1
                };
                writeUnsignedInteger((uint)count);
                Point64 last = new(p[0]);
                for (int i = 1; i <= count; i++)
                {
                    PointD h = GeoWrangler.distanceBetweenPoints_point(p[i], last);
                    write1Delta(h, i % 2 == 1);
                    last = new(p[i]);
                }
            }
            else
            {
                if (type2)
                {
                    writeUnsignedInteger(2);
                    int count = excludeImplicid switch
                    {
                        true => p.Count - 2,
                        _ => p.Count - 1
                    };
                    writeUnsignedInteger((uint)count);
                    Point64 last = new(p[0]);
                    for (int i = 1; i <= count; i++)
                    {
                        PointD h = GeoWrangler.distanceBetweenPoints_point(p[i], last);
                        write2Delta(h);
                        last = new(p[i]);
                    }
                }
                else
                {
                    if (type3)
                    {
                        writeUnsignedInteger(3);
                        int count = excludeImplicid switch
                        {
                            true => p.Count - 2,
                            _ => p.Count - 1
                        };
                        writeUnsignedInteger((uint)count);
                        Point64 last = new(p[0]);
                        for (int i = 1; i <= count; i++)
                        {
                            PointD h = GeoWrangler.distanceBetweenPoints_point(p[i], last);
                            write3Delta(h);
                            last = new(p[i]);
                        }
                    }
                    else
                    {
                        writeUnsignedInteger(4);
                        int count = excludeImplicid switch
                        {
                            true => p.Count - 2,
                            _ => p.Count - 1
                        };
                        writeUnsignedInteger((uint)count);
                        Point64 last = new(p[0]);
                        for (int i = 1; i <= count; i++)
                        {
                            PointD h = GeoWrangler.distanceBetweenPoints_point(p[i], last);
                            writeGDelta(h);
                            last = new(p[i]);
                        }
                    }
                }
            }
        }
    }

    private void write1Delta(PointD p, bool dir)
    {
        int w;
        if (dir)
        {
            w = (int)p.y;
        }
        else
        {
            w = (int)p.x;
        }

        if (w > 0)
        {
            w <<= 1;
        }
        else
        {
            w = (-w << 1) + 1;
        }

        writeUnsignedInteger((uint)w);
    }

    private void write2Delta(PointD p)
    {
        int w;
        if (p.x == 0 && p.y > 0)
        {
            w = ((int)p.y << 2) + 1;
        }
        else if (p.x == 0)
        {
            w = ((int)-p.y << 2) + 3;
        }
        else if (p.x > 0)
        {
            w = ((int)p.x << 2) + 0;
        }
        else
        {
            w = ((int)-p.x << 2) + 2;
        }

        writeUnsignedInteger((uint)w);
    }

    private void write3Delta(PointD p)
    {
        int w;
        if (p.x == 0 && p.y > 0)
        {
            w = ((int)p.y << 3) + 1;
        }
        else if (p.x == 0)
        {
            w = ((int)-p.y << 3) + 3;
        }
        else
        {
            if (p.y == 0 && p.x > 0)
            {
                w = ((int)p.x << 3) + 0;
            }
            else if (p.y == 0)
            {
                w = ((int)-p.x << 3) + 2;
            }
            else
            {
                if (Math.Abs(p.y - p.x) <= double.Epsilon && p.y > 0)
                {
                    w = ((int)p.y << 3) + 4;
                }
                else if (Math.Abs(p.y - p.x) <= double.Epsilon)
                {
                    w = ((int)-p.y << 3) + 6;
                }
                else
                {
                    w = p.x switch
                    {
                        > 0 => ((int)p.x << 3) + 7,
                        _ => ((int)-p.x << 3) + 5
                    };
                }
            }
        }

        writeUnsignedInteger((uint)w);
    }

    public void writeGDelta(PointD p)
    {
        if (p.y == 0)
        {
            int i;
            if (p.x >= 0)
            {
                i = ((int)p.x << 4);
            }
            else
            {
                i = 2;
                i = ((int)-p.x << 4) + 2 * i;
            }

            writeUnsignedInteger((uint)i);
        }
        else
        {
            if (p.x == 0)
            {
                int i;
                if (p.y >= 0)
                {
                    i = ((int)p.y << 4) + 2;
                }
                else
                {
                    i = 3;
                    i = ((int)-p.y << 4) + 2 * i;
                }

                writeUnsignedInteger((uint)i);
            }
            else
            {
                if (Math.Abs(p.y - p.x) <= double.Epsilon)
                {
                    int i;
                    if (p.x >= 0)
                    {
                        i = 4;
                        i = ((int)p.x << 4) + 2 * i;
                    }
                    else
                    {
                        i = 6;
                        i = ((int)-p.x << 4) + 2 * i;
                    }

                    writeUnsignedInteger((uint)i);
                }
                else
                {
                    if (Math.Abs(p.y - -p.x) <= double.Epsilon)
                    {
                        int i;
                        if (p.x >= 0)
                        {
                            i = 7;
                            i = ((int)p.x << 4) + 2 * i;
                        }
                        else
                        {
                            i = 5;
                            i = ((int)-p.x << 4) + 2 * i;
                        }

                        writeUnsignedInteger((uint)i);
                    }
                    else
                    {
                        int k;
                        if (p.x <= 0)
                        {
                            const int i = 2;
                            k = ((int)-p.x << 2) + i + 1;
                        }
                        else
                        {
                            k = ((int)p.x << 2) + 1;
                        }

                        writeUnsignedInteger((uint)k);
                        if (p.y <= 0)
                        {
                            int j = 1;
                            k = ((int)-p.y << 1) + j;
                        }
                        else
                        {
                            k = (int)p.y << 1;
                        }

                        writeUnsignedInteger((uint)k);
                    }
                }
            }
        }
    }

    public void writeCtrapezoid(int layerNum, int type, int x, int y, int w, int h, int d)
    {
        modal.absoluteMode = modal.absoluteMode switch
        {
            false => true,
            _ => modal.absoluteMode
        };
        byte info_byte = 0;  //write point-list;
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
        if (type != 20 && type != 21)
        {
            if (w != modal.geometry_w)
            {
                info_byte += 64;
            }
        }
        if (type != 16 && type != 17 && type != 18 && type != 19 && type != 22 && type != 23 && type != 25)
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
        else
        {
            modal.absoluteMode = modal.absoluteMode;
        }

        byte info_byte = 0;  //write point-list;
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
            writeUnsignedInteger(oasValues.TRAPEZOID_B);
        }
        else
        {
            if (db == 0)
            {
                writeUnsignedInteger(oasValues.TRAPEZOID_A);
            }
            else
            {
                writeUnsignedInteger(oasValues.TRAPEZOID_AB);
            }
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
            write1Delta(new(db, 0), false);
        }
        else
        {
            if (db == 0)
            {
                write1Delta(new(da, 0), false);
            }
            else
            {
                write1Delta(new(da, 0), false);
                write1Delta(new(db, 0), false);
            }
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