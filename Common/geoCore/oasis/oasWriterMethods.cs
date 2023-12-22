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
        switch (sig)
        {
            case true:
                h++;
                break;
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
        switch (Math.Abs(d))
        {
            case >= 0.5 when Math.Abs(Math.Floor(d + 0.5) - d) < 1e-6 && Math.Abs(d) < double.MaxValue:
            {
                switch (d)
                {
                    //  whole number (negative or positive)
                    case < 0.0:
                        writeUnsignedInteger(1);
                        writeUnsignedInteger((uint)Math.Floor(-d + 0.5));
                        break;
                    default:
                        writeUnsignedInteger(0);
                        writeUnsignedInteger((uint)Math.Floor(d + 0.5));
                        break;
                }

                break;
            }
            default:
            {
                writeUnsignedInteger(7);
                byte[] a = BitConverter.GetBytes(d);
                bw.Write(a);
                break;
            }
        }
    }

    public void writePointArray(Path64 p, bool excludeImplicid)
    {
        bool type0 = true;
        bool type1 = true;
        bool type2 = true;
        bool type3 = true;
        switch (p.Count % 2)
        {
            case 0:
                type0 = false;
                type1 = false;
                break;
        }
        switch (p.Count)
        {
            case < 4:
                type0 = false;
                type1 = false;
                break;
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
        switch (type0)
        {
            case true:
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
                    last = new (p[i]);
                }

                break;
            }
            default:
            {
                switch (type1)
                {
                    case true:
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
                            last = new (p[i]);
                        }

                        break;
                    }
                    default:
                    {
                        switch (type2)
                        {
                            case true:
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
                                    last = new (p[i]);
                                }

                                break;
                            }
                            default:
                            {
                                switch (type3)
                                {
                                    case true:
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
                                            last = new (p[i]);
                                        }

                                        break;
                                    }
                                    default:
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
                                            last = new (p[i]);
                                        }

                                        break;
                                    }
                                }

                                break;
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }
    }

    private void write1Delta(PointD p, bool dir)
    {
        int w = dir switch
        {
            true => (int) p.y,
            _ => (int) p.x
        };
        switch (w)
        {
            case > 0:
                w <<= 1;
                break;
            default:
                w = (-w << 1) + 1;
                break;
        }
        writeUnsignedInteger((uint)w);
    }

    private void write2Delta(PointD p)
    {
        int w = p.x switch
        {
            0 when p.y > 0 => ((int) p.y << 2) + 1,
            0 => ((int) -p.y << 2) + 3,
            > 0 => ((int) p.x << 2) + 0,
            _ => ((int) -p.x << 2) + 2
        };
        writeUnsignedInteger((uint)w);
    }

    private void write3Delta(PointD p)
    {
        int w;
        switch (p.x)
        {
            case 0 when p.y > 0:
                w = ((int)p.y << 3) + 1;
                break;
            case 0:
                w = ((int)-p.y << 3) + 3;
                break;
            default:
            {
                switch (p.y)
                {
                    case 0 when p.x > 0:
                        w = ((int)p.x << 3) + 0;
                        break;
                    case 0:
                        w = ((int)-p.x << 3) + 2;
                        break;
                    default:
                    {
                        switch (Math.Abs(p.y - p.x))
                        {
                            case <= double.Epsilon when p.y > 0:
                                w = ((int)p.y << 3) + 4;
                                break;
                            case <= double.Epsilon:
                                w = ((int)-p.y << 3) + 6;
                                break;
                            default:
                            {
                                w = p.x switch
                                {
                                    > 0 => ((int) p.x << 3) + 7,
                                    _ => ((int) -p.x << 3) + 5
                                };

                                break;
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }
        writeUnsignedInteger((uint)w);
    }

    public void writeGDelta(PointD p)
    {
        switch (p.y)
        {
            case 0:
            {
                int i;
                switch (p.x)
                {
                    case >= 0:
                        i = ((int)p.x << 4);
                        break;
                    default:
                        i = 2;
                        i = ((int)-p.x << 4) + 2 * i;
                        break;
                }
                writeUnsignedInteger((uint)i);
                break;
            }
            default:
            {
                switch (p.x)
                {
                    case 0:
                    {
                        int i;
                        switch (p.y)
                        {
                            case >= 0:
                                i = ((int)p.y << 4) + 2;
                                break;
                            default:
                                i = 3;
                                i = ((int)-p.y << 4) + 2 * i;
                                break;
                        }
                        writeUnsignedInteger((uint)i);
                        break;
                    }
                    default:
                    {
                        switch (Math.Abs(p.y - p.x))
                        {
                            case <= double.Epsilon:
                            {
                                int i;
                                switch (p.x)
                                {
                                    case >= 0:
                                        i = 4;
                                        i = ((int)p.x << 4) + 2 * i;
                                        break;
                                    default:
                                        i = 6;
                                        i = ((int)-p.x << 4) + 2 * i;
                                        break;
                                }
                                writeUnsignedInteger((uint)i);
                                break;
                            }
                            default:
                            {
                                switch (Math.Abs(p.y - -p.x))
                                {
                                    case <= double.Epsilon:
                                    {
                                        int i;
                                        switch (p.x)
                                        {
                                            case >= 0:
                                                i = 7;
                                                i = ((int)p.x << 4) + 2 * i;
                                                break;
                                            default:
                                                i = 5;
                                                i = ((int)-p.x << 4) + 2 * i;
                                                break;
                                        }
                                        writeUnsignedInteger((uint)i);
                                        break;
                                    }
                                    default:
                                    {
                                        int k;
                                        switch (p.x)
                                        {
                                            case <= 0:
                                            {
                                                const int i = 2;
                                                k = ((int)-p.x << 2) + i + 1;
                                                break;
                                            }
                                            default:
                                                k = ((int)p.x << 2) + 1;
                                                break;
                                        }
                                        writeUnsignedInteger((uint)k);
                                        switch (p.y)
                                        {
                                            case <= 0:
                                                int j = 1;
                                                k = ((int)-p.y << 1) + j;
                                                break;
                                            default:
                                                k = (int)p.y << 1;
                                                break;
                                        }
                                        writeUnsignedInteger((uint)k);
                                        break;
                                    }
                                }

                                break;
                            }
                        }

                        break;
                    }
                }

                break;
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
        switch (info_byte & 1)
        {
            case > 0:
                modal.layer = layerNum;
                writeUnsignedInteger((uint)modal.layer);
                break;
        }
        switch (info_byte & 2)
        {
            case > 0:
                modal.datatype = d;
                writeUnsignedInteger((uint)d);
                break;
        }
        switch (info_byte & 128)
        {
            case > 0:
                modal.ctrapezoid_type = type;
                writeUnsignedInteger((uint)modal.ctrapezoid_type);
                break;
        }
        switch (info_byte & 64)
        {
            case > 0:
                modal.geometry_w = w;
                writeUnsignedInteger((uint)modal.geometry_w);
                break;
        }
        switch (info_byte & 32)
        {
            case > 0:
                modal.geometry_h = h;
                writeUnsignedInteger((uint)modal.geometry_h);
                break;
        }
        switch (info_byte & 16)
        {
            case > 0:
                modal.geometry_x = x;
                writeSignedInteger(modal.geometry_x);
                break;
        }
        switch (info_byte & 8)
        {
            case > 0:
                modal.geometry_y = y;
                writeSignedInteger(modal.geometry_y);
                break;
        }
    }

    public void writeTrapezoid(int layerNum, int type, int x, int y, int w, int h, int da, int db, int d)
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
        switch (type)
        {
            case 1:
                info_byte += 128;
                break;
        }
        if (w != modal.geometry_w)
        {
            info_byte += 64;
        }
        if (h != modal.geometry_h)
        {
            info_byte += 32;
        }
        switch (da)
        {
            case 0:
                writeUnsignedInteger(25);
                break;
            default:
            {
                switch (db)
                {
                    case 0:
                        writeUnsignedInteger(24);
                        break;
                    default:
                        writeUnsignedInteger(23);
                        break;
                }

                break;
            }
        }
        writeRaw(info_byte);
        switch (info_byte & 1)
        {
            case > 0:
                modal.layer = layerNum;
                writeUnsignedInteger((uint)modal.layer);
                break;
        }
        switch (info_byte & 2)
        {
            case > 0:
                modal.datatype = d;
                writeUnsignedInteger((uint)d);
                break;
        }
        switch (info_byte & 64)
        {
            case > 0:
                modal.geometry_w = w;
                writeUnsignedInteger((uint)modal.geometry_w);
                break;
        }
        switch (info_byte & 32)
        {
            case > 0:
                modal.geometry_h = h;
                writeUnsignedInteger((uint)modal.geometry_h);
                break;
        }
        switch (da)
        {
            case 0:
                write1Delta(new (db, 0), false);
                break;
            default:
            {
                switch (db)
                {
                    case 0:
                        write1Delta(new (da, 0), false);
                        break;
                    default:
                        write1Delta(new (da, 0), false);
                        write1Delta(new (db, 0), false);
                        break;
                }

                break;
            }
        }
        switch (info_byte & 16)
        {
            case > 0:
                modal.geometry_x = x;
                writeSignedInteger(modal.geometry_x);
                break;
        }
        switch (info_byte & 8)
        {
            case > 0:
                modal.geometry_y = y;
                writeSignedInteger(modal.geometry_y);
                break;
        }
    }
}