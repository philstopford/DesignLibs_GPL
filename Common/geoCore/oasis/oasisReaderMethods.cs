using geoLib;
using System;
using System.Linq;
using System.Text;

namespace oasis;

internal partial class oasReader
{
    private int readUnsignedInteger()
    {
        byte help;
        int result = 0;
        int pos = 0;
        do
        {
            help = readRaw();
            int h = help & 127;
            result += h << pos;
            pos += 7;
        }
        while (help >= 128);

        return pos switch
        {
            > 35 => throw new Exception("Integer with more than 32 Bits."),
            _ => result
        };
    }

    private byte readRaw()
    {
        switch (zLibUsed)
        {
            case true:
                return zLibReadRaw();
            default:
            {
                byte ret = br.ReadByte();
                return ret;
            }
        }
    }

    private void readProperty()
    {
        switch (readUnsignedInteger())
        {
            case 0:
                readUnsignedInteger();
                break;
            case 1:
                readUnsignedInteger();
                break;
            case 2:
                readUnsignedInteger();
                break;
            case 3:
                readUnsignedInteger();
                break;
            case 4:
                readUnsignedInteger();
                readUnsignedInteger();
                break;
            case 5:
                readUnsignedInteger();
                readUnsignedInteger();
                break;
            case 6:
            {
                //union {
                //float d;
                byte[] a = new byte[4];
                //};
                for (uint z = 0; z < 4; z++)
                {
                    a[z] = readRaw();
                }
            }
                break;
            case 7:
            {
                //union {
                //double d;
                byte[] a = new byte[8];
                //};
                for (uint z = 0; z < 8; z++)
                {
                    a[z] = readRaw();
                }
            }
                break;
            case 8:
                readUnsignedInteger();
                break;
            case 9:
                readSignedInteger();
                break;
            default:
                readString();
                break;
        }
    }

    private double readReal()
    {
        int i = readUnsignedInteger();
        switch (i)
        {
            case 0:
                return readUnsignedInteger();
            case 1:
                return -readUnsignedInteger();
            case 2:
                return (double)1 / readUnsignedInteger();
            case 3:
                return -(double)1 / readUnsignedInteger();
            case 4:
                return (double)readUnsignedInteger() / readUnsignedInteger();
            case 5:
                return -(double)readUnsignedInteger() / readUnsignedInteger();
            case 6:
            {
                byte[] a = new byte[4];

                for (uint z = 0; z < 4; z++)
                {
                    a[z] = readRaw();
                }
                double d = BitConverter.ToDouble(a, 0);
                return d;
            }
            case 7:
            {
                byte[] a = new byte[8];

                for (uint z = 0; z < 8; z++)
                {
                    a[z] = readRaw();
                }
                double d = BitConverter.ToDouble(a, 0);
                return d;
            }
            default:
                throw new Exception("Unknown real Format.");
        }
    }

    private string readString()
    {
        int items = readUnsignedInteger();
        uint i;
        string s1 = "";
        for (i = 0; i < items; i++)
        {
            byte help = readRaw();
            if (help != 0)
            {
                string s = Encoding.UTF8.GetString(new [] { help });
                s1 += s;
            }
        }
        return s1;
    }

    private int readSignedInteger()
    {
        int pos = 0;
        byte help = readRaw();
        bool sig = false;
        int h = help & 127;
        sig = (h % 2) switch
        {
            1 => true,
            _ => sig
        };
        int result = h >> 1;
        pos += 6;
        while (help >= 128)
        {
            help = readRaw();
            h = help & 127;
            result += h << pos;
            pos += 7;
        }
        switch (pos)
        {
            case > 34:
                throw new Exception("Integer with more then 32 Bits.");
        }

        return sig switch
        {
            true => -result,
            _ => result
        };
    }

    private GeoLibPoint read1Delta(bool dir)
    {
        // dir true-> vertical
        int i = readUnsignedInteger();

        switch (i % 2)
        {
            case 0:
                i >>= 1;
                break;
            case 1:
                i = -(i >> 1);
                break;
        }

        return dir switch
        {
            true => new GeoLibPoint(0, i),
            _ => new GeoLibPoint(i, 0)
        };
    }

    private GeoLibPoint readGDelta()
    {
        int i = readUnsignedInteger();
        int k;
        switch (i % 2)
        {
            case 0: //form 0
                i >>= 1;
                k = i % 8;
                i >>= 3;
                switch (k)
                {
                    case 0:
                        return new GeoLibPoint(i, 0);
                    case 1:
                        return new GeoLibPoint(0, i);
                    case 2:
                        return new GeoLibPoint(-i, 0);
                    case 3:
                        return new GeoLibPoint(0, -i);
                    case 4:
                        return new GeoLibPoint(i, i);
                    case 5:
                        return new GeoLibPoint(-i, i);
                    case 6:
                        return new GeoLibPoint(-i, -i);
                    case 7:
                        return new GeoLibPoint(i, -i);
                }
                break;
            case 1: //form 1
                i >>= 1;
                switch (i % 2)
                {
                    case 0:
                        i >>= 1;
                        break;
                    default:
                        i = -(i >> 1);
                        break;
                }
                k = readUnsignedInteger();
                switch (k % 2)
                {
                    case 0:
                        k >>= 1;
                        break;
                    default:
                        k = -(k >> 1);
                        break;
                }
                return new GeoLibPoint(i, k);
        }
        return new GeoLibPoint(0, 0);
    }

    private GeoLibPoint read2Delta()
    {
        int i = readUnsignedInteger();
        switch (i % 4)
        {
            case 0:
                i >>= 2;
                return new GeoLibPoint(i, 0);
            case 1:
                i >>= 2;
                return new GeoLibPoint(0, i);
            case 2:
                i = -(i >> 2);
                return new GeoLibPoint(i, 0);
            case 3:
                i = -(i >> 2);
                return new GeoLibPoint(0, i);
        }
        return new GeoLibPoint(0, 0);
    }

    private GeoLibPoint read3Delta()
    {
        int i = readUnsignedInteger();
        switch (i % 8)
        {
            case 0:
                i >>= 3;
                return new GeoLibPoint(i, 0);
            case 1:
                i >>= 3;
                return new GeoLibPoint(0, i);
            case 2:
                i >>= 3;
                return new GeoLibPoint(-i, 0);
            case 3:
                i >>= 3;
                return new GeoLibPoint(0, -i);
            case 4:
                i >>= 3;
                return new GeoLibPoint(i, i);
            case 5:
                i >>= 3;
                return new GeoLibPoint(-i, i);
            case 6:
                i >>= 3;
                return new GeoLibPoint(-i, -i);
            case 7:
                i >>= 3;
                return new GeoLibPoint(i, -i);
        }
        return new GeoLibPoint(0, 0);
    }

    private void readRepetition()
    {
        int i = readUnsignedInteger();
        switch (i)
        {
            case 0: break;
            case 1:
                modal.x_dimension = readUnsignedInteger() + 2;
                modal.y_dimension = readUnsignedInteger() + 2;
                modal.x_space = readUnsignedInteger();
                modal.y_space = readUnsignedInteger();
                break;
            case 2:
                modal.x_dimension = readUnsignedInteger() + 2;
                modal.x_space = readUnsignedInteger();
                break;
            case 3:
                modal.y_dimension = readUnsignedInteger() + 2;
                modal.y_space = readUnsignedInteger();
                break;
            case 4:
            {
                int k = readUnsignedInteger() + 2;
                int j = 0;
                modal.repArray.Clear();
                modal.repArray.Add(new GeoLibPoint(0, 0));
                for (int z = 1; z < k; z++)
                {
                    j += readUnsignedInteger();
                    modal.repArray.Add(new GeoLibPoint(j, 0));
                }
            }
                break;
            case 5:
            {
                int k = readUnsignedInteger() + 2;
                int j = 0;
                int grid = readUnsignedInteger();
                modal.repArray.Clear();
                modal.repArray.Add(new GeoLibPoint(0, 0));
                for (int z = 1; z < k; z++)
                {
                    j += readUnsignedInteger() * grid;
                    modal.repArray.Add(new GeoLibPoint(j, 0));
                }
            }
                break;
            case 6:
            {
                int k = readUnsignedInteger() + 2;
                int j = 0;
                modal.repArray.Clear();
                modal.repArray.Add(new GeoLibPoint(0, 0));
                for (int z = 1; z < k; z++)
                {
                    j += readUnsignedInteger();
                    modal.repArray.Add(new GeoLibPoint(0, j));
                }
            }
                break;
            case 7:
            {
                int k = readUnsignedInteger() + 2;
                int j = 0;
                int grid = readUnsignedInteger();
                modal.repArray.Clear();
                modal.repArray.Add(new GeoLibPoint(0, 0));
                for (int z = 1; z < k; z++)
                {
                    j += readUnsignedInteger() * grid;
                    modal.repArray.Add(new GeoLibPoint(0, j));
                }
            }
                break;
            case 8:
            {
                int anz1 = readUnsignedInteger() + 2;
                int anz2 = readUnsignedInteger() + 2;
                GeoLibPoint p1 = readGDelta();
                GeoLibPoint p2 = readGDelta();
                modal.repArray.Clear();
                for (int x1 = 0; x1 < anz1; x1++)
                {
                    for (int x2 = 0; x2 < anz2; x2++)
                    {
                        modal.repArray.Add(new GeoLibPoint());
                    }
                }
                for (int x1 = 0; x1 < anz1; x1++)
                {
                    for (int x2 = 0; x2 < anz2; x2++)
                    {
                        modal.repArray[x1 + x2 * anz1] = new GeoLibPoint(x1 * p1.X + x2 * p2.X, x1 * p1.Y + x2 * p2.Y);
                    }
                }
            }
                break;
            case 9:
            {
                int anz1 = readUnsignedInteger() + 2;
                GeoLibPoint p1 = readGDelta();
                modal.repArray.Clear();
                for (int x1 = 0; x1 < anz1; x1++)
                {
                    modal.repArray.Add(new GeoLibPoint(x1 * p1.X, x1 * p1.Y));
                }
            }
                break;
            case 10:
            {
                int k = readUnsignedInteger() + 2;
                GeoLibPoint j = new(0, 0);
                modal.repArray.Clear();
                modal.repArray.Add(new GeoLibPoint(0, 0));
                for (int z = 1; z < k; z++)
                {
                    j.Offset(readGDelta());
                    modal.repArray.Add(new GeoLibPoint(j));
                }
            }
                break;
            case 11:
            {
                int k = readUnsignedInteger() + 2;
                GeoLibPoint j = new(0, 0);
                int grid = readUnsignedInteger();
                modal.repArray.Clear();
                modal.repArray.Add(new GeoLibPoint(0, 0));
                for (int z = 1; z < k; z++)
                {
                    GeoLibPoint tmpPoint = readGDelta();
                    j.Offset(new GeoLibPoint(tmpPoint.X * grid, tmpPoint.Y * grid));
                    modal.repArray.Add(new GeoLibPoint(j));
                }
            }
                break;
            default: throw new Exception("Repetition unknown.");
        }
        if (i != 0)
        {
            modal.repetition = i;
        }
    }

    private void readExtension()
    {
        int i = readUnsignedInteger();
        switch (i & 3)
        {
            case 3:
                modal.path_start_extension_value = readSignedInteger();
                modal.path_start_extension = 4;
                break;
        }

        modal.path_start_extension = (i & 3) switch
        {
            2 => 2,
            _ => (i & 3) switch
            {
                1 => 0,
                _ => modal.path_start_extension
            }
        };
        switch (i & 12)
        {
            case 12:
                modal.path_end_extension_value = readSignedInteger();
                modal.path_end_extension = 4;
                break;
        }

        modal.path_end_extension = (i & 12) switch
        {
            8 => 2,
            _ => (i & 12) switch
            {
                4 => 0,
                _ => modal.path_end_extension
            }
        };
    }

    private void readPointList(bool addImplecid)
    {
        int type = readUnsignedInteger();
        int count = readUnsignedInteger();
        GeoLibPoint[] pointlist = new GeoLibPoint[count + 1];
        pointlist[0] = new GeoLibPoint(0, 0);
        GeoLibPoint last = new(0, 0);
        int i;
        bool dir = true;
        switch (type)
        {
            case 0:
                dir = false;
                for (i = 1; i <= count; i++)
                {
                    last.Offset(read1Delta(dir));
                    pointlist[i] = new GeoLibPoint(last);
                    dir = !dir;
                }
                switch (addImplecid)
                {
                    case true:
                        Array.Resize(ref pointlist, pointlist.Length + 2);
                        pointlist[^2] = new GeoLibPoint(0, pointlist[^3].Y);
                        pointlist[^1] = new GeoLibPoint(pointlist[0]);
                        break;
                }
                break;
            case 1:
                for (i = 1; i <= count; i++)
                {
                    last.Offset(read1Delta(dir));
                    pointlist[i] = new GeoLibPoint(last);
                    dir = !dir;
                }
                switch (addImplecid)
                {
                    case true:
                        Array.Resize(ref pointlist, pointlist.Length + 2);
                        pointlist[^2] = new GeoLibPoint(pointlist[^3].X, 0);
                        pointlist[^1] = new GeoLibPoint(pointlist[0]);
                        break;
                }
                break;
            case 2:
                for (i = 1; i <= count; i++)
                {
                    last.Offset(read2Delta());
                    pointlist[i] = new GeoLibPoint(last);
                }
                switch (addImplecid)
                {
                    case true:
                        Array.Resize(ref pointlist, pointlist.Length + 1);
                        pointlist[^1] = new GeoLibPoint(pointlist[0]);
                        break;
                }
                break;
            case 3:
                for (i = 1; i <= count; i++)
                {
                    last.Offset(read3Delta());
                    pointlist[i] = new GeoLibPoint(last);
                }
                switch (addImplecid)
                {
                    case true:
                        Array.Resize(ref pointlist, pointlist.Length + 1);
                        pointlist[^1] = new GeoLibPoint(pointlist[0]);
                        break;
                }
                break;
            case 4:
                for (i = 1; i <= count; i++)
                {
                    last.Offset(readGDelta());
                    pointlist[i] = new GeoLibPoint(last);
                }
                switch (addImplecid)
                {
                    case true:
                        Array.Resize(ref pointlist, pointlist.Length + 1);
                        pointlist[^1] = new GeoLibPoint(pointlist[0]);
                        break;
                }
                break;
            case 5:
            {
                GeoLibPoint l = new(0, 0);
                for (i = 1; i <= count; i++)
                {
                    last.Offset(readGDelta());
                    l.Offset(last);
                    pointlist[i] = new GeoLibPoint(l);
                }
                switch (addImplecid)
                {
                    case true:
                        Array.Resize(ref pointlist, pointlist.Length + 1);
                        pointlist[^1] = new GeoLibPoint(pointlist[0]);
                        break;
                }
            }
                break;
        }
        modal.polygon_point_list = pointlist.ToList();
    }

    private void zLibInit(uint before, uint after)
    {
        zLibUsed = true;
        byte[] data = br.ReadBytes((int)(after - before));
        zLibOut = utility.Utils.decompress(data);
        zlibOutPos = 0;
    }

    private byte zLibReadRaw()
    {
        byte get = zLibOut[zlibOutPos];
        zlibOutPos++;
        if (zlibOutPos == zLibOut.Length)
        {
            zLibUsed = false;
        }
        return get;
    }
}