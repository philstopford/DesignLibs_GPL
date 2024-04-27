using System;
using System.Text;
using Clipper2Lib;
using geoCoreLib;
using geoWrangler;

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
            if (help == 0)
            {
                continue;
            }

            string s = Encoding.UTF8.GetString(new [] { help });
            s1 += s;
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

    private Point64 read1Delta(bool dir)
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
            true => new Point64(0, i),
            _ => new Point64(i, 0)
        };
    }

    private Point64 readGDelta()
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
                        return new (i, 0);
                    case 1:
                        return new (0, i);
                    case 2:
                        return new (-i, 0);
                    case 3:
                        return new (0, -i);
                    case 4:
                        return new (i, i);
                    case 5:
                        return new (-i, i);
                    case 6:
                        return new (-i, -i);
                    case 7:
                        return new (i, -i);
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
                return new (i, k);
        }
        return new (0, 0);
    }

    private Point64 read2Delta()
    {
        int i = readUnsignedInteger();
        switch (i % 4)
        {
            case 0:
                i >>= 2;
                return new (i, 0);
            case 1:
                i >>= 2;
                return new (0, i);
            case 2:
                i = -(i >> 2);
                return new (i, 0);
            case 3:
                i = -(i >> 2);
                return new (0, i);
        }
        return new (0, 0);
    }

    private Point64 read3Delta()
    {
        int i = readUnsignedInteger();
        switch (i % 8)
        {
            case 0:
                i >>= 3;
                return new (i, 0);
            case 1:
                i >>= 3;
                return new (0, i);
            case 2:
                i >>= 3;
                return new (-i, 0);
            case 3:
                i >>= 3;
                return new (0, -i);
            case 4:
                i >>= 3;
                return new (i, i);
            case 5:
                i >>= 3;
                return new (-i, i);
            case 6:
                i >>= 3;
                return new (-i, -i);
            case 7:
                i >>= 3;
                return new (i, -i);
        }
        return new (0, 0);
    }

    private void readRepetition()
    {
        modal.repetition = new();
        int i = readUnsignedInteger();
        switch (i)
        {
            case 0:
                break;
            case 1:
                modal.repetition.type = Repetition.RepetitionType.Rectangular;
                modal.repetition.columns = readUnsignedInteger() + 2;
                modal.repetition.rows = readUnsignedInteger() + 2;
                modal.repetition.spacing.X = readUnsignedInteger();
                modal.repetition.spacing.Y = readUnsignedInteger();
                break;
            case 2:
                modal.repetition.type = Repetition.RepetitionType.Rectangular;
                modal.repetition.columns = readUnsignedInteger() + 2;
                modal.repetition.rows = 1;
                modal.repetition.spacing.X = readUnsignedInteger();
                modal.repetition.spacing.Y = 0;
                break;
            case 3:
                modal.repetition.type = Repetition.RepetitionType.Rectangular;
                modal.repetition.columns = 1;
                modal.repetition.rows = readUnsignedInteger() + 2;
                modal.repetition.spacing.X = 0;
                modal.repetition.spacing.Y = readUnsignedInteger();
                break;
            case 4:
            case 5:
            {
                modal.repetition.type = Repetition.RepetitionType.ExplicitX;
                int count = readUnsignedInteger() + 1;
                double grid_factor = 1.0; // REVIEW: scaling
                if (i == 5)
                {
                    grid_factor *= readUnsignedInteger();
                }

                while (count > 0)
                {
                    modal.repetition.coords.Add(grid_factor + readUnsignedInteger());
                    count--;
                }

                break;
            }
            case 6:
            case 7:
            {
                modal.repetition.type = Repetition.RepetitionType.ExplicitY;
                int count = readUnsignedInteger() + 1;
                int grid_factor = 1; // REVIEW: scaling
                if (i == 7)
                {
                    grid_factor *= readUnsignedInteger();
                }

                while (count > 0)
                {
                    modal.repetition.coords.Add(grid_factor + readUnsignedInteger());
                    count--;
                }

                break;
            }
            case 8:
                modal.repetition.type = Repetition.RepetitionType.Regular;
                modal.repetition.columns = readUnsignedInteger() + 2;
                modal.repetition.rows = readUnsignedInteger() + 2;
                modal.repetition.colVector = readGDelta();
                modal.repetition.rowVector = readGDelta();
                break;
            case 9:
                modal.repetition.type = Repetition.RepetitionType.Regular;
                modal.repetition.columns = readUnsignedInteger() + 2;
                modal.repetition.rows = 1;
                modal.repetition.colVector = readGDelta();
                modal.repetition.rowVector = new(modal.repetition.colVector.Y, modal.repetition.colVector.X);
                break;
            case 10:
            case 11:
            {
                modal.repetition.type = Repetition.RepetitionType.Explicit;
                int count = readUnsignedInteger() + 2;
                int grid_factor = 1;
                if (i == 11)
                {
                    grid_factor = readUnsignedInteger();
                }

                while (count > 0)
                {
                    Point64 p = readGDelta();
                    modal.repetition.offsets.Add(new (p.X * grid_factor, p.Y * grid_factor));
                    count--;
                }
                break;
            }
            default:
                throw new Exception("Repetition unknown.");
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
        Path64 pointlist = Helper.initedPath64(count + 1);
        pointlist[0] = new (0, 0);
        Point64 last = new(0, 0);
        Point64 oPt;
        int i;
        bool dir = true;
        switch (type)
        {
            case 0:
                dir = false;
                for (i = 1; i <= count; i++)
                {
                    oPt = read1Delta(dir);
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    pointlist[i] = new (last);
                    dir = !dir;
                }
                switch (addImplecid)
                {
                    case true:
                        /*
                        pointlist.Add(new());
                        pointlist.Add(new());
                        pointlist[^2] = new (pointlist[^3].X, 0);
                        pointlist[^1] = new (pointlist[0]);
                        */
                        pointlist.Add(new (pointlist[^1].X, 0));
                        pointlist.Add(new (pointlist[0]));
                        break;
                }
                break;
            case 1:
                for (i = 1; i <= count; i++)
                {
                    oPt = read1Delta(dir);
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    pointlist[i] = new (last);
                    dir = !dir;
                }
                switch (addImplecid)
                {
                    case true:
                        /*
                        pointlist.Add(new());
                        pointlist.Add(new());
                        pointlist[^2] = new (pointlist[^3].X, 0);
                        pointlist[^1] = new (pointlist[0]);
                        */
                        pointlist.Add(new (pointlist[^1].X, 0));
                        pointlist.Add(new (pointlist[0]));
                        break;
                }
                break;
            case 2:
                for (i = 1; i <= count; i++)
                {
                    oPt = read2Delta();
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    pointlist[i] = new (last);
                }
                switch (addImplecid)
                {
                    case true:
                        pointlist.Add(new (pointlist[0]));
                        break;
                }
                break;
            case 3:
                for (i = 1; i <= count; i++)
                {
                    oPt = read3Delta();
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    pointlist[i] = new (last);
                }
                switch (addImplecid)
                {
                    case true:
                        pointlist.Add(new (pointlist[0]));
                        break;
                }
                break;
            case 4:
                for (i = 1; i <= count; i++)
                {
                    oPt = readGDelta();
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    pointlist[i] = new (last);
                }
                switch (addImplecid)
                {
                    case true:
                        pointlist.Add(new (pointlist[0]));
                        break;
                }
                break;
            case 5:
            {
                Point64 l = new(0, 0);
                for (i = 1; i <= count; i++)
                {
                    oPt = readGDelta();
                    last = GeoWrangler.move(last, oPt.X, oPt.Y);
                    l = GeoWrangler.move(l, last.X, last.Y);
                    pointlist[i] = new (l);
                }
                switch (addImplecid)
                {
                    case true:
                        pointlist.Add(new (pointlist[0]));
                        break;
                }
            }
                break;
        }
        modal.polygon_point_list = new(pointlist);
    }

    private void zLibInit(uint before, uint after)
    {
        zLibUsed = true;
        Int64 br_pos = br.BaseStream.Position;
        byte[] data = br.ReadBytes((int)before);
        try
        {
            zLibOut = utility.Utils.decompress(data);
            if (zLibOut.Length != after)
            {
                throw new ("Incorrect number of bytes after decompression");
            }
        }
        catch (Exception e)
        {
            Console.WriteLine(e);
            throw;
        }
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