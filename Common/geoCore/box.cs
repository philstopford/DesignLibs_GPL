using gds;
using geoLib;
using oasis;
using System.Collections.Generic;
using Clipper2Lib;

namespace geoCoreLib;

public class GCBox : GCElement
{
    private GeoLibRectangle rect;

    public GCBox()
    {
    }

    public GCBox(GeoLibRectangle r, int layer, int datatype)
    {
        pGCBox(r, layer, datatype);
    }

    private void pGCBox(GeoLibRectangle r, int layer, int datatype)
    {
        rect = r;
        layer_nr = layer;
        datatype_nr = datatype;
    }

    public GCBox(int x_, int y_, int b_, int h_, int layer, int datatype)
    {
        pGCBox(x_, y_, b_, h_, layer, datatype);
    }

    private void pGCBox(int x_, int y_, int b_, int h_, int layer, int datatype)
    {
        rect = new GeoLibRectangle(x_, y_, b_, h_);
        pGCBox(rect, layer, datatype);
    }

    public override void clean()
    {
        // Nothing needed since the rectangle primitive does everything.
    }

    public override bool correct()
    {
        return pCorrect();
    }

    private bool pCorrect()
    {
        if (rect.Left == rect.Right)
        {
            return false;
        }
        return rect.Top != rect.Bottom;
    }

    public override void maximum(Point64 p)
    {
        pMaximum(p);
    }

    private void pMaximum(Point64 p)
    {
        p.X = rect.Right;
        p.Y = rect.Top;
    }

    public override void minimum(Point64 p)
    {
        pMinimum(p);
    }

    private void pMinimum(Point64 p)
    {
        p.X = rect.Left;
        p.Y = rect.Bottom;
    }

    public override void move(Point64 pos)
    {
        pMove(pos);
    }

    private void pMove(Point64 pos)
    {
        rect.Offset(pos);
    }

    public override void moveSelect(Point64 p)
    {
        pMoveSelect(p);
    }

    private void pMoveSelect(Point64 p)
    {
        switch (select)
        {
            case true:
                rect.Offset(p);
                break;
        }
    }

    public override void resize(double size)
    {
        pResize(size);
    }

    private void pResize(double size)
    {
        rect = new GeoLibRectangle(
            (int)(rect.X * size),
            (int)(rect.Y * size),
            (int)(rect.Width * size),
            (int)(rect.Height * size)
        );
    }

    public override List<GCPolygon> convertToPolygons()
    {
        return pConvertToPolygons();
    }

    private List<GCPolygon> pConvertToPolygons()
    {
        List<GCPolygon> ret = new();
        Path64 points = Helper.initedPath64(5);
        points[0] = new (rect.Left, rect.Top);
        points[1] = new (rect.Right, rect.Top);
        points[2] = new (rect.Right, rect.Bottom);
        points[3] = new (rect.Left, rect.Bottom);
        points[4] = new (rect.Left, rect.Top);
        ret.Add(new GCPolygon(points, layer_nr, datatype_nr));
        return ret;
    }

    public override void saveGDS(gdsWriter gw)
    {
        pSaveGDS(gw);
    }

    private void pSaveGDS(gdsWriter gw)
    {
        // box
        gw.bw.Write((ushort)4);
        gw.bw.Write((byte)0x2D);
        gw.bw.Write((byte)0);
        //layer
        gw.bw.Write((ushort)6);
        gw.bw.Write((byte)0x0D);
        gw.bw.Write((byte)2);
        gw.bw.Write((short)layer_nr);
        //boxtype
        gw.bw.Write((ushort)6);
        gw.bw.Write((byte)0x2E);
        gw.bw.Write((byte)2);
        gw.bw.Write((short)datatype_nr);
        //xy 
        gw.bw.Write((ushort)(5 * 2 * 4 + 4));
        gw.bw.Write((byte)0x10);
        gw.bw.Write((byte)3);
        gw.bw.Write(rect.Left);
        gw.bw.Write(rect.Top);
        gw.bw.Write(rect.Right);
        gw.bw.Write(rect.Top);
        gw.bw.Write(rect.Right);
        gw.bw.Write(rect.Bottom);
        gw.bw.Write(rect.Left);
        gw.bw.Write(rect.Bottom);
        gw.bw.Write(rect.Left);
        gw.bw.Write(rect.Top);
        // endel
        gw.bw.Write((ushort)4);
        gw.bw.Write((byte)0x11);
        gw.bw.Write((byte)0);
    }

    public override void saveOASIS(oasWriter ow)
    {
        pSaveOASIS(ow);
    }

    private void pSaveOASIS(oasWriter ow)
    {
        ow.modal.absoluteMode = ow.modal.absoluteMode switch
        {
            false => true,
            _ => ow.modal.absoluteMode
        };

        byte info_byte = 0;
        if (layer_nr != ow.modal.layer)
        {
            info_byte += 1;
        }
        if (datatype_nr != ow.modal.datatype)
        {
            info_byte += 2;
        }
        if (rect.Left != ow.modal.geometry_x)
        {
            info_byte += 16;
        }
        if (rect.Bottom != ow.modal.geometry_y)
        {
            info_byte += 8;
        }
        if (rect.Height != ow.modal.geometry_h)
        {
            info_byte += 32;
        }
        if (rect.Width != ow.modal.geometry_w)
        {
            info_byte += 64;
        }
        if (rect.Width == rect.Height)
        {
            info_byte += 128;
            info_byte = (byte)(info_byte & (255 - 32));
        }
        ow.writeUnsignedInteger(20);
        ow.writeRaw(info_byte);
        switch (info_byte & 1)
        {
            case > 0:
                ow.modal.layer = layer_nr;
                ow.writeUnsignedInteger((uint)ow.modal.layer);
                break;
        }
        switch (info_byte & 2)
        {
            case > 0:
                ow.modal.datatype = datatype_nr;
                ow.writeUnsignedInteger((uint)datatype_nr);
                break;
        }
        switch (info_byte & 64)
        {
            case > 0:
            {
                ow.modal.geometry_w = rect.Width;
                ow.writeUnsignedInteger((uint)ow.modal.geometry_w);

                break;
            }
        }

        switch (info_byte & 128)
        {
            case > 0:
                ow.modal.geometry_h = ow.modal.geometry_w;
                break;
            default:
                ow.modal.geometry_h = rect.Height;
                break;
        }

        switch (info_byte & 32)
        {
            case > 0:
                ow.writeUnsignedInteger((uint)ow.modal.geometry_h);
                break;
        }
        switch (info_byte & 16)
        {
            case > 0:
                ow.modal.geometry_x = rect.Left; //rect.X;
                ow.writeSignedInteger(ow.modal.geometry_x);
                break;
        }
        switch (info_byte & 8)
        {
            case > 0:
                ow.modal.geometry_y = rect.Bottom; // rect.Y;
                ow.writeSignedInteger(ow.modal.geometry_y);
                break;
        }
    }
}