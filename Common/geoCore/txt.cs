using gds;
using oasis;
using System;
using System.Collections.Generic;
using Clipper2Lib;

namespace geoCoreLib;

public class GCTxt : GCElement
{
    private string name;
    private Point64 point;
    private GCStrans trans;
    private int width;
    private int presentation;

    public GCTxt()
    {
        pGCTxt();
    }

    private void pGCTxt()
    {
        trans = new GCStrans();
    }

    public GCTxt(int l, int d, Point64 p, string s)
    {
        pGCTxt(l, d, p, s);
    }

    private void pGCTxt(int l, int d, Point64 p, string s)
    {
        name = s;
        point = p;
        layer_nr = l;
        datatype_nr = d;
        trans = new GCStrans();
        presentation = 0;
        width = 10;
    }

    public override void minimum(ref Point64 pos)
    {
        pMinimum(ref pos);
    }

    private void pMinimum(ref Point64 pos)
    {
        if (point.X < pos.X)
        {
            pos.X = point.X;
        }
        if (point.Y < pos.Y)
        {
            pos.Y = point.Y;
        }
    }

    public override void maximum(ref Point64 pos)
    {
        pMaximum(ref pos);
    }

    private void pMaximum(ref Point64 pos)
    {
        if (point.X > pos.X)
        {
            pos.X = point.X;
        }
        if (point.Y > pos.Y)
        {
            pos.Y = point.Y;
        }
    }

    public override void moveSelect(Point64 pos)
    {
        pMoveSelect(pos);
    }

    private void pMoveSelect(Point64 pos)
    {
        point = select switch
        {
            true => new (point.X + pos.X, point.Y + pos.Y),
            _ => point
        };
    }

    public override void move(Point64 pos)
    {
        pMove(pos);
    }

    private void pMove(Point64 pos)
    {
        point = new (point.X + pos.X, point.Y + pos.Y);
    }

    public override void setMirrorx()
    {
        pSetMirrorx();
    }

    private void pSetMirrorx()
    {
        trans.setMirror_x();
    }

    public override void clearMirrorx()
    {
        pClearMirrorx();
    }

    private void pClearMirrorx()
    {
        trans.clearMirror_x();
    }

    public override void toggleMirrorx()
    {
        pToggleMirrorx();
    }

    private void pToggleMirrorx()
    {
        trans.toggleMirror_x();
    }

    public override void rotate(double angle)
    {
        pRotate(angle);
    }

    private void pRotate(double angle)
    {
        trans.rotate(-angle);
    }

    public override void resize(double factor)
    {
        pResize(factor);
    }

    private void pResize(double factor)
    {
        point.X = (int)(point.X * factor);
        point.Y = (int)(point.Y * factor);
        width = (int)(width * factor);
    }

    public override void scale(double scale)
    {
        pScale(scale);
    }

    private void pScale(double scale)
    {
        trans.scale(scale);
    }

    public override bool isText()
    {
        return pIsText();
    }

    private bool pIsText()
    {
        return true;
    }

    public override string getName()
    {
        return pGetName();
    }

    private string pGetName()
    {
        return name;
    }

    public override void saveGDS(gdsWriter gw)
    {
        pSaveGDS(gw);
    }

    private void pSaveGDS(gdsWriter gw)
    {
        //text
        gw.bw.Write((ushort)4);
        gw.bw.Write(gdsValues.sTEXT);
        //layer
        gw.bw.Write((ushort)6);
        gw.bw.Write(gdsValues.sLAYER);
        gw.bw.Write((short)layer_nr);
        //datatype
        gw.bw.Write((ushort)6);
        gw.bw.Write(gdsValues.sDATATYPE);
        gw.bw.Write((short)datatype_nr);
        //presentation
        gw.bw.Write((ushort)6);
        gw.bw.Write(gdsValues.sPRESENTATION);
        gw.bw.Write((short)presentation);
        //width
        gw.bw.Write((ushort)8);
        gw.bw.Write(gdsValues.sWIDTH);
        gw.bw.Write(width);
        int strans_ = 0;
        switch (trans.mirror_x)
        {
            case true:
                strans_ |= 0x8000;
                break;
        }
        //STRANS
        gw.bw.Write((ushort)6);
        gw.bw.Write(gdsValues.sSTRANS);
        gw.bw.Write((short)strans_);
        //mag
        gw.bw.Write((ushort)12);
        gw.bw.Write(gdsValues.sMAG);
        gw.write8ByteReal(trans.mag);
        //angle
        gw.bw.Write((ushort)12);
        gw.bw.Write(gdsValues.sANGLE);
        switch (trans.mirror_x)
        {
            case true when trans.angle != 0:
                gw.write8ByteReal(360 - trans.angle);
                break;
            default:
                gw.write8ByteReal(trans.angle);
                break;
        }
        //xy
        const int val = 1 * 2 * 4 + 4;
        gw.bw.Write((ushort)val);
        gw.bw.Write(gdsValues.sXY);
        gw.bw.Write((int)point.X);
        gw.bw.Write((int)point.Y);
        gw.writeString(name, 0x19);
        // endel
        gw.bw.Write((ushort)4);
        gw.bw.Write(gdsValues.sENDEL);
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

        byte info_byte = 64;  //explicid text
        if (layer_nr != ow.modal.textlayer)
        {
            info_byte += 1;
        }
        if (datatype_nr != ow.modal.texttype)
        {
            info_byte += 2;
        }
        if (point.X != ow.modal.text_x)
        {
            info_byte += 16;
        }
        if (point.Y != ow.modal.text_y)
        {
            info_byte += 8;
        }
        ow.writeUnsignedInteger(19);
        ow.writeRaw(info_byte);
        ow.writeString(name);
        switch (info_byte & 1)
        {
            case > 0:
                ow.modal.textlayer = layer_nr;
                ow.writeUnsignedInteger((uint)ow.modal.textlayer);
                break;
        }
        switch (info_byte & 2)
        {
            case > 0:
                ow.modal.texttype = datatype_nr;
                ow.writeUnsignedInteger((uint)datatype_nr);
                break;
        }
        switch (info_byte & 16)
        {
            case > 0:
                ow.modal.text_x = (int)point.X;
                ow.writeSignedInteger(ow.modal.text_x);
                break;
        }
        switch (info_byte & 8)
        {
            case > 0:
                ow.modal.text_y = (int)point.Y;
                ow.writeSignedInteger(ow.modal.text_y);
                break;
        }
        if (presentation != GCSetup.defaultTextPresentation)
        {
            //o->error->addItem("Text presentation can not be saved in OASIS files.", 4);
        }
        switch (Math.Abs(trans.mag - 1))
        {
            case > double.Epsilon:
                //o->error->addItem("Scaled text can not be saved in OASIS files.", 4);
                break;
        }
        if (trans.angle != 0)
        {
            //o->error->addItem("Rotated text can not be saved in OASIS files.", 4);
        }
    }

    public override List<GCPolygon> convertToPolygons()
    {
        return pConvertToPolygons();
    }

    private List<GCPolygon> pConvertToPolygons()
    {
        List<GCPolygon> ret = new();
        Path64 points = Helper.initedPath64(5);
        points[0] = new (point.X - 1, point.Y - 1);
        points[1] = new (point.X - 1, point.Y + 1);
        points[2] = new (point.X + 1, point.Y + 1);
        points[3] = new (point.X + 1, point.Y - 1);
        points[4] = new (points[0]);
        ret.Add(new GCPolygon(points, layer_nr, datatype_nr));

        ret[0].text = true;
        ret[0].name = name;

        return ret;
    }
}