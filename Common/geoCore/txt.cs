using gds;
using geoLib;
using oasis;
using System;

namespace geoCoreLib
{
    public class GCTxt : GCElement
    {
        string name;
        GeoLibPoint point;
        GCStrans trans;
        int width;
        int presentation;

        public GCTxt()
        {
            pGCTxt();
        }

        void pGCTxt()
        {
            trans = new GCStrans();
        }

        public GCTxt(int l, int d, GeoLibPoint p, string s)
        {
            pGCTxt(l, d, p, s);
        }

        void pGCTxt(int l, int d, GeoLibPoint p, string s)
        {
            name = s;
            point = p;
            layer_nr = l;
            datatype_nr = d;
            trans = new GCStrans();
            presentation = 0;
            width = -10;
        }

        public override void minimum(GeoLibPoint pos)
        {
            pMinimum(pos);
        }

        void pMinimum(GeoLibPoint pos)
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

        public override void maximum(GeoLibPoint pos)
        {
            pMaximum(pos);
        }

        void pMaximum(GeoLibPoint pos)
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

        public override void moveSelect(GeoLibPoint pos)
        {
            pMoveSelect(pos);
        }

        void pMoveSelect(GeoLibPoint pos)
        {
            if (select)
            {
                point = new GeoLibPoint(point.X + pos.X, point.Y + pos.Y);
            }
        }

        public override void move(GeoLibPoint pos)
        {
            pMove(pos);
        }

        void pMove(GeoLibPoint pos)
        {
            point = new GeoLibPoint(point.X + pos.X, point.Y + pos.Y);
        }

        public override void setMirrorx()
        {
            pSetMirrorx();
        }

        void pSetMirrorx()
        {
            trans.setMirror_x();
        }

        public override void clearMirrorx()
        {
            pClearMirrorx();
        }

        void pClearMirrorx()
        {
            trans.clearMirror_x();
        }

        public override void toggleMirrorx()
        {
            pToggleMirrorx();
        }

        void pToggleMirrorx()
        {
            trans.toggleMirror_x();
        }

        public override void rotate(double angle)
        {
            pRotate(angle);
        }

        void pRotate(double angle)
        {
            trans.rotate(-angle);
        }

        public override void scale(double scale)
        {
            pScale(scale);
        }

        void pScale(double scale)
        {
            trans.scale(scale);
        }

        public override bool isText()
        {
            return pIsText();
        }

        bool pIsText()
        {
            return true;
        }

        public override string getName()
        {
            return pGetName();
        }

        string pGetName()
        {
            return name;
        }

        public override void saveGDS(gdsWriter gw)
        {
            pSaveGDS(gw);
        }

        void pSaveGDS(gdsWriter gw)
        {
            //text
            gw.bw.Write((UInt16)4);
            gw.bw.Write((byte)0x0C);
            gw.bw.Write((byte)0);
            //layer
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x0D);
            gw.bw.Write((byte)2);
            gw.bw.Write((Int16)layer_nr);
            //datatype
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x16);
            gw.bw.Write((byte)2);
            gw.bw.Write((Int16)(datatype_nr));
            //presentation
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x17);
            gw.bw.Write((byte)1);
            gw.bw.Write((Int16)(presentation));
            //width
            gw.bw.Write((UInt16)8);
            gw.bw.Write((byte)0x0F);
            gw.bw.Write((byte)3);
            gw.bw.Write(width);
            Int32 strans_ = 0;
            if (trans.mirror_x)
            {
                strans_ |= 0x8000;
            };
            //STRANS
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x1A);
            gw.bw.Write((byte)1);
            gw.bw.Write((Int16)strans_);
            //mag
            gw.bw.Write((UInt16)12);
            gw.bw.Write((byte)0x1B);
            gw.bw.Write((byte)5);
            gw.write8ByteReal(trans.mag);
            //angle
            gw.bw.Write((UInt16)12);
            gw.bw.Write((byte)0x1C);
            gw.bw.Write((byte)5);
            if ((trans.mirror_x) && (trans.angle != 0))
            {
                gw.write8ByteReal(360 - trans.angle);
            }
            else
            {
                gw.write8ByteReal(trans.angle);
            }
            //xy
            int val = (1 * 2 * 4) + 4;
            gw.bw.Write((UInt16)val);
            gw.bw.Write((byte)0x10);
            gw.bw.Write((byte)3);
            gw.bw.Write(point.X);
            gw.bw.Write(point.Y);
            gw.writeString(name, 0x19);
            // endel
            gw.bw.Write((UInt16)4);
            gw.bw.Write((byte)0x11);
            gw.bw.Write((byte)0);
        }

        public override void saveOASIS(oasWriter ow)
        {
            pSaveOASIS(ow);
        }

        void pSaveOASIS(oasWriter ow)
        {
            if (!ow.modal.absoluteMode)
            {
                ow.modal.absoluteMode = true;
            }

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
            if ((info_byte & 1) > 0)
            {
                ow.modal.textlayer = layer_nr;
                ow.writeUnsignedInteger((uint)ow.modal.textlayer);
            }
            if ((info_byte & 2) > 0)
            {
                ow.modal.texttype = datatype_nr;
                ow.writeUnsignedInteger((uint)datatype_nr);
            }
            if ((info_byte & 16) > 0)
            {
                ow.modal.text_x = point.X;
                ow.writeSignedInteger(ow.modal.text_x);
            }
            if ((info_byte & 8) > 0)
            {
                ow.modal.text_y = point.Y;
                ow.writeSignedInteger(ow.modal.text_y);
            }
            if (presentation != GCSetup.defaultTextPresentation)
            {
                //o->error->addItem("Text presentation can not be saved in OASIS files.", 4);
            }
            if (trans.mag != 1)
            {
                //o->error->addItem("Scaled text can not be saved in OASIS files.", 4);
            }
            if (trans.angle != 0)
            {
                //o->error->addItem("Rotated text can not be saved in OASIS files.", 4);
            }
        }

        public override GCPolygon convertToPolygon()
        {
            return pConvertToPolygon();
        }

        GCPolygon pConvertToPolygon()
        {
            GeoLibPoint[] points = new GeoLibPoint[5];
            points[0] = new GeoLibPoint(point.X - 1, point.Y - 1);
            points[1] = new GeoLibPoint(point.X - 1, point.Y + 1);
            points[2] = new GeoLibPoint(point.X + 1, point.Y + 1);
            points[3] = new GeoLibPoint(point.X + 1, point.Y - 1);
            points[4] = new GeoLibPoint(points[0]);
            return new GCPolygon(points, layer_nr, datatype_nr);
        }
    }
}
