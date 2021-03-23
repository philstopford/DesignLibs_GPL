using geoCoreLib;
using geoLib;
using MiscUtil.Conversion;
using MiscUtil.IO;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using utility;

namespace gds
{
    partial class gdsReader
    {
        public delegate void StatusUpdateUI(string text);
        public StatusUpdateUI statusUpdateUI { get; set; }
        public delegate void ProgressUpdateUI(double progress);
        public ProgressUpdateUI progressUpdateUI { get; set; }

        public List<string> error_msgs;
        public bool valid { get; set; }
        public bool abortLoad_ { get; set; }
        GCDrawingfield drawing_;
        string filename;
        EndianBinaryReader br;
        GCCell cell_;

        UInt16 reclen;
        byte record;
        byte type, help;

        Int32 items, int32, int32x, int32y;
        Int16 int16, help16;
        double double_;

        public struct modals
        {
            public Int32 beginExt { get; set; }
            public Int32 endExt { get; set; }
            public Int32 width { get; set; }
            public Int32 cap { get; set; }
            public Int32 elementmode { get; set; }
            public Int32 anzx { get; set; }
            public Int32 anzy { get; set; }
            public Int32 layer { get; set; }
            public Int32 datatype { get; set; }
            public Int32 textType { get; set; }
            public Int32 presentation { get; set; }
            public Int32 boxType { get; set; }
            public double angle { get; set; }
            public double mag { get; set; }
            public bool mirror_x { get; set; }
            public bool rotate { get; set; }
            public bool mag_ { get; set; }
            public string sname { get; set; }

            public GeoLibPoint[] point_array { get; set; }
        }

        modals modal;

        public void reset()
        {
            pReset();
        }

        void pReset()
        {
            drawing_.reset();
            modal.point_array = null;
            cell_ = null;
            resetModal();
            valid = false;
            error_msgs.Clear();
        }

        void resetModal()
        {
            modal.sname = "";

            modal.anzx = 0;
            modal.anzy = 0;
            modal.layer = 1;
            modal.datatype = 1;
            modal.textType = 0;
            modal.presentation = 0;
            modal.boxType = 0;

            modal.angle = 0;
            modal.mag = 1;
            modal.mirror_x = false;
            modal.rotate = false;
            modal.mag_ = false;
            modal.beginExt = 0;
            modal.endExt = 0;
            modal.width = 0;
            modal.cap = 0;
            modal.elementmode = 0;
        }

        public Dictionary<string, string> layerNames { get; set; }

        public gdsReader(String filename)
        {
            pGDSReader(filename);
        }

        void pGDSReader(String filename)
        {
            drawing_ = new GCDrawingfield(filename);
            this.filename = filename;
            error_msgs = new List<string>();
        }

        public bool load(ref GCDrawingfield drawing)
        {
            return pLoad(ref drawing);
        }

        bool pLoad(ref GCDrawingfield drawing)
        {
            drawing_ = drawing;
            layerNames = new Dictionary<string, string>();
            valid = true;
            resetModal();
            try
            {
                statusUpdateUI?.Invoke("Loading");
                progressUpdateUI?.Invoke(0);
                Stream s;
                s = File.OpenRead(filename);
                if (filename.ToLower().EndsWith("gz"))
                {
                    using (GZipStream gzs = new GZipStream(s, CompressionMode.Decompress))
                    {
                        MemoryStream ms = new MemoryStream();
                        gzs.CopyTo(ms);
                        ms.Seek(0, SeekOrigin.Begin);
                        br = new EndianBinaryReader(EndianBitConverter.Big, ms);
                    }
                }
                else
                {
                    br = new EndianBinaryReader(EndianBitConverter.Big, s);
                }
                System.Text.Encoding ascii = System.Text.Encoding.ASCII;

                List<GCCell> firstcellhelp = drawing_.cellList;

                cell_ = null;

                long pos = br.BaseStream.Position;
                long len = br.BaseStream.Length;

                long updateInterval = len / 100;

                double progress = 0;

                while (pos < len)
                {
                    if (pos % updateInterval == 0)
                    {
                        progressUpdateUI?.Invoke(progress);
                        progress += 0.01f;
                    }
                    reclen = br.ReadUInt16();
                    record = br.ReadByte();
                    type = br.ReadByte();

                    switch (type)
                    {
                        case 0:
                            break; // no data
                        case 1:
                            items = (reclen - 4) / 2; // bit array
                            break;
                        case 2:
                            items = (reclen - 4) / 2; // two byte signed int
                            break;
                        case 3:
                            items = (reclen - 4) / 4; // four byte signed int
                            break;
                        case 4:
                            items = (reclen - 4) / 4; // four byte real
                            break;
                        case 5:
                            items = (reclen - 4) / 8; // eight byte real
                            break;
                        case 6:
                            items = (reclen - 4); // ASC
                            break;

                    }

                    switch (record)
                    {
                        case 0: // header
                            help16 = br.ReadInt16();
                            break;
                        case 1: // bgnlib
                            drawing_.modyear = br.ReadInt16();
                            drawing_.modmonth = br.ReadInt16();
                            drawing_.modday = br.ReadInt16();
                            drawing_.modhour = br.ReadInt16();
                            drawing_.modmin = br.ReadInt16();
                            drawing_.modsec = br.ReadInt16();
                            // drawing mod time above
                            drawing_.accyear = br.ReadInt16();
                            drawing_.accmonth = br.ReadInt16();
                            drawing_.accday = br.ReadInt16();
                            drawing_.acchour = br.ReadInt16();
                            drawing_.accmin = br.ReadInt16();
                            drawing_.accsec = br.ReadInt16();
                            // drawing access time above
                            break;
                        case 2:  // LIBNAM
                            byte[] tmp_lbnbyte = br.ReadBytes(items);

                            char[] tmplbnchars = new char[ascii.GetCharCount(tmp_lbnbyte, 0, tmp_lbnbyte.Length)];
                            ascii.GetChars(tmp_lbnbyte, 0, tmp_lbnbyte.Length, tmplbnchars, 0);

                            drawing_.libname = new string(tmplbnchars); //  br.ReadString();
                            // Some files have null termination issues.
                            string[] lNTokens = drawing_.libname.Split(new char[] { '\0' });
                            drawing_.libname = lNTokens[0];
                            break;
                        case 3: // UNITS
                            double_ = read8ByteReal();
                            drawing_.userunits = double_;
                            double_ = read8ByteReal();
                            drawing_.databaseunits = double_;
                            break;
                        case 4: // ENDLIB
                            break;
                        case 5: // BGNSTR
                            cell_ = drawing_.addCell();

                            cell_.modyear = br.ReadInt16();
                            cell_.modmonth = br.ReadInt16();
                            cell_.modday = br.ReadInt16();
                            cell_.modhour = br.ReadInt16();
                            cell_.modmin = br.ReadInt16();
                            cell_.modsec = br.ReadInt16();
                            cell_.accyear = br.ReadInt16();
                            cell_.accmonth = br.ReadInt16();
                            cell_.accday = br.ReadInt16();
                            cell_.acchour = br.ReadInt16();
                            cell_.accmin = br.ReadInt16();
                            cell_.accsec = br.ReadInt16();

                            break;
                        case 6: // STRNAM
                            byte[] tmp_strnbyte = br.ReadBytes(items);

                            char[] tmpstrnchars = new char[ascii.GetCharCount(tmp_strnbyte, 0, tmp_strnbyte.Length)];
                            ascii.GetChars(tmp_strnbyte, 0, tmp_strnbyte.Length, tmpstrnchars, 0);

                            cell_.cellName = new string(tmpstrnchars);
                            // Some files have null termination issues.
                            string[] tokens = cell_.cellName.Split(new char[] { '\0' });
                            cell_.cellName = tokens[0];
                            break;
                        case 7: // ENDSTR
                            break;
                        case 8: // BONDRY -> Polygon
                            modal.elementmode = 110;
                            break;
                        case 9: //PATH 
                            modal.elementmode = 150;
                            modal.angle = 0;
                            modal.mag = 1;
                            modal.mirror_x = false;
                            modal.width = 0;
                            modal.cap = 0;
                            modal.endExt = 0;
                            modal.beginExt = 0;
                            break;
                        case 10: //SREF
                            modal.elementmode = 120;
                            modal.angle = 0;
                            modal.mag = 1;
                            modal.mirror_x = false;
                            modal.rotate = false;
                            modal.mag_ = false;
                            break;
                        case 11: //AREF
                            modal.elementmode = 130;
                            modal.angle = 0;
                            modal.mag = 1;
                            modal.mirror_x = false;
                            modal.rotate = false;
                            modal.mag_ = false;
                            break;
                        case 12: //TEXT
                            modal.elementmode = 140;
                            modal.angle = 0;
                            modal.mag = 1;
                            modal.width = 0;
                            modal.mirror_x = false;
                            modal.rotate = false;
                            modal.mag_ = false;
                            modal.presentation = 0;
                            break;
                        case 13: //LAYER
                            modal.layer = br.ReadInt16();
                            break;
                        case 14: //DTATYP
                            modal.datatype = br.ReadInt16();
                            try
                            {
                                layerNames.Add("L" + modal.layer + "D" + modal.datatype, "L" + modal.layer + "D" + modal.datatype);
                            }
                            catch (Exception)
                            {

                            }
                            break;
                        case 15: //WIDTH 
                            modal.width = br.ReadInt32();
                            if (modal.width < 0)
                            {
                                modal.width = -modal.width;
                            }
                            // Zero width is problematic for Variance and Quilt.
                            if (modal.width == 0)
                            {
                                modal.width = 10;
                            }
                            break;
                        case 16: //XY
                            modal.point_array = new GeoLibPoint[(items / 2)];
                            for (Int32 g = 0; g < modal.point_array.Length; g++)
                            {
                                int32x = br.ReadInt32();
                                int32y = br.ReadInt32();
                                modal.point_array[g] = new GeoLibPoint(int32x, int32y);
                            }
                            break;
                        case 17: //ENDEL
                            // Looks like some cases, we don't register the associated LD (e.g. layers with only text). Workaround for now until cause is understood.
                            try
                            {
                                layerNames.Add("L" + modal.layer + "D" + modal.datatype, "L" + modal.layer + "D" + modal.datatype);
                            }
                            catch (Exception)
                            {

                            }
                            switch (modal.elementmode)
                            {
                                case 100:
                                    addBox();
                                    break;
                                case 110:
                                    addPolygon();
                                    break;
                                case 120:
                                    addCellRef();
                                    break;
                                case 130:
                                    addCellRefArray();
                                    break;
                                case 140:
                                    addText();
                                    break;
                                case 150:
                                    addPath();
                                    break;
                            }
                            resetModal();
                            break;
                        case 18: //SNAME
                            byte[] tmp_snamebyte = br.ReadBytes(items);

                            char[] tmpnamechars = new char[ascii.GetCharCount(tmp_snamebyte, 0, tmp_snamebyte.Length)];
                            ascii.GetChars(tmp_snamebyte, 0, tmp_snamebyte.Length, tmpnamechars, 0);

                            modal.sname = new string(tmpnamechars);
                            // Some files have null termination issues.
                            string[] sNtokens = modal.sname.Split(new char[] { '\0' });
                            modal.sname = sNtokens[0];
                            break;
                        case 19: //COLROW
                            modal.anzx = br.ReadInt16();
                            modal.anzy = br.ReadInt16();
                            break;
                        /*case 20: //TXTNOD
							break;
						case 21: //NODE
							break;*/
                        case 22: //TXTTYP
                            modal.textType = br.ReadInt16();
                            modal.datatype = modal.textType; // we don't treat text differently.
                            break;
                        case 23: //PRSTTN
                            modal.presentation = br.ReadInt16() & 0x000F;
                            break;
                        /*case 24: //SPACNG
							break;*/
                        case 25: //STRING
                            byte[] tmp_strbyte = br.ReadBytes(items);

                            char[] tmpstrchars = new char[ascii.GetCharCount(tmp_strbyte, 0, tmp_strbyte.Length)];
                            ascii.GetChars(tmp_strbyte, 0, tmp_strbyte.Length, tmpstrchars, 0);
                            modal.sname = new string(tmpstrchars);
                            // Some files have null termination issues.
                            string[] sN2tokens = modal.sname.Split(new char[] { '\0' });
                            modal.sname = sN2tokens[0];
                            break;
                        case 26: //STRANS
                            int16 = br.ReadInt16();
                            if ((int16 & 0x8000) != 0)
                            {
                                modal.mirror_x = true;
                            }
                            else
                            {
                                modal.mirror_x = false;
                            }
                            if ((int16 & 0x0002) != 0)
                            {
                                modal.rotate = true;
                            }
                            else
                            {
                                modal.rotate = false;
                            }
                            if ((int16 & 0x0004) != 0)
                            {
                                modal.mag_ = true;
                            }
                            else
                            {
                                modal.mag_ = false;
                            };
                            break;
                        case 27: //MAG
                            modal.mag = read8ByteReal();
                            break;
                        case 28: //ANGLE 
                            modal.angle = read8ByteReal();
                            break;
                        /*case 29: //UINTEG
							break;
						case 30: //USTRNG
							break;
						case 31: //REFLIB
							s=readString(&ts,items);
							if (layout::debug){
								printf("REFLIB %s\n",s.latin1());
							}
							break;
						case 32: //FONTS 
							break;*/
                        case 33: //PTHTYP
                            modal.cap = br.ReadInt16();
                            break;
                        /*case 34: //GENRTS
							break;
						case 35: //ATRTBL
							break;
						case 36: //STPTBL
							break;
						case 37: //STRTYP
							break;
						case 38: //EFLAGS
							break;
						case 39: //ELKEY 
							break;
						case 40: //LNKTYP
							break;
						case 41: //LNKKEYa
							break;
						case 42: //NODTYP
							break;
						case 43: //PROATR 
							break;
						case 44: //PROVAL
							break;*/
                        case 45: //BOX 
                            modal.elementmode = 100;
                            break;
                        case 46: //BOXTYP
                            modal.boxType = br.ReadInt16();
                            break;
                        /*case 47: //PLEX  
							break;*/
                        case 48: //BGNEXTN
                            modal.beginExt = br.ReadInt32();
                            break;
                        case 49: //ENDEXTN
                            modal.endExt = br.ReadInt32();
                            break;
                        /*case 50: //TAPNUM
							break;
						case 51: //TAPCOD
							break;
						case 52: //STRCLS
							break;
						case 53: //RESRVD
							break;
						case 54: //FORMAT
							break;
						case 55: //MASK
							break;
						case 56: //ENDMSK
							break;
						case 57: //LDIRSZ
							break;
						case 58: //SRFNAM
							break;
						case 59: //LIBSCR
							break;*/
                        default:
                            for (int i = 0; i < items; i++)
                            {
                                switch (type)
                                {
                                    case 0:
                                        // No Data
                                        break;
                                    case 1:
                                        // Bit Array
                                        int16 = br.ReadInt16();
                                        break;
                                    case 2:
                                        // Two Byte signed int
                                        int16 = br.ReadInt16();
                                        break;
                                    case 3:
                                        // Four Byte signed int
                                        int32 = br.ReadInt32();
                                        break;
                                    case 4:
                                        // Four byte real
                                        int32 = br.ReadInt32();
                                        break;
                                    case 5:
                                        // Eight byte real
                                        int32 = br.ReadInt32();
                                        break;
                                    case 6:
                                        // ASC
                                        help = br.ReadByte();
                                        break;
                                }
                            }
                            break;
                    }
                    pos = br.BaseStream.Position;
                }

                try
                {
                    drawing_.active_cell = drawing.findCellIndex(cell_.cellName);
                    drawing_.resize((drawing_.userunits / 1E-3));
                }
                catch (Exception)
                {
                    string err = "Unable to find any cells. Is this file legal GDS?";
                    error_msgs.Add(err);
                    throw new Exception(err);
                }

                statusUpdateUI?.Invoke("Done");
                progressUpdateUI?.Invoke(1.0f);
            }
            catch (Exception e)
            {
                valid = false;
                error_msgs.Add(e.Message);
            }
            return valid;
        }

        double read8ByteReal()
        {
            int exp; // help 
            int sig; // help
            int i;
            byte help;   //help
            double double_;
            help = br.ReadByte();
            exp = help & 0x7f;
            if ((help & 0x80) != 0)
            {
                sig = -1;
            }
            else
            {
                sig = 1;
            }
            double_ = 0;
            for (i = 0; i < 7; i++)
            {
                help = br.ReadByte();
                double_ = double_ * 256 + help;
            }
            double_ = sig * double_ * Utils.myPow(16.0, (exp - 64)) / 256 / 256 / 256 / 256 / 256 / 256 / 256;
            return double_;
        }
    }
}