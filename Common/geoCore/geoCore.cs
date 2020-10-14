using gds;
using geoLib;
using geoWrangler;
using oasis;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;

namespace geoCoreLib
{
    public partial class GeoCore
    {
        public double scaling = 1.0f;
        GCDrawingfield drawingField;

        public static UInt16 maxLayers = UInt16.MaxValue;

        public double baseScale = 1.0;

        Boolean valid;

        public void setValid(bool val)
        {
            pSetValid(val);
        }

        void pSetValid(bool val)
        {
            valid = val;
        }

        public bool isValid()
        {
            return pIsValid();
        }

        bool pIsValid()
        {
            return valid;
        }

        string filename;
        int fileFormat;
        public enum fileType { gds, oasis }

        Dictionary<string, string> layerNames;

        List<string> pStructureList;

        public List<string> getStructureList()
        {
            return pGetStructureList();
        }

        List<string> pGetStructureList()
        {
            return pStructureList;
        }

        List<string> pActiveStructure_LDList;

        public List<string> getActiveStructureLDList()
        {
            return pGetActiveStructureLDList();
        }

        List<string> pGetActiveStructureLDList()
        {
            return pActiveStructure_LDList;
        }

        public ObservableCollection<string> structureList_ { get; set; }
        public ObservableCollection<string> activeStructure_LayerDataTypeList_ { get; set; }

        public int activeStructure { get; set; }
        public int activeLD { get; set; }
        List<List<string>> structure_LayerDataTypeList;
        public List<List<string>> getStructureLayerDataTypeList()
        {
            return pGetStructureLayerDataTypeList();
        }

        List<List<string>> pGetStructureLayerDataTypeList()
        {
            return structure_LayerDataTypeList;
        }

        class Structure
        {
            public List<string> layerDataTypes;

            public List<Element> elements;

            List<string> bakedGeo_LD;
            List<BakedGeo> bakedGeo;

            public class BakedGeo
            {
                public string LD { get; set; }
                public List<GeoLibPointF[]> fgeo { get; set; }
                public List<bool> isText { get; set; }

                public BakedGeo()
                {
                    fgeo = new List<GeoLibPointF[]>();
                    isText = new List<bool>();
                }

                public BakedGeo(List<GeoLibPointF[]> source, List<bool> text, string ld)
                {
                    fgeo = source;
                    isText = text;
                    LD = ld;
                }
            }

            public void addBakedGeo(List<GeoLibPointF[]> source, List<bool> text, string ld)
            {
                for (int i = 0; i < source.Count; i++)
                {
                    addBakedGeo(source[i], text[i], ld);
                }
            }

            public void addBakedGeo(GeoLibPointF[] source, bool text, string ld)
            {
                int index = bakedGeo_LD.IndexOf(ld);
                if (index == -1)
                {
                    bakedGeo.Add(new BakedGeo(new List<GeoLibPointF[]>() { source }, new List<bool>() { text }, ld));
                    bakedGeo_LD.Add(ld);
                }
                else
                {
                    bakedGeo[index].fgeo.Add(source);
                    bakedGeo[index].isText.Add(text);
                }
            }

            public List<GeoLibPointF[]> getBakedGeo(string ld)
            {
                int index = bakedGeo_LD.IndexOf(ld);
                if (index == -1)
                {
                    return new List<GeoLibPointF[]>();
                }
                else
                {
                    return bakedGeo[index].fgeo;
                }
            }

            public List<GeoLibArray> getArrayData(string ld)
            {
                List<GeoLibArray> ret = new List<GeoLibArray>();
                for (int i = 0; i < elements.Count; i++)
                {
                    //if (elements[i].LD == ld)
                    {
                        ret.Add(elements[i].arrayData);
                    }
                }

                return ret;
            }

            public List<bool> isText(string ld)
            {
                List<bool> ret = new List<bool>();
                for (int i = 0; i < bakedGeo_LD.Count; i++)
                {
                    if (bakedGeo_LD[i] == ld)
                    {
                        ret.AddRange(bakedGeo[i].isText);
                    }
                }

                return ret;
            }

            public class Element
            {
                public string LD; // might be useful to track originating LD.

                public bool isText { get; set; }
                public string isCellRefArray { get; set; }

                public List<GeoLibPointF> geometry { get; set; }

                public GeoLibArray arrayData { get; set; }

                public string name { get; set; }

                public Element()
                {
                    init();
                }

                void init()
                {
                    geometry = new List<GeoLibPointF>();
                    isText = false;
                    name = "";
                    isCellRefArray = "";
                    arrayData = null;
                }

                public Element(List<GeoLibPointF> sourceGeo, bool text)
                {
                    init(sourceGeo, text);
                }

                void init(List<GeoLibPointF> sourceGeo, bool text)
                {
                    geometry = sourceGeo.ToList();
                    isText = text;
                }
            }

            public void addPoly(List<GeoLibPointF> poly, string ldString)
            {
                pAddPoly(poly, ldString);
            }

            void pAddPoly(List<GeoLibPointF> poly, string ldString)
            {
                elements.Add(new Element(poly, false));
                elements[elements.Count - 1].LD = ldString;
                if (ldString == "L-1D-1")
                {
                    int xx = 2;
                }
                if (layerDataTypes.IndexOf(ldString) == -1)
                {
                    layerDataTypes.Add(ldString);
                }
            }

            public void addText(string name, List<GeoLibPointF> poly, string ldString)
            {
                pAddText(name, poly, ldString);
            }

            void pAddText(string text, List<GeoLibPointF> poly, string ldString)
            {
                elements.Add(new Element(poly, true));
                elements[elements.Count - 1].name = text;
                elements[elements.Count - 1].LD = ldString;
                if (ldString == "L-1D-1")
                {
                    int xx = 2;
                }
                if (layerDataTypes.IndexOf(ldString) == -1)
                {
                    layerDataTypes.Add(ldString);
                }
            }

            public Structure()
            {
                init();
            }

            void init()
            {
                elements = new List<Element>();
                layerDataTypes = new List<string>();
                bakedGeo = new List<BakedGeo>();
                bakedGeo_LD = new List<string>();
            }

            public void addElement()
            {
                pAddElement();
            }

            void pAddElement()
            {
                elements.Add(new Element());
            }
        }

        List<Structure> structures;

        public void updateCollections()
        {
            pUpdateCollections();
        }

        void pUpdateCollections()
        {
            structureList_.Clear();
            for (int i = 0; i < pStructureList.Count; i++)
            {
                structureList_.Add(pStructureList[i]);
            }
            activeStructure_LayerDataTypeList_.Clear();
            for (int i = 0; i < pActiveStructure_LDList.Count; i++)
            {
                activeStructure_LayerDataTypeList_.Add(pActiveStructure_LDList[i]);
            }
        }

        public void updateGeoCore(string filename, fileType type)
        {
            pUpdateGeoCore(filename, type);
        }

        void pUpdateGeoCore(string filename, fileType type)
        {
            reset();
            this.filename = filename;

            if (filename == "")
            {
                return;
            }

            if (!File.Exists(filename))
            {
                throw new Exception("File not accessible: " + filename);
            }

            /*
			structureList.Clear();
			structure_LayerDataTypeList[0].Clear();
			structure_Element_Poly_PointList.Clear();
			*/
            switch (type)
            {
                case fileType.gds:
                    gdsReader gdsFileData = new gdsReader(filename);
                    valid = gdsFileData.load(ref drawingField);
                    fileFormat = (int)fileType.gds;
                    if (valid)
                    {
                        processGeometry(ref drawingField, gdsFileData.layerNames);
                        layerNames = gdsFileData.layerNames;
                        activeStructure = drawingField.active_cell;
                    }
                    break;
                case fileType.oasis:
                    oasReader oasFileData = new oasReader(filename);
                    valid = oasFileData.load(ref drawingField);
                    fileFormat = (int)fileType.oasis;
                    if (valid)
                    {
                        processGeometry(ref drawingField, oasFileData.layerNames);
                        layerNames = oasFileData.layerNames;
                        activeStructure = drawingField.active_cell;
                    }
                    break;
                default:
                    valid = false;
                    break;
            }
            genLDList();
        }

        public GCDrawingfield getDrawing()
        {
            return pGetDrawing();
        }

        GCDrawingfield pGetDrawing()
        {
            if (valid)
            {
                return drawingField;
            }

            return null;
        }

        public void setDrawing(GCDrawingfield source)
        {
            drawingField = source.copy();
        }

        public Dictionary<string, string> getLayerNames()
        {
            return pGetLayerNames();
        }

        Dictionary<string, string> pGetLayerNames()
        {
            return layerNames;
        }

        public void addLayerName(string key, string value)
        {
            pAddLayerName(key, value);
        }

        void pAddLayerName(string key, string value)
        {
            layerNames.Add(key, value);
        }

        public GeoCore()
        {
            pGeoCore();
        }

        void pGeoCore()
        {
            drawingField = new GCDrawingfield("");
            pStructureList = new List<string>();
            pStructureList.Add("");
            pActiveStructure_LDList = new List<string>();
            pActiveStructure_LDList.Add("");
            structureList_ = new ObservableCollection<string>();
            activeStructure_LayerDataTypeList_ = new ObservableCollection<string>();
            updateCollections();
            reset();
            genLDList();
        }

        public void reset()
        {
            pReset();
        }

        void pReset()
        {
            drawingField.reset();

            valid = false;
            activeStructure = 0;
            filename = "";

            layerNames = new Dictionary<string, string>();

            pStructureList.Clear();
            pStructureList.Add("");
            pActiveStructure_LDList.Clear();
            pActiveStructure_LDList.Add("");

            structure_LayerDataTypeList = new List<List<string>>();
            structure_LayerDataTypeList.Add(new List<string>());
            structure_LayerDataTypeList[0].Add("");

            structures = new List<Structure>();
            structures.Add(new Structure());
            //structures[0].addElement();

            pUpdateCollections();
        }

        void genLDList()
        {
            pActiveStructure_LDList.Clear();
            for (int ld = 0; ld < structure_LayerDataTypeList[activeStructure].Count; ld++)
            {
                pActiveStructure_LDList.Add(structure_LayerDataTypeList[activeStructure][ld].ToString());
            }
        }

        public void readValues(GeoCore sourceGeoCore)
        {
            pReadValues(sourceGeoCore);
        }

        void pReadValues(GeoCore sourceGeoCore)
        {
            reset();
            filename = sourceGeoCore.filename;
            int source = 0;

            pStructureList[source] = sourceGeoCore.pStructureList[source];
            structure_LayerDataTypeList = sourceGeoCore.structure_LayerDataTypeList.ToList();
            source++;

            while (source < sourceGeoCore.pStructureList.Count())
            {
                pStructureList.Add(sourceGeoCore.pStructureList[source]);
                source++;
            }

            structures = sourceGeoCore.structures.ToList();
            valid = sourceGeoCore.valid;
            activeStructure = sourceGeoCore.activeStructure;
            genLDList();
        }

        void processGeometry(ref GCDrawingfield drawing_, Dictionary<string, string> layerNames)
        {
            // Should have been reset before this call. Remove the defaults.
            pStructureList.Clear();
            structure_LayerDataTypeList[0].Clear();
            structures.Clear();

            scaling = 1.0f;

            // Set up the scaling for the conversion
            if (valid && (fileFormat == (int)fileType.gds))
            {
                scaling = drawing_.userunits;
                drawing_.databaseunits = 1 / drawing_.userunits;
            }
            else if (valid && (fileFormat == (int)fileType.oasis))
            {
                scaling = 1.0 / drawing_.databaseunits;
                drawing_.userunits = 1 / drawing_.databaseunits;
            }

            scaling *= baseScale; // this is used to compensate for any reader application needing scaled geometry.

            for (int cell = 0; cell < drawing_.cellList.Count; cell++)
            {
                if (drawing_.cellList[cell].elementList != null)
                {
                    // Check whether our cell is known already.
                    string cellName = drawing_.cellList[cell].cellName;

                    Int32 cellIndex = -1;

                    try
                    {
                        cellIndex = pStructureList.IndexOf(cellName);
                    }
                    catch (Exception)
                    {
                    }

                    if (cellIndex == -1)
                    {
                        pStructureList.Add(cellName); // new cell.
                        cellIndex = pStructureList.Count - 1; // index.

                        // Cell marker is null so we need to skip in this case.
                        structure_LayerDataTypeList.Add(new List<string>());

                        structures.Add(new Structure());
                    }

                    // Check whether our layer / datatype combination is known.
                    // We have to be careful here - Oasis has the ability to have strings, GDS is numeric.

                    List<string> hashList = new List<string>();
                    for (int element = 0; element < drawing_.cellList[cell].elementList.Count; element++)
                    {
                        Int32 layer = drawing_.cellList[cell].elementList[element].layer_nr;
                        Int32 datatype = drawing_.cellList[cell].elementList[element].datatype_nr;

                        // See if our layer/datatype combination is known to us already.

                        string searchString = "L" + layer.ToString() + "D" + datatype.ToString();

                        if (searchString == "L-1D-1")
                        {
                            // We get into trouble here as the cell reference needs to be followed.... Ugh.
                            int xx = 2;
                        }

                        // Query our dictionary.
                        string resultString = "";
                        try
                        {
                            if (layerNames.TryGetValue(searchString, out resultString))
                            {
                                searchString = resultString;
                            }
                        }
                        catch (Exception)
                        {

                        }

                        Int32 ldIndex = -1;

                        try
                        {
                            ldIndex = structure_LayerDataTypeList[cellIndex].IndexOf(searchString);
                        }
                        catch (Exception)
                        {
                        }

                        if ((ldIndex == -1) && (searchString != "L-1D-1"))
                        {
                            structure_LayerDataTypeList[cellIndex].Add(searchString);
                            //structures[cellIndex].addElement();
                            ldIndex = structure_LayerDataTypeList[cellIndex].Count - 1;
                        }

                        getGeometry(ref drawing_, cell, element, hashList, cellIndex, searchString);
                    }
                }
            }

            // Final clean-up
            structure_LayerDataTypeList.RemoveAt(structure_LayerDataTypeList.Count - 1);
            // Need to remove any layer/datatypes that don't have any entries.

            for (int structure = 0; structure < pStructureList.Count; structure++)
            {
                List<Int32> indicesToRemove = new List<Int32>();
                for (int ld = structure_LayerDataTypeList[structure].Count - 1; ld > -1; ld--)
                {
                    if (structure_LayerDataTypeList[structure][ld] == "L-1D-1")
                    {
                        indicesToRemove.Add(ld);
                    }
                }

                for (int index = 0; index < indicesToRemove.Count(); index++)
                {
                    structure_LayerDataTypeList[structure].RemoveAt(index);
                    structures[structure].elements.RemoveAt(index);
                }
            }
        }

        void getGeometry(ref GCDrawingfield drawing_, int cell, int element, List<string> hashList, int cellIndex, string searchString)
        {
            getGeometry(ref drawing_, drawing_.cellList[cell], element, hashList, cellIndex, searchString);
        }

        GeoData getGeometry_simple(GCCell gcCell, int element, List<string> hashList, int cellIndex)
        {
            List<List<GeoLibPointF>> ret = new List<List<GeoLibPointF>>();
            List<string> lds = new List<string>();
            List<GCPolygon> lp = gcCell.elementList[element].convertToPolygons();
            List<bool> isText = new List<bool>();
            for (int poly = 0; poly < lp.Count; poly++)
            {
                GCPolygon p = lp[poly];
                isText.Add(p.isText());

                string ldString = "L" + p.layer_nr.ToString() + "D" + p.datatype_nr.ToString();
                if (ldString == "L-1D-1")
                {
                    int xx = 2;
                }
                // We should remove identical polygons here in case of doubled-up input geometry.
                string crP_Hash = utility.Utils.GetMD5Hash(p.pointarray);

                if (hashList.IndexOf(crP_Hash) == -1)
                {
                    hashList.Add(crP_Hash);
                    List<GeoLibPointF> t = new List<GeoLibPointF>();
                    for (int pt = 0; pt < p.pointarray.Length; pt++)
                    {
                        t.Add(new GeoLibPointF(p.pointarray[pt].X * scaling, p.pointarray[pt].Y * scaling));
                    }
                    if (gcCell.elementList[element].isText())
                    {
                        string text = gcCell.elementList[element].getName();
                        structures[cellIndex].addText(text, t, ldString);
                    }
                    else
                    {
                        structures[cellIndex].addPoly(t, ldString);
                    }
                    ret.Add(t);
                    lds.Add(ldString);
                }
            }

            return new GeoData(ret, isText, lds);
        }

        GeoData getGeometry_complex(GCCell gcCell, int element, List<string> hashList, int cellIndex)
        {
            GeoData ret = new GeoData();
            // Need to de-reference these cases.
            double angle;
            GeoLibPoint point;
            double mag;
            GCCell tmpCel;
            double xSpace = 0;
            double ySpace = 0;
            int xCount = 1;
            int yCount = 1;

            if (gcCell.elementList[element].isCellref())
            {
                GCCellref refCell = (GCCellref)gcCell.elementList[element];
                point = refCell.getPos();
                mag = refCell.trans.mag;
                angle = refCell.trans.angle;
                tmpCel = refCell.cell_ref;
            }
            else
            {
                GCCellRefArray refCell = (GCCellRefArray)gcCell.elementList[element];
                point = refCell.getPos();
                mag = refCell.trans.mag;
                angle = refCell.trans.angle;
                tmpCel = refCell.cell_ref;
                xSpace = refCell.pitch.X;
                ySpace = refCell.pitch.Y;
                xCount = refCell.count_x;
                yCount = refCell.count_y;
            }
            if (tmpCel != null) // guard against broken cellref
            {
                for (int referenceElement = 0; referenceElement < tmpCel.elementList.Count; referenceElement++)
                {
                    if (!tmpCel.elementList[referenceElement].isCellref() && (!tmpCel.elementList[referenceElement].GetType().Equals(typeof(GCCellRefArray))))
                    {
                        ret = getGeometry_2(gcCell, element, hashList, tmpCel, cellIndex, referenceElement, point, xCount, yCount, xSpace, ySpace, angle, mag);
                    }
                    else
                    {
                        ret = getGeometry_complex(tmpCel, element, hashList, cellIndex);
                    }    
                }
            }

            // Need to construct the array of our geometry here!
            GeoData postArray = new GeoData();
            for (int i = 0; i < ret.geo.Count; i++)
            {
                postArray.geo = GeoWrangler.makeArray2(ret.geo[i], xCount, xSpace, yCount, ySpace);
                postArray.ld = new List<string>();
                postArray.isText = new List<bool>();
                for (int j = 0; j < postArray.geo.Count; j++)
                {
                    postArray.ld.Add(ret.ld[i]);
                    postArray.isText.Add(ret.isText[i]);
                }
            }
            // return GeoWrangler.makeArray(ret, xCount, xSpace, yCount, ySpace);

            return postArray;

            // Might need to track nested  array configurations here to handle recursive settings.
        }


        GeoData getGeometry_2(GCCell gcCell, int element, List<string> hashList, GCCell tmpCel, int cellIndex, int referenceElement, GeoLibPoint point, int xCount, int yCount, double xSpace, double ySpace, double angle, double mag)
        {
            List<List<GeoLibPointF>> ret = new List<List<GeoLibPointF>>();
            Int32 crLayer = tmpCel.elementList[referenceElement].layer_nr;
            Int32 crDatatype = tmpCel.elementList[referenceElement].datatype_nr;

            List<string> lds = new List<string>();

            List<bool> text = new List<bool>();

            // See if our layer/datatype combination is known to us already.

            string crSearchString = "L" + crLayer.ToString() + "D" + crDatatype.ToString();

            if (crSearchString == "L-1D-1")
            {
                int xx = 2;
            }

            Int32 crLDIndex = -1;

            try
            {
                crLDIndex = structure_LayerDataTypeList[cellIndex].IndexOf(crSearchString);
            }
            catch (Exception)
            {
            }

            if (crLDIndex == -1)
            {
                structure_LayerDataTypeList[cellIndex].Add(crSearchString);
                structures[cellIndex].addElement();
                int adIndex = structures[cellIndex].elements.Count - 1;

                if (gcCell.elementList[element].GetType().Equals(typeof(GCCellRefArray)))
                {
                    structures[cellIndex].elements[adIndex].isCellRefArray = gcCell.elementList[element].depend().cellName;// gcCell.elementList[element].isCellrefArray());
                }

                GeoLibArray tmpArray = new GeoLibArray();
                tmpArray.count = new GeoLibPoint(xCount, yCount);
                tmpArray.point = new GeoLibPoint(point);
                tmpArray.pitch = new GeoLibPoint(xSpace, ySpace);
                structures[cellIndex].elements[adIndex].arrayData = tmpArray;
            }

            try
            {
                List<GCPolygon> lCRP = tmpCel.elementList[referenceElement].convertToPolygons();
                for (int p = 0; p < lCRP.Count; p++)
                {
                    GCPolygon crP = lCRP[p];
                    crP.move(point);
                    crP.rotate(angle, point);
                    crP.scale(point, mag);

                    string ldString = "L" + crP.layer_nr.ToString() + "D" + crP.datatype_nr.ToString();

                    // We should remove identical polygons here in case of doubled-up input geometry.
                    string crP_Hash = utility.Utils.GetMD5Hash(crP.pointarray);

                    if (hashList.IndexOf(crP_Hash) == -1)
                    {
                        hashList.Add(crP_Hash);

                        int x = 0;
                        int y = 0;
                        List<GeoLibPointF> t = new List<GeoLibPointF>();
                        for (int pt = 0; pt < crP.pointarray.Length; pt++)
                        {
                            t.Add(new GeoLibPointF((crP.pointarray[pt].X + (x * xSpace)) * scaling, (crP.pointarray[pt].Y + (y * ySpace)) * scaling));
                        }
                        structures[cellIndex].addPoly(t, ldString);
                        ret.Add(t);
                        text.Add(crP.isText());
                        lds.Add(ldString);
                    }
                }
            }
            catch (Exception)
            {
                // Exception is not a big deal.
            }

            return new GeoData(ret, text, lds);
        }


        public class GeoData
        {
            public List<List<GeoLibPointF>> geo;
            public List<string> ld;
            public List<bool> isText;

            public GeoData()
            {
                geo = new List<List<GeoLibPointF>>();
                ld = new List<string>();
                isText = new List<bool>();
            }

            public GeoData(List<List<GeoLibPointF>> poly, List<bool> text, List<string> LDs)
            {
                geo = poly.ToList();
                ld = LDs.ToList();
                isText = text.ToList();
            }
        }

        void getGeometry(ref GCDrawingfield drawing_, GCCell gcCell, int element, List<string> hashList, int cellIndex, string searchString)
        {
            GeoData ret;
            // Now we have to process our geometry into the right place.
            if (!gcCell.elementList[element].isCellref() && !gcCell.elementList[element].isCellrefArray())
            {
                ret = getGeometry_simple(gcCell, element, hashList, cellIndex);
            }
            else
            {
                ret = getGeometry_complex(gcCell, element, hashList, cellIndex);
            }

            if (searchString == "L-1D-1")
            {
                int x = 2;
            }

            for (int i = 0; i < ret.geo.Count; i++)
            {
                structures[cellIndex].addBakedGeo(ret.geo[i].ToArray(), ret.isText[i], ret.ld[i]);
            }
        }

        public List<GeoLibArray> getArrayParameters()
        {
            return pGetArrayParameters();
        }

        List<GeoLibArray> pGetArrayParameters()
        {
            return structures[activeStructure].getArrayData(pActiveStructure_LDList[activeLD]);
        }

        public List<bool> isText()
        {
            return pIsText();
        }

        List<bool> pIsText()
        {
            string text = pActiveStructure_LDList[activeLD];
            return structures[activeStructure].isText(text);
        }

        public List<GeoLibPointF[]> points(bool flatten)
        {
            return pPoints(flatten);
        }

        List<GeoLibPointF[]> pPoints(bool flatten)
        {
            if (flatten)
            {
                return structures[activeStructure].getBakedGeo(pActiveStructure_LDList[activeLD]);
            }

            // Do we ever get here?

            List<GeoLibPointF[]> points = new List<GeoLibPointF[]>();
            List<GeoLibPoint> array_count = new List<GeoLibPoint>();
            List<GeoLibPointF> array_pitch = new List<GeoLibPointF>();

            points.Add(structures[activeStructure].elements[activeLD].geometry.ToArray());
            if (structures[activeStructure].elements[activeLD].arrayData != null)
            {
                array_count.Add(new GeoLibPoint(structures[activeStructure].elements[activeLD].arrayData.count));
                array_pitch.Add(new GeoLibPointF(structures[activeStructure].elements[activeLD].arrayData.pitch));
            }
            else
            {
                array_count.Add(new GeoLibPoint(1, 1));
                array_pitch.Add(new GeoLibPointF(0, 0));
            }

            /*
            // Any arrays?
            if (flatten)
            {
                List<GeoLibPointF[]> arrayed = new List<GeoLibPointF[]>();
                for (int i = 0; i < points.Count; i++)
                {
                    List<GeoLibPointF[]> fa = GeoWrangler.makeArray(points[i], array_count[i].X, array_pitch[i].X, array_count[i].Y, array_pitch[i].Y);
                    arrayed.AddRange(fa);
                }

                points = arrayed;
            }
            */
            // Update the geo.
            // structures[activeLD].fgeo[activeLD] = points.ToList();

            return points;
        }

        public string getName(int index)
        {
            return pGetName(index);
        }

        string pGetName(int index)
        {
            return structures[activeStructure].elements[index].name;
        }

        public void updateGeometry(Int32 structure, Int32 layerdatatype)
        {
            pUpdateGeometry(structure, layerdatatype);
        }

        void pUpdateGeometry(Int32 structure, Int32 layerdatatype)
        {
            if (valid)
            {
                drawingField.active_cell = structure;
                activeStructure = structure;
                activeLD = layerdatatype;
                genLDList();
                if (structures[structure].elements[layerdatatype] == null)
                {
                    return;
                }
            }
        }
    }
}
