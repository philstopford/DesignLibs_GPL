using gds;
using geoLib;
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
            public List<Element> elements { get; set; }

            public class Element
            {
                public List<bool> isText { get; set; }

                public List<bool> isCellRefArray { get; set; }

                public List<GeoLibArray> arrayData { get; set; }

                public List<string> name { get; set; }

                public List<List<GeoLibPointF>> geometry { get; set; }

                public Element()
                {
                    init();
                }

                void init()
                {
                    geometry = new List<List<GeoLibPointF>>();
                    isText = new List<bool>();
                    name = new List<string>();
                    isCellRefArray = new List<bool>();
                    arrayData = new List<GeoLibArray>();
                }

                public Element(List<List<GeoLibPointF>> sourceGeo)
                {
                    init(sourceGeo);
                }

                void init(List<List<GeoLibPointF>> sourceGeo)
                {
                    geometry = sourceGeo.ToList();
                }

                public void addPoly(List<GeoLibPointF> poly)
                {
                    pAddPoly(poly);
                }

                void pAddPoly(List<GeoLibPointF> poly)
                {
                    geometry.Add(poly.ToList());
                    isText.Add(false);
                    name.Add("");
                }

                public void addText(string name, List<GeoLibPointF> poly)
                {
                    pAddText(name, poly);
                }

                void pAddText(string text, List<GeoLibPointF> poly)
                {
                    geometry.Add(poly.ToList());
                    isText.Add(true);
                    name.Add(text);
                }
            }

            public Structure()
            {
                init();
            }

            void init()
            {
                elements = new List<Element>();
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
            structures[0].addElement();
            structures[0].elements[0].addPoly(new List<GeoLibPointF>() { new GeoLibPointF(0, 0) });

            

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

                        if (ldIndex == -1)
                        {
                            structure_LayerDataTypeList[cellIndex].Add(searchString);
                            structures[cellIndex].addElement();
                            ldIndex = structure_LayerDataTypeList[cellIndex].Count - 1;
                        }

                        getGeometry(ref drawing_, cell, element, hashList, cellIndex, ldIndex);
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

        void getGeometry(ref GCDrawingfield drawing_, int cell, int element, List<string> hashList, int cellIndex, int ldIndex)
        {
            getGeometry(ref drawing_, drawing_.cellList[cell], element, hashList, cellIndex, ldIndex);
        }

        void getGeometry_simple(GCCell gcCell, int element, List<string> hashList, int cellIndex, int ldIndex)
        {
            GCPolygon p = gcCell.elementList[element].convertToPolygon();
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
                    structures[cellIndex].elements[ldIndex].addText(text, t);
                }
                else
                {
                    structures[cellIndex].elements[ldIndex].addPoly(t);
                }
                structures[cellIndex].elements[ldIndex].isCellRefArray.Add(false);
                structures[cellIndex].elements[ldIndex].arrayData.Add(null);
            }
        }

        void getGeometry_complex(GCCell gcCell, int element, List<string> hashList, int cellIndex, int ldIndex)
        {
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
                    if (!tmpCel.elementList[referenceElement].isCellref() && !tmpCel.elementList[referenceElement].isCellrefArray())
                    {
                        getGeometry_2(gcCell, element, hashList, tmpCel, cellIndex, referenceElement, point, xCount, yCount, xSpace, ySpace, angle, mag);
                    }
                    else
                    {
                        getGeometry_complex(tmpCel, element, hashList, cellIndex, ldIndex);
                    }    
                }
            }
        }


        void getGeometry_2(GCCell gcCell, int element, List<string> hashList, GCCell tmpCel, int cellIndex, int referenceElement, GeoLibPoint point, int xCount, int yCount, double xSpace, double ySpace, double angle, double mag)
        {
            Int32 crLayer = tmpCel.elementList[referenceElement].layer_nr;
            Int32 crDatatype = tmpCel.elementList[referenceElement].datatype_nr;

            // See if our layer/datatype combination is known to us already.

            string crSearchString = "L" + crLayer.ToString() + "D" + crDatatype.ToString();

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

                structures[cellIndex].elements[adIndex].isCellRefArray.Add(gcCell.elementList[element].isCellrefArray());

                GeoLibArray tmpArray = new GeoLibArray();
                tmpArray.count = new GeoLibPoint(xCount, yCount);
                tmpArray.point = new GeoLibPoint(point);
                tmpArray.pitch = new GeoLibPoint(xSpace, ySpace);
                structures[cellIndex].elements[adIndex].arrayData.Add(tmpArray);

                crLDIndex = structure_LayerDataTypeList[cellIndex].Count - 1;
            }

            try
            {
                GCPolygon crP = tmpCel.elementList[referenceElement].convertToPolygon();
                crP.move(point);
                crP.rotate(angle, point);
                crP.scale(point, mag);

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
                    structures[cellIndex].elements[crLDIndex].addPoly(t);
                }
            }
            catch (Exception)
            {
                // Exception is not a big deal.
            }
        }

        void getGeometry(ref GCDrawingfield drawing_, GCCell gcCell, int element, List<string> hashList, int cellIndex, int ldIndex)
        {
            // Now we have to process our geometry into the right place.
            if (!gcCell.elementList[element].isCellref() && !gcCell.elementList[element].isCellrefArray())
            {
                getGeometry_simple(gcCell, element, hashList, cellIndex, ldIndex);
            }
            else
            {
                getGeometry_complex(gcCell, element, hashList, cellIndex, ldIndex);
            }
        }


        public List<GeoLibArray> getArrayParameters()
        {
            return pGetArrayParameters();
        }

        List<GeoLibArray> pGetArrayParameters()
        {
            return structures[activeStructure].elements[activeLD].arrayData;
        }

        public List<GeoLibPointF[]> points()
        {
            return pPoints();
        }

        List<GeoLibPointF[]> pPoints()
        {
            List<GeoLibPointF[]> points = new List<GeoLibPointF[]>();
            for (Int32 poly = 0; poly < structures[activeStructure].elements[activeLD].geometry.Count; poly++)
            {
                points.Add(structures[activeStructure].elements[activeLD].geometry[poly].ToArray());
            }
            return points;
        }

        public List<bool> text()
        {
            return pText();
        }

        public List<bool> pText()
        {
            return structures[activeStructure].elements[activeLD].isText.ToList();
        }

        public string getName(int index)
        {
            return pGetName(index);
        }

        string pGetName(int index)
        {
            return structures[activeStructure].elements[activeLD].name[index];
        }

        public void updateGeometry(Int32 structure, Int32 layerdatatype)
        {
            pUpdateGeometry(structure, layerdatatype);
        }

        void pUpdateGeometry(Int32 structure, Int32 layerdatatype)
        {
            if (valid)
            {
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
