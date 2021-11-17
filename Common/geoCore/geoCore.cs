using gds;
using geoLib;
using geoWrangler;
using oasis;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;

namespace geoCoreLib;

public class GeoCore
{
    public double scaling = 1.0f;
    private GCDrawingfield drawingField;

    public static ushort maxLayers = ushort.MaxValue;

    public double baseScale = 1.0;

    public List<string> error_msgs;
    private bool valid;

    public void setValid(bool val)
    {
        pSetValid(val);
    }

    private void pSetValid(bool val)
    {
        valid = val;
    }

    public bool isValid()
    {
        return pIsValid();
    }

    private bool pIsValid()
    {
        return valid;
    }

    private string filename;
    private int fileFormat;
    public enum fileType { gds, oasis }

    private Dictionary<string, string> layerNames;

    private List<string> pStructureList;

    public List<string> getStructureList()
    {
        return pGetStructureList();
    }

    private List<string> pGetStructureList()
    {
        return pStructureList;
    }

    private List<string> pActiveStructure_LDList;

    public List<string> getActiveStructureLDList()
    {
        return pGetActiveStructureLDList();
    }

    private List<string> pGetActiveStructureLDList()
    {
        return pActiveStructure_LDList;
    }

    public ObservableCollection<string> structureList_ { get; set; }
    public ObservableCollection<string> activeStructure_LayerDataTypeList_ { get; set; }

    public int activeStructure { get; set; }
    public int activeLD { get; set; }
    private List<List<string>> structure_LayerDataTypeList;
    public List<List<string>> getStructureLayerDataTypeList()
    {
        return pGetStructureLayerDataTypeList();
    }

    private List<List<string>> pGetStructureLayerDataTypeList()
    {
        return structure_LayerDataTypeList;
    }

    private class Structure
    {
        public List<string> layerDataTypes;

        public List<Element> elements;

        private List<string> bakedGeo_LD;
        private List<BakedGeo> bakedGeo;

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
            switch (index)
            {
                case -1:
                    bakedGeo.Add(new BakedGeo(new List<GeoLibPointF[]> { source }, new List<bool> { text }, ld));
                    bakedGeo_LD.Add(ld);
                    break;
                default:
                    bakedGeo[index].fgeo.Add(source);
                    bakedGeo[index].isText.Add(text);
                    break;
            }
        }

        public List<GeoLibPointF[]> getBakedGeo(string ld)
        {
            int index = bakedGeo_LD.IndexOf(ld);
            return index switch
            {
                -1 => new List<GeoLibPointF[]>(),
                _ => bakedGeo[index].fgeo
            };
        }

        public List<GeoLibArray> getArrayData(string ld)
        {
            List<GeoLibArray> ret = new();
            foreach (Element t in elements)
            {
                //if (elements[i].LD == ld)
                {
                    ret.Add(t.arrayData);
                }
            }

            return ret;
        }

        public List<bool> isText(string ld)
        {
            List<bool> ret = new();
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

            private void init()
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

            private void init(List<GeoLibPointF> sourceGeo, bool text)
            {
                geometry = sourceGeo.ToList();
                isText = text;
            }
        }

        public void addPoly(List<GeoLibPointF> poly, string ldString)
        {
            pAddPoly(poly, ldString);
        }

        private void pAddPoly(List<GeoLibPointF> poly, string ldString)
        {
            elements.Add(new Element(poly, false));
            elements[^1].LD = ldString;
            switch (layerDataTypes.IndexOf(ldString))
            {
                case -1:
                    layerDataTypes.Add(ldString);
                    break;
            }
        }

        public void addText(string name, List<GeoLibPointF> poly, string ldString)
        {
            pAddText(name, poly, ldString);
        }

        private void pAddText(string text, List<GeoLibPointF> poly, string ldString)
        {
            elements.Add(new Element(poly, true));
            elements[^1].name = text;
            elements[^1].LD = ldString;
            switch (layerDataTypes.IndexOf(ldString))
            {
                case -1:
                    layerDataTypes.Add(ldString);
                    break;
            }
        }

        public Structure()
        {
            init();
        }

        private void init()
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

        private void pAddElement()
        {
            elements.Add(new Element());
        }
    }

    private List<Structure> structures;

    public void updateCollections()
    {
        pUpdateCollections();
    }

    private void pUpdateCollections()
    {
        structureList_.Clear();
        foreach (string t in pStructureList)
        {
            structureList_.Add(t);
        }
        activeStructure_LayerDataTypeList_.Clear();
        foreach (string t in pActiveStructure_LDList)
        {
            activeStructure_LayerDataTypeList_.Add(t);
        }
    }

    public void updateGeoCore(string filename_, fileType type)
    {
        pUpdateGeoCore(filename_, type);
    }

    private void pUpdateGeoCore(string filename_, fileType type)
    {
        reset();
        filename = filename_;

        switch (filename)
        {
            case "":
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
                gdsReader gdsFileData = new(filename);
                valid = gdsFileData.load(ref drawingField);
                fileFormat = (int)fileType.gds;
                error_msgs = gdsFileData.error_msgs;
                switch (valid)
                {
                    case true:
                        processGeometry(ref drawingField, gdsFileData.layerNames);
                        layerNames = gdsFileData.layerNames;
                        activeStructure = drawingField.active_cell;
                        break;
                }
                break;
            case fileType.oasis:
                oasReader oasFileData = new(filename);
                valid = oasFileData.load(ref drawingField);
                fileFormat = (int)fileType.oasis;
                error_msgs = oasFileData.error_msgs;
                switch (valid)
                {
                    case true:
                        processGeometry(ref drawingField, oasFileData.layerNames);
                        layerNames = oasFileData.layerNames;
                        activeStructure = drawingField.active_cell;
                        break;
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

    private GCDrawingfield pGetDrawing()
    {
        return valid switch
        {
            true => drawingField,
            _ => null
        };
    }

    public void setDrawing(GCDrawingfield source)
    {
        drawingField = source.copy();
    }

    public Dictionary<string, string> getLayerNames()
    {
        return pGetLayerNames();
    }

    private Dictionary<string, string> pGetLayerNames()
    {
        return layerNames;
    }

    public void addLayerName(string key, string value)
    {
        pAddLayerName(key, value);
    }

    private void pAddLayerName(string key, string value)
    {
        layerNames.Add(key, value);
    }

    public GeoCore()
    {
        pGeoCore();
    }

    private void pGeoCore()
    {
        drawingField = new GCDrawingfield("");
        pStructureList = new List<string> {""};
        pActiveStructure_LDList = new List<string> {""};
        structureList_ = new ObservableCollection<string>();
        activeStructure_LayerDataTypeList_ = new ObservableCollection<string>();
        error_msgs = new List<string>();
        updateCollections();
        reset();
        genLDList();
    }

    public void reset()
    {
        pReset();
    }

    private void pReset()
    {
        drawingField.reset();

        valid = false;
        error_msgs.Clear();
        activeStructure = 0;
        filename = "";

        layerNames = new Dictionary<string, string>();

        pStructureList.Clear();
        pStructureList.Add("");
        pActiveStructure_LDList.Clear();
        pActiveStructure_LDList.Add("");

        structure_LayerDataTypeList = new List<List<string>> {new()};
        structure_LayerDataTypeList[0].Add("");

        structures = new List<Structure> {new()};
        //structures[0].addElement();

        pUpdateCollections();
    }

    private void genLDList()
    {
        pActiveStructure_LDList.Clear();
        for (int ld = 0; ld < structure_LayerDataTypeList[activeStructure].Count; ld++)
        {
            pActiveStructure_LDList.Add(structure_LayerDataTypeList[activeStructure][ld]);
        }
    }

    public void readValues(GeoCore sourceGeoCore)
    {
        pReadValues(sourceGeoCore);
    }

    private void pReadValues(GeoCore sourceGeoCore)
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

    private void processGeometry(ref GCDrawingfield drawing_, Dictionary<string, string> layerNames_)
    {
        // Should have been reset before this call. Remove the defaults.
        pStructureList.Clear();
        structure_LayerDataTypeList[0].Clear();
        structures.Clear();

        scaling = 1.0f;
            
        for (int cell = 0; cell < drawing_.cellList.Count; cell++)
        {
            if (drawing_.cellList[cell].elementList != null)
            {
                // Check whether our cell is known already.
                string cellName = drawing_.cellList[cell].cellName;

                int cellIndex = -1;

                try
                {
                    cellIndex = pStructureList.IndexOf(cellName);
                }
                catch (Exception)
                {
                }

                switch (cellIndex)
                {
                    case -1:
                        pStructureList.Add(cellName); // new cell.
                        cellIndex = pStructureList.Count - 1; // index.

                        // Cell marker is null so we need to skip in this case.
                        structure_LayerDataTypeList.Add(new List<string>());

                        structures.Add(new Structure());
                        break;
                }

                // Check whether our layer / datatype combination is known.
                // We have to be careful here - Oasis has the ability to have strings, GDS is numeric.

                List<string> hashList = new();
                for (int element = 0; element < drawing_.cellList[cell].elementList.Count; element++)
                {
                    int layer = drawing_.cellList[cell].elementList[element].layer_nr;
                    int datatype = drawing_.cellList[cell].elementList[element].datatype_nr;

                    // See if our layer/datatype combination is known to us already.

                    string searchString = "L" + layer + "D" + datatype;

                    // Query our dictionary.
                    try
                    {
                        string resultString;
                        if (layerNames_.TryGetValue(searchString, out resultString))
                        {
                            searchString = resultString;
                        }
                    }
                    catch (Exception)
                    {

                    }

                    int ldIndex = -1;

                    try
                    {
                        ldIndex = structure_LayerDataTypeList[cellIndex].IndexOf(searchString);
                    }
                    catch (Exception)
                    {
                    }

                    switch (ldIndex)
                    {
                        case -1 when searchString != "L-1D-1":
                            structure_LayerDataTypeList[cellIndex].Add(searchString);
                            //structures[cellIndex].addElement();
                            ldIndex = structure_LayerDataTypeList[cellIndex].Count - 1;
                            break;
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
            List<int> indicesToRemove = new();
            for (int ld = structure_LayerDataTypeList[structure].Count - 1; ld > -1; ld--)
            {
                switch (structure_LayerDataTypeList[structure][ld])
                {
                    case "L-1D-1":
                        indicesToRemove.Add(ld);
                        break;
                }
            }

            for (int index = 0; index < indicesToRemove.Count(); index++)
            {
                structure_LayerDataTypeList[structure].RemoveAt(index);
                structures[structure].elements.RemoveAt(index);
            }
        }
    }

    private void getGeometry(ref GCDrawingfield drawing_, int cell, int element, List<string> hashList, int cellIndex, string searchString)
    {
        getGeometry(drawing_.cellList[cell], element, hashList, cellIndex, searchString);
    }

    private GeoData getGeometry_simple(GCCell gcCell, int element, List<string> hashList, int cellIndex)
    {
        List<List<GeoLibPointF>> ret = new();
        List<string> lds = new();
        List<GCPolygon> lp = gcCell.elementList[element].convertToPolygons();
        List<bool> isText = new();
        List<string> names = new();
        foreach (GCPolygon p in lp)
        {
            isText.Add(p.isText());
            names.Add(p.getName());

            string ldString = "L" + p.layer_nr + "D" + p.datatype_nr;

            // We should remove identical polygons here in case of doubled-up input geometry.
            string crP_Hash = utility.Utils.GetMD5Hash(p.pointarray);

            switch (hashList.IndexOf(crP_Hash))
            {
                case -1:
                {
                    hashList.Add(crP_Hash);
                    List<GeoLibPointF> t = new();
                    foreach (GeoLibPoint t1 in p.pointarray)
                    {
                        t.Add(new GeoLibPointF(t1.X * scaling, t1.Y * scaling));
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
                    break;
                }
            }
        }

        return new GeoData(ret, isText, lds, names);
    }

    private GeoData getGeometry_complex(GCCell gcCell, int element, List<string> hashList, int cellIndex)
    {
        GeoData ret = new();
        // Need to de-reference these cases.
        double angle = 0.0f;
        GeoLibPoint point = null;
        double mag = 1.0f;
        GCCell tmpCel = gcCell;
        double xSpace = 0;
        double ySpace = 0;
        int xCount = 1;
        int yCount = 1;

        try
        {
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
        }
        catch (Exception)
        {
        }

        if (tmpCel != null) // guard against broken cellref
        {
            for (int referenceElement = 0; referenceElement < tmpCel.elementList.Count; referenceElement++)
            {
                // ReSharper disable once CheckForReferenceEqualityInstead.1
                // ReSharper disable once CheckForReferenceEqualityInstead.3
                if (!tmpCel.elementList[referenceElement].isCellref() && !tmpCel.elementList[referenceElement].GetType().Equals(typeof(GCCellRefArray)))
                {
                    ret = getGeometry_2(gcCell, element, hashList, tmpCel, cellIndex, referenceElement, point, xCount, yCount, xSpace, ySpace, angle, mag);
                }
                else
                {
                    ret = getGeometry_complex(tmpCel, referenceElement, hashList, cellIndex);
                }
            }
        }

        // Need to construct the array of our geometry here!
        GeoData postArray = new();
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

        return postArray;

        // Might need to track nested  array configurations here to handle recursive settings.
    }


    private GeoData getGeometry_2(GCCell gcCell, int element, List<string> hashList, GCCell tmpCel, int cellIndex, int referenceElement, GeoLibPoint point, int xCount, int yCount, double xSpace, double ySpace, double angle, double mag)
    {
        List<List<GeoLibPointF>> ret = new();
        int crLayer = tmpCel.elementList[referenceElement].layer_nr;
        int crDatatype = tmpCel.elementList[referenceElement].datatype_nr;

        List<string> lds = new();

        List<bool> text = new();

        List<string> names = new();

        // See if our layer/datatype combination is known to us already.

        string crSearchString = "L" + crLayer + "D" + crDatatype;
            
        int crLDIndex = -1;

        try
        {
            crLDIndex = structure_LayerDataTypeList[cellIndex].IndexOf(crSearchString);
        }
        catch (Exception)
        {
        }

        switch (crLDIndex)
        {
            case -1:
            {
                structure_LayerDataTypeList[cellIndex].Add(crSearchString);
                structures[cellIndex].addElement();
                int adIndex = structures[cellIndex].elements.Count - 1;

                // ReSharper disable once CheckForReferenceEqualityInstead.1
                if (gcCell.elementList[element].GetType().Equals(typeof(GCCellRefArray)))
                {
                    structures[cellIndex].elements[adIndex].isCellRefArray = gcCell.elementList[element].depend().cellName;
                }

                GeoLibArray tmpArray = new()
                {
                    count = new GeoLibPoint(xCount, yCount),
                    point = new GeoLibPoint(point),
                    pitch = new GeoLibPoint(xSpace, ySpace)
                };
                structures[cellIndex].elements[adIndex].arrayData = tmpArray;
                break;
            }
        }

        try
        {
            List<GCPolygon> lCRP = tmpCel.elementList[referenceElement].convertToPolygons();
            foreach (GCPolygon crP in lCRP)
            {
                crP.move(point);
                crP.rotate(angle, point);
                crP.scale(point, mag);

                string ldString = "L" + crP.layer_nr + "D" + crP.datatype_nr;

                // We should remove identical polygons here in case of doubled-up input geometry.
                string crP_Hash = utility.Utils.GetMD5Hash(crP.pointarray);

                switch (hashList.IndexOf(crP_Hash))
                {
                    case -1:
                    {
                        hashList.Add(crP_Hash);

                        int x = 0;
                        int y = 0;
                        List<GeoLibPointF> t = new();
                        foreach (GeoLibPoint t1 in crP.pointarray)
                        {
                            t.Add(new GeoLibPointF((t1.X + x * xSpace) * scaling, (t1.Y + y * ySpace) * scaling));
                        }
                        structures[cellIndex].addPoly(t, ldString);
                        ret.Add(t);
                        text.Add(crP.isText());
                        lds.Add(ldString);
                        names.Add(crP.getName());
                        break;
                    }
                }
            }
        }
        catch (Exception)
        {
            // Exception is not a big deal.
        }

        return new GeoData(ret, text, lds, names);
    }


    public class GeoData
    {
        public List<List<GeoLibPointF>> geo;
        public List<string> ld;
        public List<bool> isText;
        public List<string> name;

        public GeoData()
        {
            geo = new List<List<GeoLibPointF>>();
            ld = new List<string>();
            isText = new List<bool>();
            name = new List<string>();
        }

        public GeoData(List<List<GeoLibPointF>> poly, List<bool> text, List<string> LDs, List<string> names)
        {
            geo = poly.ToList();
            ld = LDs.ToList();
            isText = text.ToList();
            name = names.ToList();
        }
    }

    private void getGeometry(GCCell gcCell, int element, List<string> hashList, int cellIndex, string searchString)
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
            
        for (int i = 0; i < ret.geo.Count; i++)
        {
            structures[cellIndex].addBakedGeo(ret.geo[i].ToArray(), ret.isText[i], ret.ld[i]);
        }
    }

    public List<GeoLibArray> getArrayParameters()
    {
        return pGetArrayParameters();
    }

    private List<GeoLibArray> pGetArrayParameters()
    {
        return structures[activeStructure].getArrayData(pActiveStructure_LDList[activeLD]);
    }

    public List<bool> isText()
    {
        return pIsText();
    }

    private List<bool> pIsText()
    {
        //string text = pActiveStructure_LDList[activeLD];
        //return structures[activeStructure].isText(text);

        List<bool> ret = new();

        List<GCPolygon> tmp = convertToPolygons(true);
        foreach (GCPolygon t in tmp)
        {
            ret.Add(t.isText());
        }

        return ret;

    }

    public List<string> names()
    {
        return pNames();
    }

    private List<string> pNames()
    {
        //string text = pActiveStructure_LDList[activeLD];
        //return structures[activeStructure].isText(text);

        List<string> ret = new();

        List<GCPolygon> tmp = convertToPolygons(true);
        foreach (GCPolygon t in tmp)
        {
            ret.Add(t.getName());
        }

        return ret;

    }

    public List<GeoLibPointF[]> points(bool flatten)
    {
        return pPoints(flatten);
    }

    private List<GeoLibPointF[]> pPoints(bool flatten)
    {
        switch (flatten)
        {
            case true:
            {
                List<GeoLibPointF[]> ret = new();
                List<GCPolygon> tmp = convertToPolygons(true);
                foreach (GCPolygon t in tmp)
                {
                    GeoLibPointF[] p = GeoWrangler.pointFsFromPoint(t.pointarray, 1);
                    ret.Add(p);
                }
                return ret;
            }
        }

        // Do we ever get here?

        List<GeoLibPointF[]> points = new();
        List<GeoLibPoint> array_count = new();
        List<GeoLibPointF> array_pitch = new();

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
            
        return points;
    }

    public string getName(int index)
    {
        return pGetName(index);
    }

    private string pGetName(int index)
    {
        return structures[activeStructure].elements[index].name;
    }

    public void updateGeometry(int structure, int layerdatatype)
    {
        pUpdateGeometry(structure, layerdatatype);
    }

    private void pUpdateGeometry(int structure, int layerdatatype)
    {
        switch (valid)
        {
            case true:
                drawingField.active_cell = structure;
                activeStructure = structure;
                activeLD = layerdatatype;
                genLDList();
                break;
        }
    }

    public List<GCPolygon> convertToPolygons(int layer, int datatype)
    {
        return pConvertToPolygons(layer: layer, datatype: datatype);
    }

    public List<GCPolygon> convertToPolygons(bool activeLDOnly = false)
    {
        return pConvertToPolygons(activeLDOnly);
    }

    private List<GCPolygon> pConvertToPolygons(bool activeLDOnly = false)
    {
        int layer = -1;
        int datatype = -1;
        switch (activeLDOnly)
        {
            case true:
            {
                // Recover our layer and datatype from the string representation.
                string ld = layerNames.FirstOrDefault(x => x.Value == pActiveStructure_LDList[activeLD]).Key;
                string[] temp = ld.Split(new [] { 'L' })[1].Split(new [] { 'D' });
                layer = Convert.ToInt32(temp[0]);
                datatype = Convert.ToInt32(temp[1]);
                break;
            }
        }
        return pConvertToPolygons(layer: layer, datatype: datatype);
    }

    private List<GCPolygon> pConvertToPolygons(int layer = -1, int datatype = -1)
    {
        return drawingField.cellList[activeStructure].convertToPolygons(layer: layer, datatype: datatype);
    }

    public bool nestedCellRef(int cellIndex)
    {
        bool ret = false;

        for (int i = 0; i < drawingField.cellList[activeStructure].elementList.Count; i++)
        {
            ret = nestedCellRef(cellIndex, i);

            if (ret)
            {
                break;
            }
        }

        return ret;
    }

    public bool nestedCellRef(int cellIndex, int elementIndex)
    {
        return nestedCellRef(drawingField.cellList[activeStructure], elementIndex);
    }

    private bool nestedCellRef(GCCell cell, int elementIndex)
    {
        bool ret = false;
        if (cell.elementList[elementIndex].isCellref() || cell.elementList[elementIndex].isCellrefArray())
        {
            GCCell rCell = cell.elementList[elementIndex].getCellref();
            foreach (GCElement t in rCell.elementList)
            {
                if (t.isCellref() || t.isCellrefArray())
                {
                    ret = true;
                    break;
                }
            }
        }

        return ret;
    }

}