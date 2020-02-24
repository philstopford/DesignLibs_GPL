using Eto.Drawing;
using Eto.Forms;
using Eto.OpenTK;
using OpenTK;
using OpenTK.Graphics.OpenGL;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace etoViewport
{
    public class TestViewport : GLSurface
    {
        public delegate void updateHost();
        public updateHost updateHostFunc { get; set; }

        public bool ok { get; set; }
        bool immediateMode;
        public bool savedLocation_valid { get; set; }
        PointF savedLocation;

        Vector3[] polyArray;
        Vector4[] polyColorArray;
        int[] first;
        int[] count;
        int poly_vbo_size;

        Vector3[] tessPolyArray;
        Vector4[] tessPolyColorArray;
        int[] tessFirst;
        int[] tessCount;
        int tessPoly_vbo_size;

        Vector3[] lineArray;
        Vector4[] lineColorArray;
        int[] lineFirst;
        int[] lineCount;
        int line_vbo_size;

        Vector3[] gridArray;
        Vector3[] gridColorArray;
        int grid_vbo_size;

        Vector3[] axesArray;
        Vector3[] axesColorArray;
        int axes_vbo_size;

        public OVPSettings ovpSettings { get; set; } // note that this is a reference to the real settings.
        float axisZ;
        float gridZ;

        // Use for drag handling.
        bool dragging;
        float x_orig;
        float y_orig;

        ContextMenu menu;

        object polyLock = new object();

        Point WorldToScreen(float x, float y)
        {
            return new Point((int)((x - ovpSettings.cameraPosition.X / (ovpSettings.zoomFactor * ovpSettings.base_zoom)) + Width / 2),
                    (int)((y - ovpSettings.cameraPosition.Y / (ovpSettings.zoomFactor * ovpSettings.base_zoom)) + Height / 2));
        }

        Point WorldToScreen(PointF pt)
        {
            return WorldToScreen(pt.X, pt.Y);
        }

        Size WorldToScreen(SizeF pt)
        {
            Point pt1 = WorldToScreen(0, 0);
            Point pt2 = WorldToScreen(pt.Width, pt.Height);
            return new Size(pt2.X - pt1.X, pt2.Y - pt1.Y);
        }

        PointF ScreenToWorld(int x, int y)
        {
            return new PointF((x - Width / 2) * (ovpSettings.zoomFactor * ovpSettings.base_zoom) + ovpSettings.cameraPosition.X,
                     (y - Height / 2) * (ovpSettings.zoomFactor * ovpSettings.base_zoom) + ovpSettings.cameraPosition.Y);
        }

        PointF ScreenToWorld(Point pt)
        {
            return ScreenToWorld(pt.X, pt.Y);
        }

        RectangleF getViewPort()
        {
            PointF bl = ScreenToWorld(Location.X - Width / 2, Location.Y - Height / 2);
            PointF tr = ScreenToWorld(Location.X + Width / 2, Location.Y + Height / 2);
            return new RectangleF(bl.X, bl.Y, tr.X - bl.X, tr.Y - bl.Y);
        }

        void setViewPort(float x1, float y1, float x2, float y2)
        {
            float h = Math.Abs(y1 - y2);
            float w = Math.Abs(x1 - x2);
            ovpSettings.cameraPosition = new PointF((x1 + x2) / 2, (y1 + y2) / 2);
            if ((Height != 0) && (Width != 0))
            {
                ovpSettings.zoomFactor = Math.Max(h / Height, w / Width);
            }
            else
            {
                ovpSettings.zoomFactor = 1;
            }
        }

        void downHandler(object sender, MouseEventArgs e)
        {
            if (e.Buttons == MouseButtons.Primary)
            {
                if (!dragging && !ovpSettings.lockedViewport) // might not be needed, but seemed like a safe approach to avoid re-setting these in a drag event.
                {
                    x_orig = e.Location.X;
                    y_orig = e.Location.Y;
                    dragging = true;
                }
            }
            //e.Handled = true;
        }

        public void saveLocation()
        {
            savedLocation = new PointF(ovpSettings.cameraPosition.X, ovpSettings.cameraPosition.Y);
            savedLocation_valid = true;
        }

        public void zoomExtents()
        {
            getExtents();

            if (((ovpSettings.polyList.Count == 0) && (ovpSettings.bgPolyList.Count == 0 && (ovpSettings.lineList.Count == 0))) ||
                ((ovpSettings.minX == 0) && (ovpSettings.maxX == 0)) ||
                ((ovpSettings.minY == 0) && (ovpSettings.maxY == 0)))
            {
                reset();
                return;
            }

            // Locate camera at center of the polygon field.
            float dX = ovpSettings.maxX - ovpSettings.minX;
            float dY = ovpSettings.maxY - ovpSettings.minY;
            float cX = (dX / 2.0f) + ovpSettings.minX;
            float cY = (dY / 2.0f) + ovpSettings.minY;

            // Now need to get the zoom level organized.
            float zoomLevel_x = dX / Width;
            float zoomLevel_y = dY / Height;

            if (zoomLevel_x > zoomLevel_y)
            {
                ovpSettings.zoomFactor = zoomLevel_x / ovpSettings.base_zoom;
            }
            else
            {
                ovpSettings.zoomFactor = zoomLevel_y / ovpSettings.base_zoom;
            }

            goToLocation(cX, cY);
        }

        public void loadLocation()
        {
            if (savedLocation_valid)
            {
                ovpSettings.cameraPosition = new PointF(savedLocation.X, savedLocation.Y);
                updateViewport();
            }
        }

        public void goToLocation(float x, float y)
        {
            ovpSettings.cameraPosition = new PointF(x, y);
            updateViewport();
        }

        void dragHandler(object sender, MouseEventArgs e)
        {
            if (ovpSettings.lockedViewport)
            {
                return;
            }
            if (e.Buttons == MouseButtons.Primary)
            {
                object locking = new object();
                lock (locking)
                {
                    // Scaling factor is arbitrary - just based on testing to avoid insane panning speeds.
                    float new_X = ovpSettings.cameraPosition.X - ((e.Location.X - x_orig) * (ovpSettings.base_zoom * ovpSettings.zoomFactor));
                    float new_Y = ovpSettings.cameraPosition.Y + ((e.Location.Y - y_orig) * (ovpSettings.base_zoom * ovpSettings.zoomFactor));
                    ovpSettings.cameraPosition = new PointF(new_X, new_Y);
                    x_orig = e.Location.X;
                    y_orig = e.Location.Y;
                }
            }
            //e.Handled = true;
            updateViewport();
        }

        public void freeze_thaw()
        {
            ovpSettings.lockedViewport = !ovpSettings.lockedViewport;
            updateHostFunc?.Invoke();
        }

        void upHandler(object sender, MouseEventArgs e)
        {
            if (e.Buttons == MouseButtons.Alternate)
            {
                menu.Show(this);
            }
            if (ovpSettings.lockedViewport)
            {
                return;
            }
            if (e.Buttons == MouseButtons.Primary)
            {
                dragging = false;
            }
            //e.Handled = true
        }

        public void zoomIn(float delta)
        {
            if (ovpSettings.lockedViewport)
            {
                return;
            }
            ovpSettings.zoomFactor += (ovpSettings.zoomStep * 0.01f * delta);
        }

        public void zoomOut(float delta)
        {
            if (ovpSettings.lockedViewport)
            {
                return;
            }
            ovpSettings.zoomFactor -= (ovpSettings.zoomStep * 0.01f * delta);
            if (ovpSettings.zoomFactor < 0.0001)
            {
                ovpSettings.zoomFactor = 0.0001f; // avoid any chance of getting to zero.
            }
        }

        void panVertical(float delta)
        {
            if (ovpSettings.lockedViewport)
            {
                return;
            }
            ovpSettings.cameraPosition = new PointF(ovpSettings.cameraPosition.Y, ovpSettings.cameraPosition.Y + delta);// / 10;
        }

        void panHorizontal(float delta)
        {
            if (ovpSettings.lockedViewport)
            {
                return;
            }
            ovpSettings.cameraPosition = new PointF(ovpSettings.cameraPosition.X, ovpSettings.cameraPosition.X + delta);// / 10;
        }

        void addKeyHandler(object sender, EventArgs e)
        {
            KeyDown += keyHandler;
        }

        void removeKeyHandler(object sender, EventArgs e)
        {
            KeyDown -= keyHandler;
        }

        public void reset()
        {
            if (ovpSettings.lockedViewport)
            {
                return;
            }
            ovpSettings.cameraPosition = new PointF(ovpSettings.default_cameraPosition.X, ovpSettings.default_cameraPosition.Y);
            ovpSettings.zoomFactor = 1.0f;
        }

        void keyHandler(object sender, KeyEventArgs e)
        {
            if (ovpSettings.lockedViewport)
            {
                if (e.Key != Keys.F)
                {
                    return;
                }
                ovpSettings.lockedViewport = false;
                return;
            }

            if (e.Key == Keys.F)
            {
                ovpSettings.lockedViewport = true;
                return;
            }

            if (e.Key == Keys.R)
            {
                reset();
            }

            float stepping = 10.0f * (ovpSettings.zoomFactor * ovpSettings.base_zoom);

            bool doUpdate = true;
            if (e.Key == Keys.A)
            {
                panHorizontal(-stepping);
            }
            if (e.Key == Keys.D)
            {
                panHorizontal(stepping);
            }
            if (e.Key == Keys.W)
            {
                panVertical(stepping);
            }
            if (e.Key == Keys.S)
            {
                panVertical(-stepping);
            }
            if (e.Key == Keys.N)
            {
                zoomOut(-1);
            }
            if (e.Key == Keys.M)
            {
                zoomIn(-1);
            }

            if (e.Key == Keys.X)
            {
                zoomExtents();
                doUpdate = false; // update performed in extents
            }

            if (doUpdate)
            {
                updateViewport();
            }
            e.Handled = true;
        }

        void zoomHandler(object sender, MouseEventArgs e)
        {
            if (ovpSettings.lockedViewport)
            {
                return;
            }

            float wheelZoom = e.Delta.Height; // SystemInformation.MouseWheelScrollLines;
            if (wheelZoom > 0)
            {
                zoomIn(wheelZoom);
            }
            if (wheelZoom < 0)
            {
                zoomOut(-wheelZoom);
            }
            updateViewport();
            //e.Handled = true;
        }

        public void updateViewport()
        {
            if (immediateMode)
            {
                _updateVP_immediate();
            }
            else
            {
                try
                {
                    _updateVP_VBO();
                }
                catch (Exception)
                {
                    // Fallback in case VBO support blows up.
                    immediateMode = true;
                    _updateVP_immediate();
                }
            }
        }

        void _updateVP_immediate()
        {
            MakeCurrent();
            init();
            drawGrid_immediate();
            drawAxes_immediate();
            drawPolygons_immediate();
            if (ovpSettings.showDrawn)
            {
                drawLines_immediate();
            }
            SwapBuffers();
        }

        // Need this to handle the OpenTK memory violation if VBO isn't supported. Without this, the exception is managed by the runtime and the tool crashes.
        // We can, however, handle this gracefully.
        [System.Runtime.ExceptionServices.HandleProcessCorruptedStateExceptions]
        void _updateVP_VBO()
        {
            if (!IsInitialized)
                return;
            try
            {
                init();
                try
                {
                    drawGrid_VBO();
                    drawAxes_VBO();
                    drawPolygons_VBO();
                    if (ovpSettings.showDrawn)
                    {
                        drawLines_VBO();
                    }
                    else
                    {
                        lineArray = null;
                        lineColorArray = null;
                    }
                }
                catch (Exception)
                {
                    throw new Exception("VBO had an issue. Aborting.");
                }

                // Fix in case of nulls
                if (polyArray == null)
                {
                    polyArray = new Vector3[2];
                    polyColorArray = new Vector4[polyArray.Length];
                    for (int i = 0; i < polyArray.Length; i++)
                    {
                        polyArray[i] = new Vector3(0.0f);
                        polyColorArray[i] = new Vector4(1.0f);
                    }
                }
                if (tessPolyArray == null)
                {
                    tessPolyArray = new Vector3[2];
                    tessPolyColorArray = new Vector4[tessPolyArray.Length];
                    for (int i = 0; i < tessPolyArray.Length; i++)
                    {
                        tessPolyArray[i] = new Vector3(0.0f);
                        tessPolyColorArray[i] = new Vector4(1.0f);
                    }
                }
                if (lineArray == null)
                {
                    lineArray = new Vector3[2];
                    lineColorArray = new Vector4[lineArray.Length];
                    for (int i = 0; i < lineArray.Length; i++)
                    {
                        lineArray[i] = new Vector3(0.0f);
                        lineColorArray[i] = new Vector4(1.0f);
                    }
                }


                // Now we wrangle our VBOs
                grid_vbo_size = gridArray.Length; // Necessary for rendering later on
                axes_vbo_size = axesArray.Length; // Necessary for rendering later on
                poly_vbo_size = polyArray.Length; // Necessary for rendering later on
                tessPoly_vbo_size = tessPolyArray.Length;
                line_vbo_size = lineArray.Length;

                int numBuffers = 5;
                int[] vbo_id = new int[numBuffers];
                GL.GenBuffers(numBuffers, vbo_id);

                int[] col_id = new int[numBuffers];
                GL.GenBuffers(numBuffers, col_id);

                GL.BindBuffer(BufferTarget.ArrayBuffer, vbo_id[0]);
                GL.BufferData(BufferTarget.ArrayBuffer,
                          new IntPtr(gridArray.Length * BlittableValueType.StrideOf(gridArray)),
                          gridArray, BufferUsageHint.StaticDraw);

                GL.BindBuffer(BufferTarget.ArrayBuffer, col_id[0]);
                GL.BufferData(BufferTarget.ArrayBuffer,
                          new IntPtr(gridColorArray.Length * BlittableValueType.StrideOf(gridColorArray)),
                          gridColorArray, BufferUsageHint.StaticDraw);

                GL.BindBuffer(BufferTarget.ArrayBuffer, vbo_id[1]);
                GL.BufferData(BufferTarget.ArrayBuffer,
                          new IntPtr(axesArray.Length * BlittableValueType.StrideOf(axesArray)),
                          axesArray, BufferUsageHint.StaticDraw);

                GL.BindBuffer(BufferTarget.ArrayBuffer, col_id[1]);
                GL.BufferData(BufferTarget.ArrayBuffer,
                          new IntPtr(axesColorArray.Length * BlittableValueType.StrideOf(axesColorArray)),
                          axesColorArray, BufferUsageHint.StaticDraw);

                GL.BindBuffer(BufferTarget.ArrayBuffer, vbo_id[2]);
                GL.BufferData(BufferTarget.ArrayBuffer,
                          new IntPtr(polyArray.Length * BlittableValueType.StrideOf(polyArray)),
                          polyArray, BufferUsageHint.StaticDraw);

                GL.BindBuffer(BufferTarget.ArrayBuffer, col_id[2]);
                GL.BufferData(BufferTarget.ArrayBuffer,
                          new IntPtr(polyColorArray.Length * BlittableValueType.StrideOf(polyColorArray)),
                          polyColorArray, BufferUsageHint.StaticDraw);

                GL.BindBuffer(BufferTarget.ArrayBuffer, vbo_id[3]);
                GL.BufferData(BufferTarget.ArrayBuffer,
                          new IntPtr(lineArray.Length * BlittableValueType.StrideOf(lineArray)),
                          lineArray, BufferUsageHint.StaticDraw);

                GL.BindBuffer(BufferTarget.ArrayBuffer, col_id[3]);
                GL.BufferData(BufferTarget.ArrayBuffer,
                          new IntPtr(lineColorArray.Length * BlittableValueType.StrideOf(lineColorArray)),
                          lineColorArray, BufferUsageHint.StaticDraw);

                GL.BindBuffer(BufferTarget.ArrayBuffer, vbo_id[4]);
                GL.BufferData(BufferTarget.ArrayBuffer,
                          new IntPtr(tessPolyArray.Length * BlittableValueType.StrideOf(tessPolyArray)),
                          tessPolyArray, BufferUsageHint.StaticDraw);

                GL.BindBuffer(BufferTarget.ArrayBuffer, col_id[4]);
                GL.BufferData(BufferTarget.ArrayBuffer,
                          new IntPtr(tessPolyColorArray.Length * BlittableValueType.StrideOf(tessPolyColorArray)),
                          tessPolyColorArray, BufferUsageHint.StaticDraw);

                // To draw a VBO:
                // 1) Ensure that the VertexArray client state is enabled.
                // 2) Bind the vertex and element buffer handles.
                // 3) Set up the data pointers(vertex, normal, color) according to your vertex format.


                try
                {
                    if (ovpSettings.antiAlias)
                    {
                        GL.Enable(EnableCap.Multisample);
                        if (ovpSettings.drawPoints)
                        {
                            GL.Enable(EnableCap.PointSmooth); // should result in circles rather than squares. We shall see.
                        }
                        //GL.Enable(EnableCap.LineSmooth);
                        //GL.Enable(EnableCap.PolygonSmooth);
                    }

                    GL.EnableClientState(ArrayCap.VertexArray);
                    GL.EnableClientState(ArrayCap.ColorArray);
                    GL.BindBuffer(BufferTarget.ArrayBuffer, vbo_id[0]);
                    GL.VertexPointer(3, VertexPointerType.Float, Vector3.SizeInBytes, new IntPtr(0));
                    GL.BindBuffer(BufferTarget.ArrayBuffer, col_id[0]);
                    GL.ColorPointer(3, ColorPointerType.Float, Vector3.SizeInBytes, new IntPtr(0));
                    GL.DrawArrays(PrimitiveType.Lines, 0, grid_vbo_size);
                    GL.DisableClientState(ArrayCap.VertexArray);
                    GL.DisableClientState(ArrayCap.ColorArray);

                    GL.EnableClientState(ArrayCap.VertexArray);
                    GL.EnableClientState(ArrayCap.ColorArray);
                    GL.BindBuffer(BufferTarget.ArrayBuffer, vbo_id[1]);
                    GL.VertexPointer(3, VertexPointerType.Float, Vector3.SizeInBytes, new IntPtr(0));
                    GL.BindBuffer(BufferTarget.ArrayBuffer, col_id[1]);
                    GL.ColorPointer(3, ColorPointerType.Float, Vector3.SizeInBytes, new IntPtr(0));
                    GL.DrawArrays(PrimitiveType.Lines, 0, axes_vbo_size);
                    GL.DisableClientState(ArrayCap.VertexArray);
                    GL.DisableClientState(ArrayCap.ColorArray);

                    GL.EnableClientState(ArrayCap.VertexArray);
                    GL.EnableClientState(ArrayCap.ColorArray);

                    if (ovpSettings.drawPoints)
                    {
                        GL.PointSize(2.0f);
                    }

                    // Allow alpha blending
                    GL.Enable(EnableCap.Blend);
                    GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha);

                    // Poly data
                    GL.BindBuffer(BufferTarget.ArrayBuffer, vbo_id[2]);
                    GL.VertexPointer(3, VertexPointerType.Float, Vector3.SizeInBytes, new IntPtr(0));
                    GL.BindBuffer(BufferTarget.ArrayBuffer, col_id[2]);
                    GL.ColorPointer(4, ColorPointerType.Float, Vector4.SizeInBytes, new IntPtr(0));
                    // Draw the border.
                    GL.DrawArrays(PrimitiveType.Lines, 0, poly_vbo_size);
                    if (ovpSettings.drawPoints)
                    {
                        GL.DrawArrays(PrimitiveType.Points, 0, poly_vbo_size);
                    }

                    GL.DisableClientState(ArrayCap.VertexArray);
                    GL.DisableClientState(ArrayCap.ColorArray);

                    GL.EnableClientState(ArrayCap.VertexArray);
                    GL.EnableClientState(ArrayCap.ColorArray);
                    // Line data
                    GL.BindBuffer(BufferTarget.ArrayBuffer, vbo_id[3]);
                    GL.VertexPointer(3, VertexPointerType.Float, Vector3.SizeInBytes, new IntPtr(0));
                    GL.BindBuffer(BufferTarget.ArrayBuffer, col_id[3]);
                    GL.ColorPointer(4, ColorPointerType.Float, Vector4.SizeInBytes, new IntPtr(0));
                    GL.DrawArrays(PrimitiveType.Lines, 0, line_vbo_size);

                    GL.Disable(EnableCap.Blend);
                    if (ovpSettings.antiAlias)
                    {
                        GL.Disable(EnableCap.Multisample);
                        //GL.Disable(EnableCap.LineSmooth);
                        //GL.Disable(EnableCap.PolygonSmooth);
                    }
                    GL.DisableClientState(ArrayCap.VertexArray);
                    GL.DisableClientState(ArrayCap.ColorArray);

                    if (tessCount.Length > 0)
                    {
                        GL.EnableClientState(ArrayCap.VertexArray);
                        GL.EnableClientState(ArrayCap.ColorArray);

                        // Allow alpha blending
                        GL.Enable(EnableCap.Blend);
                        GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha);

                        // Tessellation Poly data
                        GL.BindBuffer(BufferTarget.ArrayBuffer, vbo_id[4]);
                        GL.VertexPointer(3, VertexPointerType.Float, Vector3.SizeInBytes, new IntPtr(0));
                        GL.BindBuffer(BufferTarget.ArrayBuffer, col_id[4]);
                        GL.ColorPointer(4, ColorPointerType.Float, Vector4.SizeInBytes, new IntPtr(0));
                        // GL.MultiDrawArrays(PrimitiveType.Triangles, tessFirst, tessCount, tessCount.Length);
                        GL.DrawArrays(PrimitiveType.Triangles, 0, tessPoly_vbo_size);

                        GL.Disable(EnableCap.Blend);
                        if (ovpSettings.antiAlias)
                        {
                            GL.Disable(EnableCap.Multisample);
                            //GL.Disable(EnableCap.LineSmooth);
                            //GL.Disable(EnableCap.PolygonSmooth);
                        }

                        GL.DisableClientState(ArrayCap.VertexArray);
                        GL.DisableClientState(ArrayCap.ColorArray);
                    }
                }
                catch (Exception)
                {
                    throw new Exception("VBO had an issue. Aborting.");
                }

                SwapBuffers();
                GL.Flush();
                GL.DeleteBuffers(numBuffers, vbo_id);
                GL.DeleteBuffers(numBuffers, col_id);
            }
            catch (Exception)
            {
                ok = false;
                throw new Exception("VBO had an issue. Aborting.");
            }
            updateHostFunc?.Invoke();
        }

        void getExtents()
        {
            float minX = 0;
            float maxX = 0;
            float minY = 0, maxY = 0;

            int polyListCount = ovpSettings.polyList.Count;
            int bgPolyListCount = ovpSettings.bgPolyList.Count;
            int lineListCount = ovpSettings.lineList.Count;

            if ((polyListCount == 0) && (bgPolyListCount == 0) && (lineListCount == 0))
            {
                ovpSettings.minX = 0;
                ovpSettings.maxX = 0;
                ovpSettings.minY = 0;
                ovpSettings.maxY = 0;
                return;
            }

            if (polyListCount != 0)
            {
                minX = ovpSettings.polyList[0].poly[0].X;
                maxX = ovpSettings.polyList[0].poly[0].X;
                minY = ovpSettings.polyList[0].poly[0].Y;
                maxY = ovpSettings.polyList[0].poly[0].Y;
                for (int poly = 0; poly < polyListCount; poly++)
                {
                    if (ovpSettings.drawnPoly[poly] && !ovpSettings.showDrawn)
                    {
                        continue;
                    }
                    float tMinX = ovpSettings.polyList[poly].poly.Min(p => p.X);
                    if (tMinX < minX)
                    {
                        minX = tMinX;
                    }
                    float tMaxX = ovpSettings.polyList[poly].poly.Max(p => p.X);
                    if (tMaxX > maxX)
                    {
                        maxX = tMaxX;
                    }
                    float tMinY = ovpSettings.polyList[poly].poly.Min(p => p.Y);
                    if (tMinY < minY)
                    {
                        minY = tMinY;
                    }
                    float tMaxY = ovpSettings.polyList[poly].poly.Max(p => p.Y);
                    if (tMaxY > maxY)
                    {
                        maxY = tMaxY;
                    }
                }
            }

            if (bgPolyListCount != 0)
            {
                if (polyListCount == 0)
                {
                    minX = ovpSettings.bgPolyList[0].poly[0].X;
                    maxX = ovpSettings.bgPolyList[0].poly[0].X;
                    minY = ovpSettings.bgPolyList[0].poly[0].Y;
                    maxY = ovpSettings.bgPolyList[0].poly[0].Y;
                }

                for (int line = 0; line < bgPolyListCount; line++)
                {
                    float tMinX = ovpSettings.bgPolyList[line].poly.Min(p => p.X);
                    if (tMinX < minX)
                    {
                        minX = tMinX;
                    }
                    float tMaxX = ovpSettings.bgPolyList[line].poly.Max(p => p.X);
                    if (tMaxX > maxX)
                    {
                        maxX = tMaxX;
                    }
                    float tMinY = ovpSettings.bgPolyList[line].poly.Min(p => p.Y);
                    if (tMinY < minY)
                    {
                        minY = tMinY;
                    }
                    float tMaxY = ovpSettings.bgPolyList[line].poly.Max(p => p.Y);
                    if (tMaxY > maxY)
                    {
                        maxY = tMaxY;
                    }
                }
            }

            if (ovpSettings.showDrawn)
            {
                if (lineListCount != 0)
                {
                    if ((polyListCount == 0) && (bgPolyListCount == 0))
                    {
                        minX = ovpSettings.lineList[0].poly[0].X;
                        maxX = ovpSettings.lineList[0].poly[0].X;
                        minY = ovpSettings.lineList[0].poly[0].Y;
                        maxY = ovpSettings.lineList[0].poly[0].Y;
                    }

                    for (int line = 0; line < lineListCount; line++)
                    {
                        float tMinX = ovpSettings.lineList[line].poly.Min(p => p.X);
                        if (tMinX < minX)
                        {
                            minX = tMinX;
                        }
                        float tMaxX = ovpSettings.lineList[line].poly.Max(p => p.X);
                        if (tMaxX > maxX)
                        {
                            maxX = tMaxX;
                        }
                        float tMinY = ovpSettings.lineList[line].poly.Min(p => p.Y);
                        if (tMinY < minY)
                        {
                            minY = tMinY;
                        }
                        float tMaxY = ovpSettings.lineList[line].poly.Max(p => p.Y);
                        if (tMaxY > maxY)
                        {
                            maxY = tMaxY;
                        }
                    }
                }
            }
            ovpSettings.minX = minX;
            ovpSettings.maxX = maxX;
            ovpSettings.minY = minY;
            ovpSettings.maxY = maxY;
        }

        void drawPolygons_VBO()
        {
            try
            {
                List<Vector3> polyList = new List<Vector3>();
                List<Vector4> polyColorList = new List<Vector4>();

                List<Vector3> tessPolyList = new List<Vector3>();
                List<Vector4> tessPolyColorList = new List<Vector4>();

                int polyListCount = ovpSettings.polyList.Count();
                int bgPolyListCount = ovpSettings.bgPolyList.Count();
                int tessPolyListCount = ovpSettings.tessPolyList.Count();

                // Carve our Z-space up to stack polygons
                int numPolys = 1;

                numPolys = polyListCount + bgPolyListCount;
                if (ovpSettings.enableFilledPolys)
                {
                    numPolys += tessPolyListCount;
                }

                float polyZStep = 1.0f / Math.Max(1, numPolys + 1); // avoid a div by zero risk; pad the poly number also to reduce risk of adding a poly beyond the clipping range

                // Create our first and count arrays for the vertex indices, to enable polygon separation when rendering.
                first = new int[numPolys];
                count = new int[numPolys];

                tessFirst = new int[tessPolyListCount];
                tessCount = new int[tessPolyListCount];

                int counter = 0; // vertex count that will be used to define 'first' index for each polygon.
                int previouscounter = 0; // will be used to derive the number of vertices in each polygon.

                float polyZ = 0;

                if (ovpSettings.enableFilledPolys)
                {
                    for (int i = 0; i < tessPolyListCount; i++)
                    {
                        tessFirst[i] = i * 3;
                        float alpha = ovpSettings.tessPolyList[i].alpha;
                        polyZ += polyZStep;
                        for (int j = 0; j < 3; j++)
                        {
                            tessPolyList.Add(new Vector3(ovpSettings.tessPolyList[i].poly[j].X, ovpSettings.tessPolyList[i].poly[j].Y, polyZ));
                            tessPolyColorList.Add(new Vector4(ovpSettings.tessPolyList[i].color.R, ovpSettings.tessPolyList[i].color.G, ovpSettings.tessPolyList[i].color.B, alpha));
                        }
                        tessCount[i] = 3;
                    }
                }

                // Pondering options here - this would make a nice border construct around the filled geometry, amongst other things.
                for (int poly = 0; poly < polyListCount; poly++)
                {
                    float alpha = ovpSettings.polyList[poly].alpha;
                    if (ovpSettings.enableFilledPolys)
                    {
                        alpha = 1.0f;
                    }
                    polyZ += polyZStep;
                    first[poly] = counter;
                    previouscounter = counter;
                    int polyLength = ovpSettings.polyList[poly].poly.Length - 1;
                    for (int pt = 0; pt < polyLength; pt++)
                    {
                        polyList.Add(new Vector3(ovpSettings.polyList[poly].poly[pt].X, ovpSettings.polyList[poly].poly[pt].Y, polyZ));
                        counter++;
                        polyColorList.Add(new Vector4(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, alpha));
                        polyList.Add(new Vector3(ovpSettings.polyList[poly].poly[pt + 1].X, ovpSettings.polyList[poly].poly[pt + 1].Y, polyZ));
                        counter++;
                        polyColorList.Add(new Vector4(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, alpha));
                    }
                    count[poly] = counter - previouscounter; // set our vertex count for the polygon.
                }

                polyZ = 0;
                for (int poly = 0; poly < bgPolyListCount; poly++)
                {
                    float alpha = ovpSettings.bgPolyList[poly].alpha;
                    polyZ += polyZStep;
                    first[poly] = counter;
                    previouscounter = counter;

                    int bgPolyLength = ovpSettings.bgPolyList[poly].poly.Length - 1;
                    for (int pt = 0; pt < bgPolyLength; pt++)
                    {
                        polyList.Add(new Vector3(ovpSettings.bgPolyList[poly].poly[pt].X, ovpSettings.bgPolyList[poly].poly[pt].Y, polyZ));
                        counter++;
                        polyColorList.Add(new Vector4(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G, ovpSettings.bgPolyList[poly].color.B, alpha));
                        polyList.Add(new Vector3(ovpSettings.bgPolyList[poly].poly[pt + 1].X, ovpSettings.bgPolyList[poly].poly[pt + 1].Y, polyZ));
                        counter++;
                        polyColorList.Add(new Vector4(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G, ovpSettings.bgPolyList[poly].color.B, alpha));
                    }
                    count[poly + polyListCount] = counter - previouscounter; // set our vertex count for the polygon.
                }

                polyArray = polyList.ToArray();
                polyColorArray = polyColorList.ToArray();

                tessPolyArray = tessPolyList.ToArray();
                tessPolyColorArray = tessPolyColorList.ToArray();
            }
            catch (Exception)
            {
                // Can ignore - not critical.
            }
        }

        void drawPolygons_immediate()
        {
            MakeCurrent();
            try
            {
                GL.LoadIdentity();

                int polyListCount = ovpSettings.polyList.Count;
                int tessPolyListCount = ovpSettings.tessPolyList.Count;
                int bgPolyListCount = ovpSettings.bgPolyList.Count;

                if (ovpSettings.enableFilledPolys)
                {
                    int numPolys = tessPolyListCount + bgPolyListCount;
                    float polyZ = 1.0f / (numPolys + 1); // push our filled polygons behind the boundary
                    for (int poly = 0; poly < tessPolyListCount; poly++)
                    {
                        GL.Color4(ovpSettings.tessPolyList[poly].color.R, ovpSettings.tessPolyList[poly].color.G, ovpSettings.tessPolyList[poly].color.B, ovpSettings.tessPolyList[poly].alpha);
                        GL.Enable(EnableCap.Blend);
                        GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha);
                        GL.Begin(PrimitiveType.Triangles);
                        GL.Vertex3(new Vector3(ovpSettings.tessPolyList[poly].poly[0].X, ovpSettings.tessPolyList[poly].poly[0].Y, polyZ));
                        GL.Vertex3(new Vector3(ovpSettings.tessPolyList[poly].poly[1].X, ovpSettings.tessPolyList[poly].poly[1].Y, polyZ));
                        GL.Vertex3(new Vector3(ovpSettings.tessPolyList[poly].poly[2].X, ovpSettings.tessPolyList[poly].poly[2].Y, polyZ));
                        GL.End();
                    }
                    for (int poly = 0; poly < bgPolyListCount; poly++)
                    {
                        GL.Color4(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G, ovpSettings.bgPolyList[poly].color.B, ovpSettings.bgPolyList[poly].alpha);
                        GL.Enable(EnableCap.Blend);
                        GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha);
                        GL.Begin(PrimitiveType.Triangles);
                        GL.Vertex3(new Vector3(ovpSettings.bgPolyList[poly].poly[0].X, ovpSettings.bgPolyList[poly].poly[0].Y, polyZ));
                        GL.Vertex3(new Vector3(ovpSettings.bgPolyList[poly].poly[1].X, ovpSettings.bgPolyList[poly].poly[1].Y, polyZ));
                        GL.Vertex3(new Vector3(ovpSettings.bgPolyList[poly].poly[2].X, ovpSettings.bgPolyList[poly].poly[2].Y, polyZ));
                        GL.End();
                    }
                }
                else
                {
                    int numPolys = polyListCount + bgPolyListCount;
                    // Carve our Z-space up to stack polygons
                    float polyZStep = 1.0f / numPolys;
                    float polyZ = 0;
                    for (int poly = 0; poly < polyListCount; poly++)
                    {
                        int polyLength = ovpSettings.polyList[poly].poly.Length - 1;
                        GL.Color4(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G, ovpSettings.polyList[poly].color.B, ovpSettings.polyList[poly].alpha);
                        GL.Begin(PrimitiveType.Lines);
                        for (int pt = 0; pt < polyLength; pt++)
                        {
                            GL.Vertex3(ovpSettings.polyList[poly].poly[pt].X, ovpSettings.polyList[poly].poly[pt].Y, polyZ);
                            GL.Vertex3(ovpSettings.polyList[poly].poly[pt + 1].X, ovpSettings.polyList[poly].poly[pt + 1].Y, polyZ);
                        }
                        GL.End();
                        polyZ += polyZStep;
                    }
                    for (int poly = 0; poly < bgPolyListCount; poly++)
                    {
                        int bgPolyLength = ovpSettings.bgPolyList[poly].poly.Length - 1;
                        GL.Color4(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G, ovpSettings.bgPolyList[poly].color.B, ovpSettings.bgPolyList[poly].alpha);
                        GL.Begin(PrimitiveType.Lines);
                        for (int pt = 0; pt < bgPolyLength; pt++)
                        {
                            GL.Vertex3(ovpSettings.bgPolyList[poly].poly[pt].X, ovpSettings.bgPolyList[poly].poly[pt].Y, polyZ);
                            GL.Vertex3(ovpSettings.bgPolyList[poly].poly[pt + 1].X, ovpSettings.bgPolyList[poly].poly[pt + 1].Y, polyZ);
                        }
                        GL.End();
                        polyZ += polyZStep;
                    }
                }
            }
            catch (Exception)
            {
                // Can ignore - not critical.
            }
        }

        void drawLines_VBO()
        {
            try
            {
                List<Vector3> polyList = new List<Vector3>();
                List<Vector4> polyColorList = new List<Vector4>();

                int lineListCount = ovpSettings.lineList.Count;

                // Carve our Z-space up to stack polygons
                float polyZStep = 1.0f / lineListCount;

                // Create our first and count arrays for the vertex indices, to enable polygon separation when rendering.
                lineFirst = new int[lineListCount];
                lineCount = new int[lineListCount];
                int counter = 0; // vertex count that will be used to define 'first' index for each polygon.
                int previouscounter = 0; // will be used to derive the number of vertices in each polygon.

                for (int poly = 0; poly < lineListCount; poly++)
                {
                    int lineListPolyLength = ovpSettings.lineList[poly].poly.Length - 1;
                    float alpha = ovpSettings.lineList[poly].alpha;
                    float polyZ = poly * polyZStep;
                    lineFirst[poly] = counter;
                    previouscounter = counter;
                    for (int pt = 0; pt < lineListPolyLength; pt++)
                    {
                        polyList.Add(new Vector3(ovpSettings.lineList[poly].poly[pt].X, ovpSettings.lineList[poly].poly[pt].Y, polyZ));
                        counter++;
                        polyColorList.Add(new Vector4(ovpSettings.lineList[poly].color.R, ovpSettings.lineList[poly].color.G, ovpSettings.lineList[poly].color.B, alpha));
                        polyList.Add(new Vector3(ovpSettings.lineList[poly].poly[pt + 1].X, ovpSettings.lineList[poly].poly[pt + 1].Y, polyZ));
                        counter++;
                        polyColorList.Add(new Vector4(ovpSettings.lineList[poly].color.R, ovpSettings.lineList[poly].color.G, ovpSettings.lineList[poly].color.B, alpha));
                    }
                    lineCount[poly] = counter - previouscounter; // set our vertex count for the polygon.
                }

                lineArray = polyList.ToArray();
                lineColorArray = polyColorList.ToArray();
            }
            catch (Exception)
            {
                // Can ignore - not critical.
            }
        }

        void drawLines_VBO_()
        {
            try
            {
                int lineListCount = ovpSettings.lineList.Count;
                int lineListPtCount = ovpSettings.lineListPtCount.Sum();

                lineArray = new Vector3[lineListPtCount * 2];
                lineColorArray = new Vector4[lineListPtCount * 2];

                // Carve our Z-space up to stack polygons
                float polyZStep = 1.0f / lineListCount;

                // Create our first and count arrays for the vertex indices, to enable polygon separation when rendering.
                lineFirst = new int[lineListCount];
                lineCount = new int[lineListCount];

                for (int poly = 0; poly < lineListCount; poly++)
                {
                    int lineListPolyLength = (ovpSettings.lineListPtCount[poly]);
                    float alpha = ovpSettings.lineList[poly].alpha;
                    float polyZ = poly * polyZStep;
                    lineFirst[poly] = ovpSettings.lineListPtCount.Take(poly).Sum() * 2;
                    for (int pt = 0; pt < (lineListPolyLength) / 2; pt++)
                    {
                        lineArray[lineFirst[poly] + (pt * 2)] = new Vector3(ovpSettings.lineList[poly].poly[pt].X, ovpSettings.lineList[poly].poly[pt].Y, polyZ);
                        lineColorArray[lineFirst[poly] + (pt * 2)] = new Vector4(ovpSettings.lineList[poly].color.R, ovpSettings.lineList[poly].color.G, ovpSettings.lineList[poly].color.B, alpha);
                        lineArray[lineFirst[poly] + (pt * 2) + 1] = new Vector3(ovpSettings.lineList[poly].poly[pt + 1].X, ovpSettings.lineList[poly].poly[pt + 1].Y, polyZ);
                        lineColorArray[lineFirst[poly] + (pt * 2) + 1] = new Vector4(ovpSettings.lineList[poly].color.R, ovpSettings.lineList[poly].color.G, ovpSettings.lineList[poly].color.B, alpha);
                    }
                    lineCount[poly] = (ovpSettings.lineListPtCount[poly]) * 2; // set our vertex count for the polygon.
                }
            }
            catch (Exception)
            {
                // Can ignore - not critical.
            }
        }

        void drawLines_immediate()
        {
            try
            {
                // Carve our Z-space up to stack polygons
                int lineListCount = ovpSettings.lineList.Count();
                float polyZStep = 1.0f / lineListCount;

                for (int poly = 0; poly < lineListCount; poly++)
                {
                    int lineListPolyLength = ovpSettings.lineList[poly].poly.Length - 1;
                    float polyZ = poly * polyZStep;
                    GL.Color4(ovpSettings.lineList[poly].color.R, ovpSettings.lineList[poly].color.G, ovpSettings.lineList[poly].color.B, ovpSettings.lineList[poly].alpha);
                    GL.Begin(PrimitiveType.Lines);
                    for (int pt = 0; pt < lineListPolyLength; pt++)
                    {
                        GL.Vertex3(ovpSettings.lineList[poly].poly[pt].X, ovpSettings.lineList[poly].poly[pt].Y, polyZ);
                        GL.Vertex3(ovpSettings.lineList[poly].poly[pt + 1].X, ovpSettings.lineList[poly].poly[pt + 1].Y, polyZ);
                    }
                    GL.End();
                }
            }
            catch (Exception)
            {
                // Can ignore - not critical.
            }
        }

        public void defaults()
        {
            MakeCurrent();
            if (ovpSettings.antiAlias)
            {
                GL.Enable(EnableCap.LineSmooth);
            }
            else
            {
                GL.Disable(EnableCap.LineSmooth);
            }
            GL.Disable(EnableCap.Lighting);
            GL.ShadeModel(ShadingModel.Flat);
            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha);
            GL.PolygonOffset(0.0f, 0.5f);
            GL.LineStipple(1, 61680);
            gridZ = -0.95f;
            axisZ = gridZ + 0.01f;
        }

        public TestViewport(ref OVPSettings svpSettings)
        {
            try
            {
                immediateMode = svpSettings.immediateMode;
                ovpSettings = svpSettings;
                MouseDown += downHandler;
                MouseMove += dragHandler;
                MouseUp += upHandler;
                MouseWheel += zoomHandler;
                GotFocus += addKeyHandler;
                // MouseHover += addKeyHandler;
                LostFocus += removeKeyHandler;
                ok = true;
            }
            catch (Exception)
            {
                //Console.WriteLine($"Error: {ex}");
                ok = false;
            }
        }

        public void setContextMenu(ref ContextMenu menu_)
        {
            menu = menu_;
        }

        public void changeSettingsRef(ref OVPSettings newSettings)
        {
            ovpSettings = newSettings;
            updateViewport();
        }

        protected override void OnDraw(EventArgs e)
        {
            base.OnDraw(e);
            updateViewport();
        }

        public void init()
        {
            MakeCurrent();
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();
            GL.Ortho(ovpSettings.cameraPosition.X - Width * (ovpSettings.zoomFactor * ovpSettings.base_zoom) / 2,
                      ovpSettings.cameraPosition.X + Width * (ovpSettings.zoomFactor * ovpSettings.base_zoom) / 2,
                      ovpSettings.cameraPosition.Y - Height * (ovpSettings.zoomFactor * ovpSettings.base_zoom) / 2,
                      ovpSettings.cameraPosition.Y + Height * (ovpSettings.zoomFactor * ovpSettings.base_zoom) / 2,
                      -1.0f, 1.0f);
            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();
            ovpSettings.bounds = getViewPort();
            GL.ClearColor(ovpSettings.backColor.R, ovpSettings.backColor.G, ovpSettings.backColor.B, ovpSettings.backColor.A);
            GL.ClearDepth(1.0);
            GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);
        }

        public void drawGrid_VBO()
        {
            try
            {
                if (ovpSettings.showGrid)
                {
                    float spacing = ovpSettings.gridSpacing;
                    if (ovpSettings.dynamicGrid)
                    {
                        while (WorldToScreen(new SizeF(spacing, 0.0f)).Width > 12.0f)
                            spacing /= 10.0f;

                        while (WorldToScreen(new SizeF(spacing, 0.0f)).Width < 4.0f)
                            spacing *= 10.0f;
                    }

                    float xLimit = Width * 0.5f;
                    float yLimit = Height * 0.5f;
                    float zoom = ovpSettings.zoomFactor * ovpSettings.base_zoom;

                    if (WorldToScreen(new SizeF(spacing, 0.0f)).Width >= 4.0f)
                    {
                        float nX_end = -(xLimit * zoom) + ovpSettings.cameraPosition.X;
                        int nX_iterations = (int)Math.Abs(Math.Ceiling(nX_end / spacing)) + 1;
                        float pX_end = (xLimit * zoom) + ovpSettings.cameraPosition.X;
                        int pX_iterations = (int)Math.Abs(Math.Ceiling(pX_end / spacing)) + 1;
                        float nY_end = -(yLimit * zoom) + ovpSettings.cameraPosition.Y;
                        int nY_iterations = (int)Math.Abs(Math.Ceiling(nY_end / spacing)) + 1;
                        float pY_end = (yLimit * zoom) + ovpSettings.cameraPosition.Y;
                        int pY_iterations = (int)Math.Abs(Math.Ceiling(pY_end / spacing)) + 1;

                        int lines = 2 * (nX_iterations + pX_iterations + nY_iterations + pY_iterations);

                        gridArray = new Vector3[lines];
                        gridColorArray = new Vector3[lines];

                        ParallelOptions po = new ParallelOptions();

                        Parallel.For(0, nX_iterations, po, (i, loopState) =>
                        {
                            float r = 0.0f;
                            float g = 0.0f;
                            float b = 0.0f;
                            if (i % 10 == 0)
                            {
                                r = ovpSettings.majorGridColor.R;
                                g = ovpSettings.majorGridColor.G;
                                b = ovpSettings.majorGridColor.B;
                            }
                            else
                            {
                                r = ovpSettings.minorGridColor.R;
                                g = ovpSettings.minorGridColor.G;
                                b = ovpSettings.minorGridColor.B;
                            }

                            int index = i * 2;
                            gridArray[index] = new Vector3(-i * spacing, ovpSettings.cameraPosition.Y + (yLimit * zoom), gridZ);
                            gridColorArray[index] = new Vector3(r, g, b);
                            gridArray[index + 1] = new Vector3(-i * spacing, ovpSettings.cameraPosition.Y - (yLimit * zoom), gridZ);
                            gridColorArray[index + 1] = new Vector3(r, g, b);
                        });

                        Parallel.For(0, pX_iterations, po, (i, loopState) =>
                        {
                            float r = 0.0f;
                            float g = 0.0f;
                            float b = 0.0f;
                            if (i % 10 == 0)
                            {
                                r = ovpSettings.majorGridColor.R;
                                g = ovpSettings.majorGridColor.G;
                                b = ovpSettings.majorGridColor.B;
                            }
                            else
                            {
                                r = ovpSettings.minorGridColor.R;
                                g = ovpSettings.minorGridColor.G;
                                b = ovpSettings.minorGridColor.B;
                            }

                            int index = nX_iterations + i;
                            index *= 2;
                            gridArray[index] = new Vector3(i * spacing, ovpSettings.cameraPosition.Y + (zoom * yLimit), gridZ);
                            gridColorArray[index] = new Vector3(r, g, b);
                            gridArray[index + 1] = new Vector3(i * spacing, ovpSettings.cameraPosition.Y + (zoom * -yLimit), gridZ);
                            gridColorArray[index + 1] = new Vector3(r, g, b);
                        });

                        Parallel.For(0, nY_iterations, po, (i, loopState) =>
                        {
                            float r = 0.0f;
                            float g = 0.0f;
                            float b = 0.0f;
                            if (i % 10 == 0)
                            {
                                r = ovpSettings.majorGridColor.R;
                                g = ovpSettings.majorGridColor.G;
                                b = ovpSettings.majorGridColor.B;
                            }
                            else
                            {
                                r = ovpSettings.minorGridColor.R;
                                g = ovpSettings.minorGridColor.G;
                                b = ovpSettings.minorGridColor.B;
                            }

                            int index = nX_iterations + pX_iterations + i;
                            index *= 2;
                            gridArray[index] = new Vector3(ovpSettings.cameraPosition.X + (zoom * xLimit), -i * spacing, gridZ);
                            gridColorArray[index] = new Vector3(r, g, b);
                            gridArray[index + 1] = new Vector3(ovpSettings.cameraPosition.X + (zoom * -xLimit), -i * spacing, gridZ);
                            gridColorArray[index + 1] = new Vector3(r, g, b);
                        });

                        Parallel.For(0, pY_iterations, po, (i, loopState) =>
                        {
                            float r = 0.0f;
                            float g = 0.0f;
                            float b = 0.0f;
                            if (i % 10 == 0)
                            {
                                r = ovpSettings.majorGridColor.R;
                                g = ovpSettings.majorGridColor.G;
                                b = ovpSettings.majorGridColor.B;
                            }
                            else
                            {
                                r = ovpSettings.minorGridColor.R;
                                g = ovpSettings.minorGridColor.G;
                                b = ovpSettings.minorGridColor.B;
                            }

                            int index = nX_iterations + pX_iterations + nY_iterations + i;
                            index *= 2;
                            gridArray[index] = new Vector3(ovpSettings.cameraPosition.X + (zoom * xLimit), i * spacing, gridZ);
                            gridColorArray[index] = new Vector3(r, g, b);
                            gridArray[index + 1] = new Vector3(ovpSettings.cameraPosition.X + (zoom * -xLimit), i * spacing, gridZ);
                            gridColorArray[index + 1] = new Vector3(r, g, b);
                        });
                    }
                }
            }
            catch (Exception)
            {

            }
        }

        public void drawAxes_VBO()
        {
            if (ovpSettings.showAxes)
            {
                float xLimit = Width * 0.5f;
                float yLimit = Height * 0.5f;
                float zoom = ovpSettings.zoomFactor * ovpSettings.base_zoom;

                axesArray = new Vector3[4];
                axesColorArray = new Vector3[4];
                for (int i = 0; i < axesColorArray.Length; i++)
                {
                    axesColorArray[i] = new Vector3(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B);
                }
                axesArray[0] = new Vector3(0.0f, ovpSettings.cameraPosition.Y + (yLimit * zoom), axisZ);
                axesArray[1] = new Vector3(0.0f, ovpSettings.cameraPosition.Y - (yLimit * zoom), axisZ);
                axesArray[2] = new Vector3(ovpSettings.cameraPosition.X + (xLimit * zoom), 0.0f, axisZ);
                axesArray[3] = new Vector3(ovpSettings.cameraPosition.X - (xLimit * zoom), 0.0f, axisZ);
            }
        }

        public void drawGrid_immediate()
        {
            MakeCurrent();
            GL.LoadIdentity();
            if (ovpSettings.showGrid)
            {
                float spacing = ovpSettings.gridSpacing;
                if (ovpSettings.dynamicGrid)
                {
                    while (WorldToScreen(new SizeF(spacing, 0.0f)).Width > 12.0f)
                        spacing /= 10.0f;

                    while (WorldToScreen(new SizeF(spacing, 0.0f)).Width < 4.0f)
                        spacing *= 10.0f;
                }

                float xLimit = Width * 0.5f;
                float yLimit = Height * 0.5f;
                float zoom = ovpSettings.zoomFactor * ovpSettings.base_zoom;

                float nX_end = -(xLimit * zoom) + ovpSettings.cameraPosition.X;
                int nX_iterations = (int)Math.Abs(Math.Ceiling(nX_end / spacing)) + 1;
                float pX_end = (xLimit * zoom) + ovpSettings.cameraPosition.X;
                int pX_iterations = (int)Math.Abs(Math.Ceiling(pX_end / spacing)) + 1;
                float nY_end = -(yLimit * zoom) + ovpSettings.cameraPosition.Y;
                int nY_iterations = (int)Math.Abs(Math.Ceiling(nY_end / spacing)) + 1;
                float pY_end = (yLimit * zoom) + ovpSettings.cameraPosition.Y;
                int pY_iterations = (int)Math.Abs(Math.Ceiling(pY_end / spacing)) + 1;

                int lines = 2 * (nX_iterations + pX_iterations + nY_iterations + pY_iterations);

                if (WorldToScreen(new SizeF(spacing, 0.0f)).Width >= 4.0f)
                {
                    int k = 0;
                    GL.Begin(PrimitiveType.Lines);

                    for (float i = 0; i < nX_iterations; i++)
                    {
                        if (k <= 1)
                        {
                            GL.Color4(ovpSettings.minorGridColor.R, ovpSettings.minorGridColor.G, ovpSettings.minorGridColor.B, ovpSettings.minorGridColor.A);
                        }
                        if (k == 10)
                        {
                            GL.Color4(ovpSettings.majorGridColor.R, ovpSettings.majorGridColor.G, ovpSettings.majorGridColor.B, ovpSettings.majorGridColor.A);
                            k = 0;
                        }
                        k++;
                        GL.Vertex3(-i * spacing, ovpSettings.cameraPosition.Y + (yLimit * zoom), gridZ);
                        GL.Vertex3(-i * spacing, ovpSettings.cameraPosition.Y - (yLimit * zoom), gridZ);
                    }
                    GL.End();
                    k = 0;
                    GL.Begin(PrimitiveType.Lines);

                    for (float i = 0; i < pX_iterations; i++)
                    {
                        if (k <= 1)
                        {
                            GL.Color4(ovpSettings.minorGridColor.R, ovpSettings.minorGridColor.G, ovpSettings.minorGridColor.B, ovpSettings.minorGridColor.A);
                        }
                        if (k == 10)
                        {
                            GL.Color4(ovpSettings.majorGridColor.R, ovpSettings.majorGridColor.G, ovpSettings.majorGridColor.B, ovpSettings.majorGridColor.A);
                            k = 0;
                        }
                        k++;
                        GL.Vertex3(i * spacing, ovpSettings.cameraPosition.Y + (yLimit * zoom), gridZ);
                        GL.Vertex3(i * spacing, ovpSettings.cameraPosition.Y - (yLimit * zoom), gridZ);
                    }
                    GL.End();
                    k = 0;

                    GL.Begin(PrimitiveType.Lines);
                    for (float i = 0; i < nY_iterations; i++)
                    {
                        if (k <= 1)
                        {
                            GL.Color4(ovpSettings.minorGridColor.R, ovpSettings.minorGridColor.G, ovpSettings.minorGridColor.B, ovpSettings.minorGridColor.A);
                        }
                        if (k == 10)
                        {
                            GL.Color4(ovpSettings.majorGridColor.R, ovpSettings.majorGridColor.G, ovpSettings.majorGridColor.B, ovpSettings.majorGridColor.A);
                            k = 0;
                        }
                        k++;
                        GL.Vertex3(ovpSettings.cameraPosition.X + (xLimit * zoom), -i * spacing, gridZ);
                        GL.Vertex3(ovpSettings.cameraPosition.X - (xLimit * zoom), -i * spacing, gridZ);
                    }
                    GL.End();
                    k = 0;

                    GL.Begin(PrimitiveType.Lines);
                    for (float i = 0; i < pY_iterations; i++)
                    {
                        if (k <= 1)
                        {
                            GL.Color4(ovpSettings.minorGridColor.R, ovpSettings.minorGridColor.G, ovpSettings.minorGridColor.B, ovpSettings.minorGridColor.A);
                        }
                        if (k == 10)
                        {
                            GL.Color4(ovpSettings.majorGridColor.R, ovpSettings.majorGridColor.G, ovpSettings.majorGridColor.B, ovpSettings.majorGridColor.A);
                            k = 0;
                        }
                        k++;
                        GL.Vertex3(ovpSettings.cameraPosition.X + (xLimit * zoom), i * spacing, gridZ);
                        GL.Vertex3(ovpSettings.cameraPosition.X - (xLimit * zoom), i * spacing, gridZ);
                    }
                    GL.End();
                }
            }
        }

        public void drawAxes_immediate()
        {
            MakeCurrent();
            GL.LoadIdentity();
            if (ovpSettings.showAxes)
            {
                float xLimit = Width * 0.5f;
                float yLimit = Height * 0.5f;
                float zoom = ovpSettings.zoomFactor * ovpSettings.base_zoom;

                GL.Color4(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, ovpSettings.axisColor.A);
                GL.Begin(PrimitiveType.Lines);
                GL.Vertex3(0.0f, ovpSettings.cameraPosition.Y + (yLimit * zoom), axisZ);
                GL.Vertex3(0.0f, ovpSettings.cameraPosition.Y - (yLimit * zoom), axisZ);
                GL.Vertex3(ovpSettings.cameraPosition.X + (xLimit * zoom), 0.0f, axisZ);
                GL.Vertex3(ovpSettings.cameraPosition.X - (xLimit * zoom), 0.0f, axisZ);
                GL.End();
            }
        }
    }
}

