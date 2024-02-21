using Eto.Drawing;
using System.Numerics;
using Veldrid;

namespace VeldridEto;

public partial class VeldridDriver
{
	private void pUpdateViewport()
	{
		if ((!ovpSettings.changed) || (Surface.GraphicsDevice == null) ||
		    (!Surface.Visible) || (Surface.Width <= 0) || (Surface.Height <= 0))
		{
			return;
		}

		drawAxes();
		drawGrid();
		drawLines();
		drawPolygons();
		updateHostFunc?.Invoke();
		Surface.Invalidate();
		ovpSettings.changed = false;
	}

	private void drawPolygons()
	{
		int polyListCount = ovpSettings.polyList.Count;
		int bgPolyListCount = ovpSettings.bgPolyList.Count;
		int tessPolyListCount = ovpSettings.tessPolyList.Count;

		List<VertexPositionColor> polyList = new();

		List<VertexPositionColor> pointsList = new();

		List<VertexPositionColor> tessPolyList = new();

		try
		{

			// Carve our Z-space up to stack polygons
			int numPolys = 1;

			numPolys = polyListCount + bgPolyListCount;
			// Create our first and count arrays for the vertex indices, to enable polygon separation when rendering.
			polyFirst = new uint[numPolys];
			polyVertexCount = new uint[numPolys];

			tessFirst = new uint[tessPolyListCount];
			tessVertexCount = new uint[tessPolyListCount];

			List<uint> tFirst = new();

			uint tCounter = 0;

			if (ovpSettings.drawFilled())
			{
				numPolys += tessPolyListCount;
			}

			float
				polyZStep = 1.0f / Math.Max(1,
					numPolys +
					1); // avoid a div by zero risk; pad the poly number also to reduce risk of adding a poly beyond the clipping range

			int counter = 0; // vertex count that will be used to define 'first' index for each polygon.
			int previouscounter = 0; // will be used to derive the number of vertices in each polygon.

			float polyZ = 0;

			if (ovpSettings.drawFilled())
			{
				for (int poly = 0; poly < tessPolyListCount; poly++)
				{
					tessFirst[poly] = (uint)(poly * 3);
					float alpha = ovpSettings.tessPolyList[poly].alpha;
					polyZ += polyZStep;
					for (int pt = 0; pt < 3; pt++)
					{
						tessPolyList.Add(new VertexPositionColor(
							new Vector3(ovpSettings.tessPolyList[poly].poly[pt].X,
								ovpSettings.tessPolyList[poly].poly[pt].Y, polyZ),
							new RgbaFloat(ovpSettings.tessPolyList[poly].color.R,
								ovpSettings.tessPolyList[poly].color.G, ovpSettings.tessPolyList[poly].color.B,
								alpha)));
					}

					tessVertexCount[poly] = 3;
				}
			}

			// Pondering options here - this would make a nice border construct around the filled geometry, amongst other things.
			for (int poly = 0; poly < polyListCount; poly++)
			{
				float alpha = ovpSettings.polyList[poly].alpha;
				if (ovpSettings.drawFilled())
				{
					alpha = 1.0f;
				}

				polyZ += polyZStep;
				polyFirst[poly] = (uint)counter;
				previouscounter = counter;
				int polyLength = ovpSettings.polyList[poly].poly.Length - 1;
				for (int pt = 0; pt < polyLength; pt++)
				{
					polyList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X, ovpSettings.polyList[poly].poly[pt].Y,
							polyZ),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					counter++;
					polyList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt + 1].X,
							ovpSettings.polyList[poly].poly[pt + 1].Y, polyZ),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					counter++;

					if (!ovpSettings.drawPoints())
					{
						continue;
					}

					tFirst.Add(tCounter);
					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;

					tFirst.Add(tCounter);
					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
				}

				polyVertexCount[poly] = (uint)(counter - previouscounter); // set our vertex count for the polygon.
			}

			polyZ = 0;
			for (int poly = 0; poly < bgPolyListCount; poly++)
			{
				float alpha = ovpSettings.bgPolyList[poly].alpha;
				polyZ += polyZStep;
				polyFirst[poly + polyListCount] = (uint)counter;
				previouscounter = counter;

				int bgPolyLength = ovpSettings.bgPolyList[poly].poly.Length - 1;
				for (int pt = 0; pt < bgPolyLength; pt++)
				{
					polyList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.bgPolyList[poly].poly[pt].X,
							ovpSettings.bgPolyList[poly].poly[pt].Y, polyZ),
						new RgbaFloat(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G,
							ovpSettings.bgPolyList[poly].color.B, alpha)));
					counter++;
					polyList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.bgPolyList[poly].poly[pt + 1].X,
							ovpSettings.bgPolyList[poly].poly[pt + 1].Y, polyZ),
						new RgbaFloat(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G,
							ovpSettings.bgPolyList[poly].color.B, alpha)));
					counter++;
				}

				polyVertexCount[poly + polyListCount] =
					(uint)(counter - previouscounter); // set our vertex count for the polygon.
			}

			pointsFirst = tFirst.ToArray();
		}
		catch (Exception)
		{
			// Can ignore - not critical.
		}

		if (polyListCount > 0 || bgPolyListCount > 0)
		{
			updateBuffer(ref PolysVertexBuffer, polyList.ToArray(), VertexPositionColor.SizeInBytes,
				BufferUsage.VertexBuffer);
		}

		if (ovpSettings.drawPoints() && polyListCount > 0)
		{
			updateBuffer(ref PointsVertexBuffer, pointsList.ToArray(), VertexPositionColor.SizeInBytes,
				BufferUsage.VertexBuffer);
		}

		if (ovpSettings.drawFilled() && tessPolyListCount > 0)
		{
			updateBuffer(ref TessVertexBuffer, tessPolyList.ToArray(), VertexPositionColor.SizeInBytes,
				BufferUsage.VertexBuffer);
		}
	}

	private void drawLines()
	{
		int tmp = ovpSettings.lineList.Count;

		switch (tmp)
		{
			// Create our first and count arrays for the vertex indices, to enable polygon separation when rendering.
			case > 0:
			{
				List<VertexPositionColor> lineList = new();

				// Carve our Z-space up to stack polygons
				float polyZStep = 1.0f / ovpSettings.lineList.Count;

				lineFirst = new uint[tmp];
				lineVertexCount = new uint[tmp];

				for (int poly = 0; poly < tmp; poly++)
				{
					float alpha = ovpSettings.lineList[poly].alpha;
					float polyZ = poly * polyZStep;
					lineFirst[poly] = (uint)lineList.Count;
					lineList.AddRange(ovpSettings.lineList[poly].poly.Select(t =>
						new VertexPositionColor(new Vector3(t.X, t.Y, polyZ),
							new RgbaFloat(ovpSettings.lineList[poly].color.R, ovpSettings.lineList[poly].color.G,
								ovpSettings.lineList[poly].color.B, alpha))));
					lineVertexCount[poly] =
						(uint)ovpSettings.lineList[poly].poly.Length; // set our vertex count for the polygon.
				}

				updateBuffer(ref LinesVertexBuffer, lineList.ToArray(), VertexPositionColor.SizeInBytes,
					BufferUsage.VertexBuffer);
				break;
			}
			default:
				LinesVertexBuffer = null;
				break;
		}
	}

	private void drawGrid()
	{
		if (!ovpSettings.drawGrid())
		{
			return;
		}

		float spacing = ovpSettings.gridSpacing();
		if (ovpSettings.isGridDynamic())
		{
			while (WorldToScreen(new SizeF(spacing, 0.0f)).Width > 12.0f)
			{
				spacing /= 10.0f;
			}

			while (WorldToScreen(new SizeF(spacing, 0.0f)).Width < 4.0f)
			{
				spacing *= 10.0f;
			}
		}

		float zoom = ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom();
		float x = ovpSettings.getCameraX();
		float y = ovpSettings.getCameraY();

		List<VertexPositionColor> grid = new();

		if (WorldToScreen(new SizeF(spacing, 0.0f)).Width >= 4.0f)
		{
			int k = 0;
			for (float i = 0; i > -(Surface.RenderWidth * zoom) + x; i -= spacing)
			{
				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;
				switch (k)
				{
					case <= 9:
						r = ovpSettings.minorGridColor.R;
						g = ovpSettings.minorGridColor.G;
						b = ovpSettings.minorGridColor.B;
						break;
					case 10:
						r = ovpSettings.majorGridColor.R;
						g = ovpSettings.majorGridColor.G;
						b = ovpSettings.majorGridColor.B;
						k = 0;
						break;
				}

				k++;
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * Surface.RenderHeight, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * -Surface.RenderHeight, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
			}

			k = 0;
			for (float i = 0; i < Surface.RenderWidth * zoom + x; i += spacing)
			{
				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;
				switch (k)
				{
					case <= 9:
						r = ovpSettings.minorGridColor.R;
						g = ovpSettings.minorGridColor.G;
						b = ovpSettings.minorGridColor.B;
						break;
					case 10:
						r = ovpSettings.majorGridColor.R;
						g = ovpSettings.majorGridColor.G;
						b = ovpSettings.majorGridColor.B;
						k = 0;
						break;
				}

				k++;
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * Surface.RenderHeight, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * -Surface.RenderHeight, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
			}

			k = 0;
			for (float i = 0; i > -(Surface.RenderHeight * zoom) + y; i -= spacing)
			{
				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;
				switch (k)
				{
					case <= 9:
						r = ovpSettings.minorGridColor.R;
						g = ovpSettings.minorGridColor.G;
						b = ovpSettings.minorGridColor.B;
						break;
					case 10:
						r = ovpSettings.majorGridColor.R;
						g = ovpSettings.majorGridColor.G;
						b = ovpSettings.majorGridColor.B;
						k = 0;
						break;
				}

				k++;
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * Surface.RenderWidth, i, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * -Surface.RenderWidth, i, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
			}

			k = 0;
			for (float i = 0; i < Surface.RenderHeight * zoom + y; i += spacing)
			{
				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;
				switch (k)
				{
					case <= 9:
						r = ovpSettings.minorGridColor.R;
						g = ovpSettings.minorGridColor.G;
						b = ovpSettings.minorGridColor.B;
						break;
					case 10:
						r = ovpSettings.majorGridColor.R;
						g = ovpSettings.majorGridColor.G;
						b = ovpSettings.majorGridColor.B;
						k = 0;
						break;
				}

				k++;
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * Surface.RenderWidth, i, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * -Surface.RenderWidth, i, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
			}
		}

		uint gridCount = (uint)grid.Count;

		switch (gridCount)
		{
			case > 0:
			{
				gridIndices = new uint[gridCount];
				for (uint i = 0; i < gridIndices.Length; i++)
				{
					gridIndices[i] = i;
				}

				updateBuffer(ref GridVertexBuffer, grid.ToArray(), VertexPositionColor.SizeInBytes,
					BufferUsage.VertexBuffer);
				updateBuffer(ref GridIndexBuffer, gridIndices, sizeof(uint), BufferUsage.IndexBuffer);
				break;
			}
			default:
				GridVertexBuffer = null;
				GridIndexBuffer = null;
				break;
		}
	}

	private void drawAxes()
	{
		if (!ovpSettings.drawAxes())
		{
			return;
		}

		float zoom = ovpSettings.getBaseZoom() * ovpSettings.getZoomFactor();
		VertexPositionColor[] axesArray = new VertexPositionColor[4];
		axesArray[0] =
			new VertexPositionColor(
				new Vector3(0.0f, ovpSettings.getCameraY() + Surface.RenderHeight * zoom, axisZ),
				new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));
		axesArray[1] =
			new VertexPositionColor(
				new Vector3(0.0f, ovpSettings.getCameraY() - Surface.RenderHeight * zoom, axisZ),
				new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));
		axesArray[2] =
			new VertexPositionColor(new Vector3(ovpSettings.getCameraX() + Surface.RenderWidth * zoom, 0.0f, axisZ),
				new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));
		axesArray[3] =
			new VertexPositionColor(new Vector3(ovpSettings.getCameraX() - Surface.RenderWidth * zoom, 0.0f, axisZ),
				new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));

		axesIndices = new uint[4] { 0, 1, 2, 3 };

		updateBuffer(ref AxesVertexBuffer, axesArray, VertexPositionColor.SizeInBytes, BufferUsage.VertexBuffer);
		updateBuffer(ref AxesIndexBuffer, axesIndices, sizeof(uint), BufferUsage.IndexBuffer);
	}

	public void Draw()
	{
		if (!Ready)
		{
			return;
		}

		CommandList.Begin();

		ModelMatrix *= Matrix4x4.CreateFromAxisAngle(
			new Vector3(0, 0, 1), 0);
		CommandList.UpdateBuffer(ModelBuffer, 0, ModelMatrix);

		float zoom = ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom();

		float left = ovpSettings.getCameraX() - (float)Surface.RenderWidth / 2 * zoom;
		float right = ovpSettings.getCameraX() + (float)Surface.RenderWidth / 2 * zoom;
		float bottom = ovpSettings.getCameraY() + (float)Surface.RenderHeight / 2 * zoom;
		float top = ovpSettings.getCameraY() - (float)Surface.RenderHeight / 2 * zoom;

		ViewMatrix = Matrix4x4.CreateOrthographicOffCenter(left, right, bottom, top, 0.0f, 1.0f);
		CommandList.UpdateBuffer(ViewBuffer, 0, ViewMatrix);

		CommandList.SetFramebuffer(Surface.Swapchain.Framebuffer);

		// These commands differ from the stock Veldrid "Getting Started"
		// tutorial in two ways. First, the viewport is cleared to pink
		// instead of black so as to more easily distinguish between errors
		// in creating a graphics context and errors drawing vertices within
		// said context. Second, this project creates its swapchain with a
		// depth buffer, and that buffer needs to be reset at the start of
		// each frame.

		RgbaFloat bgColor = new(ovpSettings.backColor.R, ovpSettings.backColor.G, ovpSettings.backColor.B, 1.0f);

		CommandList.ClearColorTarget(0, bgColor);
		CommandList.ClearDepthStencil(1.0f);

		populate_command_list();

		CommandList.End();

		try
		{
			lock (CommandList)
			{
				Surface.GraphicsDevice.SubmitCommands(CommandList);
			}

			Surface.GraphicsDevice.SwapBuffers(Surface.Swapchain);
		}
		catch (Exception ex)
		{
			Console.WriteLine("Ex: " + ex);
		}

	}

	private void populate_command_list()
	{
		drawGrid();
		drawAxes();
		drawLines();
		drawPolygons();

		if (GridVertexBuffer != null)
		{
			if (LinePipeline == null)
			{
				throw new("populate_command_list : LinePipeline not initialized!");
			}

			lock (GridVertexBuffer)
			{
				try
				{
					CommandList.SetVertexBuffer(0, GridVertexBuffer);
					CommandList.SetIndexBuffer(GridIndexBuffer, IndexFormat.UInt32);
					CommandList.SetPipeline(LinePipeline);
					CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

					CommandList.DrawIndexed(
						indexCount: (uint)gridIndices.Length,
						instanceCount: 1,
						indexStart: 0,
						vertexOffset: 0,
						instanceStart: 0);
				}
				catch (Exception ex)
				{
					Console.WriteLine("Ex: " + ex);
				}
			}
		}

		if (AxesVertexBuffer != null)
		{
			if (LinePipeline == null)
			{
				throw new("populate_command_list : LinePipeline not initialized!");
			}

			lock (AxesVertexBuffer)
			{
				try
				{
					CommandList.SetVertexBuffer(0, AxesVertexBuffer);
					CommandList.SetIndexBuffer(AxesIndexBuffer, IndexFormat.UInt32);
					CommandList.SetPipeline(LinePipeline);
					CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

					CommandList.DrawIndexed(
						indexCount: (uint)axesIndices.Length,
						instanceCount: 1,
						indexStart: 0,
						vertexOffset: 0,
						instanceStart: 0);
				}
				catch (Exception ex)
				{
					Console.WriteLine("Ex: " + ex);
				}
			}
		}

		if (ovpSettings.drawFilled())
		{
			if (TessVertexBuffer != null)
			{
				if (LinePipeline == null)
				{
					throw new("populate_command_list : FilledPipeline not initialized!");
				}

				lock (TessVertexBuffer)
				{
					try
					{
						CommandList.SetVertexBuffer(0, TessVertexBuffer);
						CommandList.SetPipeline(FilledPipeline);
						CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
						CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

						for (int l = 0; l < tessVertexCount.Length; l++)
						{
							CommandList.Draw(tessVertexCount[l], 1, tessFirst[l], 0);
						}
					}
					catch (Exception ex)
					{
						Console.WriteLine("Ex: " + ex);
					}
				}
			}
		}

		if (PolysVertexBuffer != null)
		{
			if (LinePipeline == null)
			{
				throw new("populate_command_list : LinesPipeline not initialized!");
			}

			lock (PolysVertexBuffer)
			{
				try
				{
					CommandList.SetVertexBuffer(0, PolysVertexBuffer);
					CommandList.SetPipeline(LinesPipeline);
					CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

					for (int l = 0; l < polyVertexCount.Length; l++)
					{
						CommandList.Draw(polyVertexCount[l], 1, polyFirst[l], 0);
					}
				}
				catch (Exception ex)
				{
					Console.WriteLine("Ex: " + ex);
				}
			}
		}

		if (LinesVertexBuffer != null && ovpSettings.drawDrawn())
		{
			if (LinePipeline == null)
			{
				throw new("populate_command_list : LinesPipeline not initialized!");
			}

			lock (LinesVertexBuffer)
			{
				try
				{
					CommandList.SetVertexBuffer(0, LinesVertexBuffer);
					CommandList.SetPipeline(LinesPipeline);
					CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

					for (int l = 0; l < lineVertexCount.Length; l++)
					{
						CommandList.Draw(lineVertexCount[l], 1, lineFirst[l], 0);
					}
				}
				catch (Exception ex)
				{
					Console.WriteLine("Ex: " + ex);
				}
			}
		}

		if (ovpSettings.drawPoints())
		{
			if (PointsVertexBuffer != null)
			{
				if (LinePipeline == null)
				{
					throw new("populate_command_list : FilledPipeline not initialized!");
				}

				lock (PointsVertexBuffer)
				{
					try
					{
						CommandList.SetVertexBuffer(0, PointsVertexBuffer);
						CommandList.SetPipeline(FilledPipeline);
						CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
						CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

						foreach (uint t in pointsFirst)
						{
							CommandList.Draw(3, 1, t, 0);
						}
					}
					catch (Exception ex)
					{
						Console.WriteLine("Ex: " + ex);
					}
				}
			}
		}
	}

	private void updateBuffer<T>(ref DeviceBuffer buffer, T[] data, uint elementSize, BufferUsage usage)
		where T : unmanaged
	{
		switch (data.Length)
		{
			case > 0:
			{
				buffer?.Dispose();

				ResourceFactory factory = Surface.GraphicsDevice.ResourceFactory;

				buffer = factory.CreateBuffer(new BufferDescription(elementSize * (uint)data.Length, usage));

				Surface.GraphicsDevice.UpdateBuffer(buffer, 0, data);
				break;
			}
		}
	}
}
