using Eto.Drawing;
using System.Numerics;
using Veldrid;

namespace VeldridEto;

public partial class VeldridDriver
{
	private bool drawing = false;
	private bool done_drawing = false;
	private async void pUpdateViewport()
	{
		if ((!ovpSettings.changed) || (Surface.GraphicsDevice == null) ||
		    (!Surface.Visible) || (Surface.Width <= 0) || (Surface.Height <= 0) || drawing)
		{
			return;
		}

		drawing = true;
		done_drawing = false;
		
		// Trying to push things into tasks to speed up the computation. Not sure if this is entirely robust.
		await Task.WhenAll(drawAxes(), drawGrid(), drawLines(), drawPolygons());

		bool b1 = await drawAxes();
		bool b2 = await drawGrid();
		bool b3 = await drawLines();
		bool b4 = await drawPolygons();

		done_drawing = b1 && b2 && b3 && b4;
	}
	
	private async Task<bool> drawPolygons()
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
			polyIndices = Array.Empty<uint>();
			pointsIndices = Array.Empty<uint>();
			polyVertexCount = new uint[numPolys];

			tessFirst = new uint[tessPolyListCount];
			tessIndices = Array.Empty<uint>();
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
				int tess_polyIndex = 0;
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
						tess_polyIndex++;
					}

					tessVertexCount[poly] = 3;
				}
				
				tessIndices = new uint[tess_polyIndex];
				for (int i = 0; i < tess_polyIndex; i++)
				{
					tessIndices[i] = (uint)i;
				}
			}

			// Pondering options here - this would make a nice border construct around the filled geometry, amongst other things.
			int poly_index = 0;
			int line_polyIndex = 0;
			int line_pointsIndex = 0;
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
					line_polyIndex++;

					polyList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt + 1].X,
							ovpSettings.polyList[poly].poly[pt + 1].Y, polyZ),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					counter++;
					line_polyIndex++;

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
					line_pointsIndex++;

					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					line_pointsIndex++;
					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					line_pointsIndex++;

					tFirst.Add(tCounter);
					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					line_pointsIndex++;
					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					line_pointsIndex++;
					pointsList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f,
							ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha)));
					tCounter++;
					line_pointsIndex++;
				}

				polyVertexCount[poly] = (uint)(counter - previouscounter); // set our vertex count for the polygon.
				poly_index ++;
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
					line_polyIndex++;
					polyList.Add(new VertexPositionColor(
						new Vector3(ovpSettings.bgPolyList[poly].poly[pt + 1].X,
							ovpSettings.bgPolyList[poly].poly[pt + 1].Y, polyZ),
						new RgbaFloat(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G,
							ovpSettings.bgPolyList[poly].color.B, alpha)));
					counter++;
					line_polyIndex++;
				}

				polyVertexCount[poly] = (uint)(counter - previouscounter); // set our vertex count for the polygon.
				poly_index++;
			}

			polyIndices = new uint[line_polyIndex];
			for (int i = 0; i < line_polyIndex; i++)
			{
				polyIndices[i] = (uint)i;
			}
			pointsIndices = new uint[line_pointsIndex];
			for (int i = 0; i < line_pointsIndex; i++)
			{
				pointsIndices[i] = (uint)i;
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
			updateBuffer(ref PolysIndexBuffer, polyIndices, sizeof(uint), BufferUsage.IndexBuffer);
		}

		if (ovpSettings.drawPoints() && polyListCount > 0)
		{
			updateBuffer(ref PointsVertexBuffer, pointsList.ToArray(), VertexPositionColor.SizeInBytes,
				BufferUsage.VertexBuffer);
			updateBuffer(ref PointsIndexBuffer, pointsIndices, sizeof(uint), BufferUsage.IndexBuffer);
		}

		if (ovpSettings.drawFilled() && tessPolyListCount > 0)
		{
			updateBuffer(ref TessVertexBuffer, tessPolyList.ToArray(), VertexPositionColor.SizeInBytes,
				BufferUsage.VertexBuffer);
			updateBuffer(ref TessIndexBuffer, tessIndices, sizeof(uint), BufferUsage.IndexBuffer);
		}

		return true;
	}

	private async Task<bool>  drawLines()
	{
		int linesCount = ovpSettings.lineList.Count;

		switch (linesCount)
		{
			// Create our first and count arrays for the vertex indices, to enable polygon separation when rendering.
			case > 0:
			{
				int point_count = 0;
				List<VertexPositionColor> lineList = new();
				linesIndices = new uint[linesCount];

				// Carve our Z-space up to stack polygons
				float polyZStep = 1.0f / ovpSettings.lineList.Count;

				lineFirst = new uint[linesCount];
				lineVertexCount = new uint[linesCount];

				int lineIndex = 0;
				for (int poly = 0; poly < linesCount; poly++)
				{
					float alpha = ovpSettings.lineList[poly].alpha;
					float polyZ = poly * polyZStep;
					lineFirst[poly] = (uint)lineList.Count;
					lineList.AddRange(ovpSettings.lineList[poly].poly.Select(t =>
						new VertexPositionColor(new Vector3(t.X, t.Y, polyZ),
							new RgbaFloat(ovpSettings.lineList[poly].color.R, ovpSettings.lineList[poly].color.G,
								ovpSettings.lineList[poly].color.B, alpha))));
					point_count += ovpSettings.lineList[poly].poly.Length;
					lineIndex += lineList.Count;
					lineVertexCount[poly] =
						(uint)ovpSettings.lineList[poly].poly.Length; // set our vertex count for the polygon.
				}

				linesIndices = new uint[point_count];
				for (int i = 0; i < point_count; i++)
				{
					linesIndices[i] = (uint)i;
				}

				updateBuffer(ref LinesVertexBuffer, lineList.ToArray(), VertexPositionColor.SizeInBytes,
					BufferUsage.VertexBuffer);
				updateBuffer(ref LinesIndexBuffer, linesIndices, sizeof(uint), BufferUsage.IndexBuffer);
				
				break;
			}
			default:
				LinesVertexBuffer = null;
				break;
		}

		return true;
	}

	private async Task<bool>  drawGrid()
	{
		if (!ovpSettings.drawGrid())
		{
			return true;
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

		return true;
	}

	private async Task<bool>  drawAxes()
	{
		if (!ovpSettings.drawAxes())
		{
			return true;
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

		return true;
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

		if (gridIndices.Length != 0)
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

		if (axesIndices.Length != 0)
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
			if (tessIndices.Length != 0)
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
						CommandList.SetIndexBuffer(TessIndexBuffer, IndexFormat.UInt32);
						CommandList.SetPipeline(FilledPipeline);
						CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
						CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

						CommandList.DrawIndexed(
							indexCount:(uint)tessIndices.Length,
							instanceCount:1,
							indexStart:0,
							vertexOffset:0,
							instanceStart:0);
						/*
						for (int l = 0; l < tessVertexCount.Length; l++)
						{
							CommandList.Draw(tessVertexCount[l], 1, tessFirst[l], 0);
						}
						*/
					}
					catch (Exception ex)
					{
						Console.WriteLine("Ex: " + ex);
					}
				}
			}
		}

		if (polyIndices.Length != 0)
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
					CommandList.SetIndexBuffer(PolysIndexBuffer, IndexFormat.UInt32);
					CommandList.SetPipeline(LinesPipeline);
					CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

					CommandList.DrawIndexed(
						indexCount:(uint)polyIndices.Length,
						instanceCount:1,
						indexStart:0,
						vertexOffset:0,
						instanceStart:0);
					/*
					for (int l = 0; l < polyVertexCount.Length; l++)
					{
						CommandList.Draw(polyVertexCount[l], 1, polyFirst[l], 0);
					}
					*/
				}
				catch (Exception ex)
				{
					Console.WriteLine("Ex: " + ex);
				}
			}
		}

		if (linesIndices.Length != 0 && ovpSettings.drawDrawn())
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
					CommandList.SetIndexBuffer(LinesIndexBuffer, IndexFormat.UInt32);
					CommandList.SetPipeline(LinesPipeline);
					CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

					CommandList.DrawIndexed(
						indexCount:(uint)linesIndices.Length,
						instanceCount:1,
						indexStart:0,
						vertexOffset:0,
						instanceStart:0);

					/*
					for (int l = 0; l < lineVertexCount.Length; l++)
					{
						CommandList.Draw(lineVertexCount[l], 1, lineFirst[l], 0);
					}
					*/
				}
				catch (Exception ex)
				{
					Console.WriteLine("Ex: " + ex);
				}
			}
		}

		if (ovpSettings.drawPoints())
		{
			if (polyIndices.Length != 0)
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
						CommandList.SetIndexBuffer(PointsIndexBuffer, IndexFormat.UInt32);
						CommandList.SetPipeline(FilledPipeline);
						CommandList.SetGraphicsResourceSet(0, ViewMatrixSet);
						CommandList.SetGraphicsResourceSet(1, ModelMatrixSet);

						CommandList.DrawIndexed(
							indexCount:(uint)pointsIndices.Length,
							instanceCount:1,
							indexStart:0,
							vertexOffset:0,
							instanceStart:0);
						/*
						foreach (uint t in pointsFirst)
						{
							CommandList.Draw(3, 1, t, 0);
						}
						*/
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
