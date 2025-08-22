using Eto.Drawing;
using System.Numerics;
using System.Runtime.CompilerServices;
using Veldrid;

namespace VeldridEto;

public partial class VeldridDriver
{
	private bool drawing = false;
	private bool done_drawing = false;
	private async void pUpdateViewport()
	{
		if ((!ovpSettings.changed) || (Surface!.GraphicsDevice == null) ||
			(!Surface.Visible) || (Surface.Width <= 0) || (Surface.Height <= 0) || drawing)
		{
			return;
		}

		drawing = true;
		done_drawing = false;
		
		// Trying to push things into tasks to speed up the computation. Not sure if this is entirely robust.
		await Task.WhenAll(drawAxes(), drawGrid(), drawLines(), drawPolygons());
		
		done_drawing = true;
	}

	private int fgPolyListCount;
	private int bgPolyListCount;
	private int tessPolyListCount;

	private VertexPositionColor[]? tessPolyList;
	private VertexPositionColor[]? polyList;
	private VertexPositionColor[]? pointsList;

	private Task drawPolygons()
	{
		fgPolyListCount = ovpSettings.polyList!.Count;
		bgPolyListCount = ovpSettings.bgPolyList!.Count;
		tessPolyListCount = ovpSettings.tessPolyList!.Count;

		int[] pointCountBeforeCurrentPolygon_fg;
		int totalPointCount_fg;
		if (fgPolyListCount > 0)
		{
			pointCountBeforeCurrentPolygon_fg = new int[fgPolyListCount];
			pointCountBeforeCurrentPolygon_fg[0] = 0;
			totalPointCount_fg = ovpSettings.polyList[0].poly.Length;
		}
		else
		{
			pointCountBeforeCurrentPolygon_fg = new int[0];
			totalPointCount_fg = 0;
		}

		for (int i = 1; i < fgPolyListCount; i++)
		{
			int previousPolygonPointCount = ovpSettings.polyList[i - 1].poly.Length;
			pointCountBeforeCurrentPolygon_fg[i] = pointCountBeforeCurrentPolygon_fg[i-1] + previousPolygonPointCount;
			totalPointCount_fg += ovpSettings.polyList[i].poly.Length;
		}

		int[] pointCountBeforeCurrentPolygon_bg;
		int totalPointCount_bg;
		if (bgPolyListCount > 0)
		{
			pointCountBeforeCurrentPolygon_bg = new int[bgPolyListCount];
			pointCountBeforeCurrentPolygon_bg[0] = 0;
			totalPointCount_bg = ovpSettings.bgPolyList[0].poly.Length;
		}
		else
		{
			pointCountBeforeCurrentPolygon_bg = new int[0];
			totalPointCount_bg = 0;
		}

		for (int i = 1; i < bgPolyListCount; i++)
		{
			int previousPolygonPointCount = ovpSettings.bgPolyList[i - 1].poly.Length;
			pointCountBeforeCurrentPolygon_bg[i] = pointCountBeforeCurrentPolygon_bg[i-1] + previousPolygonPointCount;
			totalPointCount_bg += ovpSettings.bgPolyList[i].poly.Length;
		}

		// Start and end points for each polygon are not duplicated.
		int totalPolyListCount = ((totalPointCount_fg) + (totalPointCount_bg)) * 2;
		polyList = new VertexPositionColor[totalPolyListCount];
		// Two triangles per point in foreground polygon.
		int pointListCount = (totalPointCount_fg) * 6;
		pointsList = new VertexPositionColor[pointListCount];

		try
		{
			// Carve our Z-space up to stack polygons
			int numPolys = totalPolyListCount;
			// Create our first and count arrays for the vertex indices, to enable polygon separation when rendering.
			polyIndices = new uint[0];
			pointsIndices = new uint[0];
			tessIndices = new uint[0];
			
			if (ovpSettings.drawFilled())
			{
				numPolys += tessPolyListCount;
			}

			float
				polyZStep = 1.0f / Math.Max(1,
				numPolys +
				1); // avoid a div by zero risk; pad the poly number also to reduce risk of adding a poly beyond the clipping range
			
			if (ovpSettings.drawFilled())
			{
				tessPolyList = new VertexPositionColor[tessPolyListCount * 3];
				for (int poly = 0; poly < tessPolyListCount; poly++)
				{
					float alpha = ovpSettings.tessPolyList[poly].alpha;
					float polyZ = poly * polyZStep;
					for (int pt = 0; pt < 3; pt++)
					{
						tessPolyList[(poly * 3) + pt] = new VertexPositionColor(
							new Vector3(ovpSettings.tessPolyList[poly].poly[pt].X,
								ovpSettings.tessPolyList[poly].poly[pt].Y, polyZ),
							new RgbaFloat(ovpSettings.tessPolyList[poly].color.R,
								ovpSettings.tessPolyList[poly].color.G, ovpSettings.tessPolyList[poly].color.B,
								alpha));
					}
				}
				
				tessIndices = new uint[tessPolyListCount * 3];
				for (int i = 0; i < tessPolyListCount * 3; i++)
				{
					tessIndices[i] = (uint)i;
				}
			}

			// Pondering options here - this would make a nice border construct around the filled geometry, amongst other things.
			for (int poly = 0; poly < fgPolyListCount; poly++)
			{
				float alpha = ovpSettings.polyList[poly].alpha;
				if (ovpSettings.drawFilled())
				{
					alpha = 1.0f;
				}

				float polyZ = (tessPolyListCount * polyZStep) + (poly * polyZStep);
				int polyLength = ovpSettings.polyList[poly].poly.Length - 1;
				int poly_index_offset = 0;
				int point_index_offset = 0;
				if (poly > 0)
				{
					// Start and end points for each polygon are not duplicated.
					poly_index_offset = ((pointCountBeforeCurrentPolygon_fg[poly] * 2) - 2);
					point_index_offset = ((pointCountBeforeCurrentPolygon_fg[poly] * 6) - 6);
				}

				for (int pt = 0; pt < polyLength; pt++)
				{
					polyList[poly_index_offset + (pt * 2)] = new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt].X, ovpSettings.polyList[poly].poly[pt].Y,
							polyZ),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha));

					polyList[poly_index_offset + ((pt * 2) + 1)] = new VertexPositionColor(
						new Vector3(ovpSettings.polyList[poly].poly[pt + 1].X,
							ovpSettings.polyList[poly].poly[pt + 1].Y, polyZ),
						new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
							ovpSettings.polyList[poly].color.B, alpha));

					if (ovpSettings.drawPoints())
					{
						try
						{
							pointsList[point_index_offset + (pt * 6)] = new VertexPositionColor(
								new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f,
								ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f),
								new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
								ovpSettings.polyList[poly].color.B, alpha));
							pointsList[point_index_offset + ((pt * 6) + 1)] = new VertexPositionColor(
								new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f,
								ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f),
								new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
								ovpSettings.polyList[poly].color.B, alpha));
							pointsList[point_index_offset + ((pt * 6) + 2)] = new VertexPositionColor(
								new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f,
								ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f),
								new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
								ovpSettings.polyList[poly].color.B, alpha));
							pointsList[point_index_offset + ((pt * 6) + 3)] = new VertexPositionColor(
								new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f,
								ovpSettings.polyList[poly].poly[pt].Y - pointWidth / 2.0f, 1.0f),
								new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
								ovpSettings.polyList[poly].color.B, alpha));
							pointsList[point_index_offset + ((pt * 6) + 4)] = new VertexPositionColor(
								new Vector3(ovpSettings.polyList[poly].poly[pt].X - pointWidth / 2.0f,
								ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f),
								new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
								ovpSettings.polyList[poly].color.B, alpha));
							pointsList[point_index_offset + ((pt * 6) + 5)] = new VertexPositionColor(
								new Vector3(ovpSettings.polyList[poly].poly[pt].X + pointWidth / 2.0f,
								ovpSettings.polyList[poly].poly[pt].Y + pointWidth / 2.0f, 1.0f),
								new RgbaFloat(ovpSettings.polyList[poly].color.R, ovpSettings.polyList[poly].color.G,
								ovpSettings.polyList[poly].color.B, alpha));
						}
						catch (Exception ex)
						{
							Console.WriteLine(ex);
						}
					}
				}

			for (int poly = 0; poly < bgPolyListCount; poly++)
			{
				float alpha = ovpSettings.bgPolyList[poly].alpha;
				float polyZ = poly * polyZStep;

				int fg_poly_index_offset = (totalPointCount_fg - 2) * 2;
				if (fg_poly_index_offset < 0)
				{
					fg_poly_index_offset = 0;
				}

				int poly_index_offset = 0;
				if (poly > 0)
				{
					// Start and end points for each polygon are not duplicated.
					poly_index_offset = ((pointCountBeforeCurrentPolygon_bg[poly] * 2) - 2);
				}

				int bgPolyLength = ovpSettings.bgPolyList[poly].poly.Length - 1;
				for (int pt = 0; pt < bgPolyLength; pt++)
				{
					try
					{
						polyList[fg_poly_index_offset + poly_index_offset + (pt * 2)] = new VertexPositionColor(
							new Vector3(ovpSettings.bgPolyList[poly].poly[pt].X,
							ovpSettings.bgPolyList[poly].poly[pt].Y, polyZ),
							new RgbaFloat(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G,
							ovpSettings.bgPolyList[poly].color.B, alpha));

						polyList[fg_poly_index_offset + poly_index_offset + ((pt * 2) + 1)] = new VertexPositionColor(
							new Vector3(ovpSettings.bgPolyList[poly].poly[pt + 1].X,
							ovpSettings.bgPolyList[poly].poly[pt + 1].Y, polyZ),
							new RgbaFloat(ovpSettings.bgPolyList[poly].color.R, ovpSettings.bgPolyList[poly].color.G,
							ovpSettings.bgPolyList[poly].color.B, alpha));
					}
					catch (Exception e)
					{
						Console.WriteLine(e);
					}
				}
			}

			polyIndices = new uint[totalPolyListCount];
			for (int i = 0; i < totalPolyListCount; i++) // (int i = 0; i < totalPolyListCount; i++)
			{
				polyIndices[i] = (uint)i;
			}
		}
		catch (Exception ex)
		{
			Console.WriteLine(ex);
			
			// Can ignore - not critical.
		}
		
		return Task.CompletedTask;
	}

	private void updatePolygonBuffers()
	{
		if (fgPolyListCount > 0 || bgPolyListCount > 0)
		{
			updateBuffer(ref PolysVertexBuffer, ref PolysVertexBufferCapacity, polyList!, VertexPositionColor.SizeInBytes,
				BufferUsage.VertexBuffer);
			updateBuffer(ref PolysIndexBuffer, ref PolysIndexBufferCapacity, polyIndices!, sizeof(uint), BufferUsage.IndexBuffer);
		}

		if (ovpSettings.drawPoints() && fgPolyListCount > 0)
		{
			updateBuffer(ref PointsVertexBuffer, ref PointsVertexBufferCapacity, pointsList!, VertexPositionColor.SizeInBytes,
				BufferUsage.VertexBuffer);
			updateBuffer(ref PointsIndexBuffer, ref PointsIndexBufferCapacity, pointsIndices!, sizeof(uint), BufferUsage.IndexBuffer);
		}

		if (!ovpSettings.drawFilled() || tessPolyListCount <= 0)
		{
			return;
		}
		updateBuffer(ref TessVertexBuffer, ref TessVertexBufferCapacity, tessPolyList!, VertexPositionColor.SizeInBytes,
			BufferUsage.VertexBuffer);
		updateBuffer(ref TessIndexBuffer, ref TessIndexBufferCapacity, tessIndices!, sizeof(uint), BufferUsage.IndexBuffer);
	}

	private int linesCount;

	private VertexPositionColor[]? lineList;

	private Task drawLines()
	{
		linesCount = ovpSettings.lineList!.Count;
		if (linesCount == 0)
		{
			lineList = new VertexPositionColor[0];
			linesIndices = new uint[0];
			return Task.CompletedTask;
		}

		try
		{
			// Experimental stuff... can we do this with parallel loops?
			int[] pointCountBeforeCurrentPolygon = new int[linesCount];
			pointCountBeforeCurrentPolygon[0] = 0;
			int totalPointCount = ovpSettings.lineList[0].poly.Length;
			for (int i = 1; i < linesCount; i++)
			{
				int previousPolygonPointCount = ovpSettings.lineList[i - 1].poly.Length;
				pointCountBeforeCurrentPolygon[i] = pointCountBeforeCurrentPolygon[i-1] + previousPolygonPointCount;
				totalPointCount += ovpSettings.lineList[i].poly.Length;
			}

			// Start and end points for each polygon are not duplicated.
			// Allow for line shape to be closed.
			lineList = new VertexPositionColor[(totalPointCount) * 2];
			for (int poly = 0; poly < linesCount; poly++)
			{
				float alpha = ovpSettings.lineList[poly].alpha;
				float polyZ = 1.0f;

				int polyLength = ovpSettings.lineList[poly].poly.Length;
				int index_offset = 0;
				if (poly > 0)
				{
					// Start and end points for each polygon are not duplicated.
					index_offset = (pointCountBeforeCurrentPolygon[poly] * 2);
				}
				for (int pt = 0; pt < polyLength - 1; pt++)
				{
					lineList[index_offset + (pt * 2)] = (new VertexPositionColor(
						new Vector3(ovpSettings.lineList[poly].poly[pt].X,
							ovpSettings.lineList[poly].poly[pt].Y,
							polyZ),
						new RgbaFloat(ovpSettings.lineList[poly].color.R, ovpSettings.lineList[poly].color.G,
							ovpSettings.lineList[poly].color.B, alpha)));

					lineList[index_offset + ((pt * 2) + 1)] = (new VertexPositionColor(
						new Vector3(ovpSettings.lineList[poly].poly[pt + 1].X,
							ovpSettings.lineList[poly].poly[pt + 1].Y, polyZ),
						new RgbaFloat(ovpSettings.lineList[poly].color.R, ovpSettings.lineList[poly].color.G,
							ovpSettings.lineList[poly].color.B, alpha)));
				}
			}
			
			int counter = lineList.Length;
			
			linesIndices = new uint[counter];
			for (int i = 0; i < counter; i++) // (int i = 0; i < counter; i++)
			{
				linesIndices[i] = (uint)i;
			}
		}
		catch (Exception ex)
		{
			Console.WriteLine(ex);
			// Can ignore - not critical.
		}
		
		return Task.CompletedTask;
	}

	private void updateLineBuffers()
	{
		if (linesCount <= 0)
		{
			return;
		}
		updateBuffer(ref LinesVertexBuffer, ref LinesVertexBufferCapacity, lineList!, VertexPositionColor.SizeInBytes,
			BufferUsage.VertexBuffer);
		updateBuffer(ref LinesIndexBuffer, ref LinesIndexBufferCapacity, linesIndices!, sizeof(uint), BufferUsage.IndexBuffer);
	}

	private List<VertexPositionColor>? grid;
	private Task  drawGrid()
	{
		if (!ovpSettings.drawGrid())
		{
			return Task.CompletedTask;
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

		grid = new List<VertexPositionColor>();

		if (WorldToScreen(new SizeF(spacing, 0.0f)).Width >= 4.0f)
		{
			int k = 0;
			for (float i = 0; i > -(Surface!.RenderWidth * zoom) + x; i -= spacing)
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
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * Surface!.RenderHeight, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * -Surface!.RenderHeight, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
			}

			k = 0;
			for (float i = 0; i < Surface!.RenderWidth * zoom + x; i += spacing)
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
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * Surface!.RenderHeight, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(i, y + zoom * -Surface!.RenderHeight, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
			}

			k = 0;
			for (float i = 0; i > -(Surface!.RenderHeight * zoom) + y; i -= spacing)
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
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * Surface!.RenderWidth, i, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * -Surface!.RenderWidth, i, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
			}

			k = 0;
			for (float i = 0; i < Surface!.RenderHeight * zoom + y; i += spacing)
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
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * Surface!.RenderWidth, i, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
				grid.Add(new VertexPositionColor(new Vector3(x + zoom * -Surface!.RenderWidth, i, gridZ),
					new RgbaFloat(r, g, b, 1.0f)));
			}
		}
		
		return Task.CompletedTask;
	}

	private void updateGridBuffers()
	{
		if (grid == null)
		{
			return;
		}
		uint gridCount = (uint)grid!.Count;

		switch (gridCount)
		{
			case > 0:
			{
				gridIndices = new uint[gridCount];
				for (int i = 0; i < gridIndices.Length; i++) //; i < gridIndices.Length; i++)
				{
					gridIndices[i] = (uint)i;
				}
				updateBuffer(ref GridVertexBuffer, ref GridVertexBufferCapacity, grid.ToArray(), VertexPositionColor.SizeInBytes,
					BufferUsage.VertexBuffer);
				updateBuffer(ref GridIndexBuffer, ref GridIndexBufferCapacity, gridIndices, sizeof(uint), BufferUsage.IndexBuffer);

				break;
			}
			default:
				GridVertexBuffer = null;
				GridIndexBuffer = null;
				break;
		}
	}

	private VertexPositionColor[]? axesArray;
	private Task  drawAxes()
	{
		if (!ovpSettings.drawAxes())
		{
			return Task.CompletedTask;
		}

		float zoom = ovpSettings.getBaseZoom() * ovpSettings.getZoomFactor();
		axesArray = new VertexPositionColor[4];
		axesArray[0] =
			new VertexPositionColor(
				new Vector3(0.0f, ovpSettings.getCameraY() + Surface!.RenderHeight * zoom, axisZ),
				new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));
		axesArray[1] =
			new VertexPositionColor(
				new Vector3(0.0f, ovpSettings.getCameraY() - Surface!.RenderHeight * zoom, axisZ),
				new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));
		axesArray[2] =
			new VertexPositionColor(new Vector3(ovpSettings.getCameraX() + Surface!.RenderWidth * zoom, 0.0f, axisZ),
				new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));
		axesArray[3] =
			new VertexPositionColor(new Vector3(ovpSettings.getCameraX() - Surface!.RenderWidth * zoom, 0.0f, axisZ),
				new RgbaFloat(ovpSettings.axisColor.R, ovpSettings.axisColor.G, ovpSettings.axisColor.B, 1.0f));

		axesIndices = new uint[4] {0,1,2,3};
		
		return Task.CompletedTask;
	}

	private void updateAxesBuffers()
	{
		updateBuffer(ref AxesVertexBuffer, ref AxesVertexBufferCapacity, axesArray!, VertexPositionColor.SizeInBytes, BufferUsage.VertexBuffer);
		updateBuffer(ref AxesIndexBuffer, ref AxesIndexBufferCapacity, axesIndices!, sizeof(uint), BufferUsage.IndexBuffer);
	}
	
	public void Draw()
	{
		if (!Ready)
		{
			return;
		}

		CommandList!.Begin();

		ModelMatrix *= Matrix4x4.CreateFromAxisAngle(
			new Vector3(0, 0, 1), 0);
		CommandList.UpdateBuffer(ModelBuffer, 0, ModelMatrix);

		float zoom = ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom();

		float left = ovpSettings.getCameraX() - (float)Surface!.RenderWidth / 2 * zoom;
		float right = ovpSettings.getCameraX() + (float)Surface!.RenderWidth / 2 * zoom;
		float bottom = ovpSettings.getCameraY() + (float)Surface!.RenderHeight / 2 * zoom;
		float top = ovpSettings.getCameraY() - (float)Surface!. RenderHeight / 2 * zoom;

		ViewMatrix = Matrix4x4.CreateOrthographicOffCenter(left, right, bottom, top, 0.0f, 1.0f);
		CommandList.UpdateBuffer(ViewBuffer, 0, ViewMatrix);

		CommandList.SetFramebuffer(Surface!.Swapchain!.Framebuffer);

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
				Surface!.GraphicsDevice!.SubmitCommands(CommandList!);
			}

			Surface!.GraphicsDevice.SwapBuffers(Surface!.Swapchain);
		}
		catch (Exception ex)
		{
			Console.WriteLine("Ex: " + ex);
		}

	}

	private void populate_command_list()
	{
		updateGridBuffers();
		updateAxesBuffers();
		updatePolygonBuffers();
		updateLineBuffers();

		if (GridVertexBuffer != null && gridIndices != null && gridIndices.Length != 0)
		{
			if (LinePipeline == null)
			{
				throw new Exception("populate_command_list : LinePipeline not initialized!");
			}

			lock (GridVertexBuffer)
			{
				try
				{
					CommandList!.SetVertexBuffer(0, GridVertexBuffer);
					CommandList!.SetIndexBuffer(GridIndexBuffer, IndexFormat.UInt32);
					CommandList!.SetPipeline(LinePipeline);
					CommandList!.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList!.SetGraphicsResourceSet(1, ModelMatrixSet);

					CommandList!.DrawIndexed(
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

		if (AxesVertexBuffer != null && axesIndices != null && axesIndices.Length != 0)
		{
			if (LinePipeline == null)
			{
				throw new Exception("populate_command_list : LinePipeline not initialized!");
			}

			lock (AxesVertexBuffer)
			{
				try
				{
					CommandList!.SetVertexBuffer(0, AxesVertexBuffer);
					CommandList!.SetIndexBuffer(AxesIndexBuffer, IndexFormat.UInt32);
					CommandList!.SetPipeline(LinePipeline);
					CommandList!.SetGraphicsResourceSet(0, ViewMatrixSet);
					CommandList!.SetGraphicsResourceSet(1, ModelMatrixSet);

					CommandList!.DrawIndexed(
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
			if (TessVertexBuffer != null && tessIndices != null && tessIndices.Length != 0)
			{
				if (LinePipeline == null)
				{
					throw new Exception("populate_command_list : FilledPipeline not initialized!");
				}

				lock (TessVertexBuffer)
				{
					try
					{
						CommandList!.SetVertexBuffer(0, TessVertexBuffer);
						CommandList!.SetIndexBuffer(TessIndexBuffer, IndexFormat.UInt32);
						CommandList!.SetPipeline(FilledPipeline);
						CommandList!.SetGraphicsResourceSet(0, ViewMatrixSet);
						CommandList!.SetGraphicsResourceSet(1, ModelMatrixSet);

						CommandList!.DrawIndexed(
							indexCount:(uint)tessIndices.Length,
							instanceCount:1,
							indexStart:0,
							vertexOffset:0,
							instanceStart:0);
					}
					catch (Exception ex)
					{
						Console.WriteLine("Ex: " + ex);
					}
				}
			}

			if (PolysVertexBuffer != null && polyIndices != null && polyIndices.Length != 0)
			{
				if (LinePipeline == null)
				{
					throw new Exception("populate_command_list : LinesPipeline not initialized!");
				}

				lock (PolysVertexBuffer)
				{
					try
					{
						CommandList!.SetVertexBuffer(0, PolysVertexBuffer);
						CommandList!.SetIndexBuffer(PolysIndexBuffer, IndexFormat.UInt32);
						CommandList!.SetPipeline(LinesPipeline);
						CommandList!.SetGraphicsResourceSet(0, ViewMatrixSet);
						CommandList!.SetGraphicsResourceSet(1, ModelMatrixSet);

						CommandList!.DrawIndexed(
							indexCount:(uint)polyIndices.Length,
							instanceCount:1,
							indexStart:0,
							vertexOffset:0,
							instanceStart:0);
					}
					catch (Exception ex)
					{
						Console.WriteLine("Ex: " + ex);
					}
				}
			}

			if (LinesVertexBuffer != null && linesIndices != null && linesIndices.Length != 0 && ovpSettings.drawDrawn())
			{
				if (LinePipeline == null)
				{
					throw new Exception("populate_command_list : LinesPipeline not initialized!");
				}

				lock (LinesVertexBuffer)
				{
					try
					{
						CommandList!.SetVertexBuffer(0, LinesVertexBuffer);
						CommandList!.SetIndexBuffer(LinesIndexBuffer, IndexFormat.UInt32);
						CommandList!.SetPipeline(LinesPipeline);
						CommandList!.SetGraphicsResourceSet(0, ViewMatrixSet);
						CommandList!.SetGraphicsResourceSet(1, ModelMatrixSet);

						CommandList!.DrawIndexed(
							indexCount:(uint)linesIndices.Length,
							instanceCount:1,
							indexStart:0,
							vertexOffset:0,
							instanceStart:0);
					}
					catch (Exception ex)
					{
						Console.WriteLine("Ex: " + ex);
					}
				}
			}

			if (ovpSettings.drawPoints())
			{
				if (PointsVertexBuffer != null && polyIndices != null && polyIndices.Length != 0)
				{
					if (LinePipeline == null)
					{
						throw new Exception("populate_command_list : FilledPipeline not initialized!");
					}

					lock (PointsVertexBuffer)
					{
						try
						{
							CommandList!.SetVertexBuffer(0, PointsVertexBuffer);
							CommandList!.SetIndexBuffer(PointsIndexBuffer, IndexFormat.UInt32);
							CommandList!.SetPipeline(FilledPipeline);
							CommandList!.SetGraphicsResourceSet(0, ViewMatrixSet);
							CommandList!.SetGraphicsResourceSet(1, ModelMatrixSet);

							CommandList!.DrawIndexed(
								indexCount:(uint)pointsIndices!.Length,
								instanceCount:1,
								indexStart:0,
								vertexOffset:0,
								instanceStart:0);
						}
						catch (Exception ex)
						{
							Console.WriteLine("Ex: " + ex);
						}
					}
				}
			}
		}
	}

	private void updateBuffer<T>(ref DeviceBuffer? buffer, ref uint bufferCapacity, T[]? data, uint elementSize, BufferUsage usage)
		where T : unmanaged
	{
		if (data == null)
		{
			return;
		}
		switch (data.Length)
		{
			case > 0:
			{
				// Reuse existing buffer unless capacity insufficient
				uint neededSize = elementSize * (uint)data.Length;
				ResourceFactory factory = Surface!.GraphicsDevice!.ResourceFactory;

				if (buffer == null || bufferCapacity < neededSize)
				{
					buffer?.Dispose();
					buffer = factory.CreateBuffer(new BufferDescription(neededSize, usage));
					bufferCapacity = neededSize;
				}

				Surface.GraphicsDevice.UpdateBuffer(buffer, 0, data);
				break;
			}
		}
	}

	// Capacity fields for vertex/index buffers to avoid reallocation each frame
	private uint PolysVertexBufferCapacity = 0;
	private uint PolysIndexBufferCapacity = 0;
	private uint PointsVertexBufferCapacity = 0;
	private uint PointsIndexBufferCapacity = 0;
	private uint TessVertexBufferCapacity = 0;
	private uint TessIndexBufferCapacity = 0;
	private uint LinesVertexBufferCapacity = 0;
	private uint LinesIndexBufferCapacity = 0;
	private uint GridVertexBufferCapacity = 0;
	private uint GridIndexBufferCapacity = 0;
	private uint AxesVertexBufferCapacity = 0;
	private uint AxesIndexBufferCapacity = 0;
}