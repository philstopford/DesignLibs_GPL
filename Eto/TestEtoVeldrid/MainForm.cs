using Eto.Drawing;
using Eto.Forms;
using Eto.Veldrid;
using System;
using System.Globalization;
using System.Threading;
using Veldrid;
using VeldridEto;

namespace TestEtoVeldrid;

public partial class MainForm : Form
{
	private VeldridSurface Surface;

	private OVPSettings ovpSettings;

	private VeldridDriver Driver;

	private bool _veldridReady;
	public bool VeldridReady
	{
		get { return _veldridReady; }
		private set
		{
			_veldridReady = value;
			SetUpVeldrid();
		}
	}

	private bool _formReady;
	public bool FormReady
	{
		get { return _formReady; }
		set
		{
			_formReady = value;
			SetUpVeldrid();
		}
	}

	public MainForm() : this(VeldridSurface.PreferredBackend)
	{
	}
	public MainForm(GraphicsBackend backend)
	{
		InitializeComponent();

		Shown += (sender, e) => FormReady = true;

		GraphicsDeviceOptions options = new(
			false,
			Veldrid.PixelFormat.R32Float,
			false,
			ResourceBindingModel.Improved);

		Surface = new VeldridSurface(backend, options);
		Surface.VeldridInitialized += (sender, e) => VeldridReady = true;

		Content = Surface;

		ovpSettings = new OVPSettings();
		ovpSettings.drawFilled(true);
		ovpSettings.drawDrawn(true);
		ovpSettings.drawPoints(true);

		addPolys();

		Driver = new VeldridDriver(ref ovpSettings, ref Surface)
		{
			Surface = Surface
		};

		// Setup progressive loading event handlers
		Driver.ProgressiveLoadingProgress += OnProgressiveLoadingProgress;
		Driver.ProgressiveLoadingStatusChanged += OnProgressiveLoadingStatusChanged;

		Surface.VeldridInitialized += (sender, e) =>
		{
			Driver.SetUpVeldrid();
			
			_veldridReady = true;
		};
		
		// TODO: Make this binding actually work both ways.
		// CmdAnimate.Bind<bool>("Checked", Driver, "Animate");
		// CmdClockwise.Bind<bool>("Checked", Driver, "Clockwise");
	}

	private void SetUpVeldrid()
	{
		switch (FormReady && VeldridReady)
		{
			case false:
				return;
		}

		Title = $"Veldrid backend: {Surface.Backend.ToString()}";

		createVPContextMenu();

		Driver.Clock.Start();
	}

	private void addPolys()
	{
		addPolys(false);
	}

	private void addPolys(bool addManyPolygons)
	{
		ovpSettings.clear();

		float r = 0.0f;

		PointF[] testPoly = new PointF[6];
		testPoly[0] = new PointF(2.0f + r, 2.0f + r);
		testPoly[1] = new PointF(15.0f + r, 12.0f + r);
		testPoly[2] = new PointF(8.0f + r, 24.0f + r);
		testPoly[3] = new PointF(8.0f + r, 15.0f + r);
		testPoly[4] = new PointF(3.0f + r, 2.0f + r);
		testPoly[5] = new PointF(2.0f + r, 2.0f + r);


		PointF[] testPoly2 = new PointF[6];
		testPoly2[0] = new PointF(12.0f + r, 2.0f + r);
		testPoly2[1] = new PointF(25.0f + r, 12.0f + r);
		testPoly2[2] = new PointF(18.0f + r, 24.0f + r);
		testPoly2[3] = new PointF(18.0f + r, 15.0f + r);
		testPoly2[4] = new PointF(13.0f + r, 2.0f + r);
		testPoly2[5] = new PointF(12.0f + r, 2.0f + r);

		ovpSettings.addPolygon(testPoly2, Color.FromArgb(0, 255, 255), 1.0f, true, 1);

		ovpSettings.addPolygon(testPoly, Color.FromArgb(255, 0, 0), 1.0f, true, 2);

		r = -30.0f;

		PointF[] testPolyF1 = new PointF[6];
		testPolyF1[0] = new PointF(2.0f + r, 2.0f + r);
		testPolyF1[1] = new PointF(15.0f + r, 12.0f + r);
		testPolyF1[2] = new PointF(8.0f + r, 24.0f + r);
		testPolyF1[3] = new PointF(8.0f + r, 15.0f + r);
		testPolyF1[4] = new PointF(3.0f + r, 2.0f + r);
		testPolyF1[5] = new PointF(2.0f + r, 2.0f + r);

		PointF[] testPolyF2 = new PointF[6];
		testPolyF2[0] = new PointF(12.0f + r, 2.0f + r);
		testPolyF2[1] = new PointF(25.0f + r, 12.0f + r);
		testPolyF2[2] = new PointF(18.0f + r, 24.0f + r);
		testPolyF2[3] = new PointF(18.0f + r, 15.0f + r);
		testPolyF2[4] = new PointF(13.0f + r, 2.0f + r);
		testPolyF2[5] = new PointF(12.0f + r, 2.0f + r);


		ovpSettings.addPolygon(testPolyF2, Color.FromArgb(0, 255, 255), 0.5f, false, 3);

		ovpSettings.addPolygon(testPolyF1, Color.FromArgb(255, 0, 0), 1.0f, false, 4);

		r = 30.0f;

		PointF[] testPolyBG1 = new PointF[6];
		testPolyBG1[0] = new PointF(2.0f + r, 2.0f + r);
		testPolyBG1[1] = new PointF(15.0f + r, 12.0f + r);
		testPolyBG1[2] = new PointF(8.0f + r, 24.0f + r);
		testPolyBG1[3] = new PointF(8.0f + r, 15.0f + r);
		testPolyBG1[4] = new PointF(3.0f + r, 2.0f + r);
		testPolyBG1[5] = new PointF(2.0f + r, 2.0f + r);

		PointF[] testPolyBG2 = new PointF[6];
		testPolyBG2[0] = new PointF(12.0f + r, 2.0f + r);
		testPolyBG2[1] = new PointF(25.0f + r, 12.0f + r);
		testPolyBG2[2] = new PointF(18.0f + r, 24.0f + r);
		testPolyBG2[3] = new PointF(18.0f + r, 15.0f + r);
		testPolyBG2[4] = new PointF(13.0f + r, 2.0f + r);
		testPolyBG2[5] = new PointF(12.0f + r, 2.0f + r);


		ovpSettings.addBGPolygon(testPolyBG2, Color.FromArgb(0, 255, 255), 1.0f, 3);

		ovpSettings.addBGPolygon(testPolyBG1, Color.FromArgb(255, 0, 0), 1.0f, 4);

		// Add performance test case with many polygons
		if (addManyPolygons)
		{
			addManyTestPolygons();
		}
	}

	private void addManyTestPolygons()
	{
		Console.WriteLine("Adding many polygons for performance testing...");
		var stopwatch = System.Diagnostics.Stopwatch.StartNew();
		
		const int polyCount = 5000; // Start with 5000 polygons for performance testing
		Random rand = new(42); // Fixed seed for reproducible results
		
		for (int i = 0; i < polyCount; i++)
		{
			// Create random polygon positions across a larger area
			float baseX = (float)(rand.NextDouble() * 2000.0 - 1000.0);
			float baseY = (float)(rand.NextDouble() * 2000.0 - 1000.0);
			float size = (float)(rand.NextDouble() * 20.0 + 5.0);
			
			// Create a simple rectangular polygon
			PointF[] poly = new PointF[5];
			poly[0] = new PointF(baseX, baseY);
			poly[1] = new PointF(baseX + size, baseY);
			poly[2] = new PointF(baseX + size, baseY + size);
			poly[3] = new PointF(baseX, baseY + size);
			poly[4] = new PointF(baseX, baseY); // Close the polygon
			
			// Use random colors
			Color color = Color.FromArgb(
				rand.Next(128, 256),
				rand.Next(128, 256),
				rand.Next(128, 256));
			
			float alpha = (float)(rand.NextDouble() * 0.5 + 0.5); // Alpha between 0.5 and 1.0
			
			ovpSettings.addPolygon(poly, color, alpha, false, i + 100);
		}
		
		stopwatch.Stop();
		Console.WriteLine($"Added {polyCount} polygons in {stopwatch.ElapsedMilliseconds}ms");
	}

	private void OnProgressiveLoadingProgress(VeldridEto.BatchLoadProgress progress)
	{
		// Update UI with loading progress
		Application.Instance.Invoke(() =>
		{
			Title = $"Veldrid backend: {Surface.Backend} - Loading: {progress.ProcessedCount}/{progress.TotalCount} ({progress.OverallProgress:P1})";
		});
	}

	private void OnProgressiveLoadingStatusChanged(string status)
	{
		// Update UI with loading status
		Application.Instance.Invoke(() =>
		{
			Console.WriteLine($"Status: {status}");
		});
	}

	private ContextMenu vp_menu;

	private void createVPContextMenu()
	{
		// Single viewport now mandates regeneration of the context menu each time, to allow for entry screening.
		vp_menu = new ContextMenu();

		int itemIndex = 0;
		vp_menu.Items.Add(new ButtonMenuItem { Text = "Reset" });
		vp_menu.Items[itemIndex].Click += delegate
		{
			Driver.reset();
			updateViewport();
		};
		itemIndex++;

		ButtonMenuItem VPMenuDisplayOptionsMenu = vp_menu.Items.GetSubmenu("Display Options");
		itemIndex++;
		int displayOptionsSubItemIndex = 0;
		VPMenuDisplayOptionsMenu.Items.Add(new ButtonMenuItem { Text = "Toggle AA" });
		VPMenuDisplayOptionsMenu.Items[displayOptionsSubItemIndex].Click += delegate
		{
			ovpSettings.aA(!ovpSettings.aA());
			updateViewport();
		};
		displayOptionsSubItemIndex++;
		VPMenuDisplayOptionsMenu.Items.Add(new ButtonMenuItem { Text = "Toggle Fill" });
		VPMenuDisplayOptionsMenu.Items[displayOptionsSubItemIndex].Click += delegate
		{
			ovpSettings.drawFilled(!ovpSettings.drawFilled());
			updateViewport();
		};
		displayOptionsSubItemIndex++;
		VPMenuDisplayOptionsMenu.Items.Add(new ButtonMenuItem { Text = "Toggle Points" });
		VPMenuDisplayOptionsMenu.Items[displayOptionsSubItemIndex].Click += delegate
		{
			ovpSettings.drawPoints(!ovpSettings.drawPoints());
			updateViewport();
		};
		displayOptionsSubItemIndex++;
		VPMenuDisplayOptionsMenu.Items.Add(new ButtonMenuItem { Text = "Test Performance" });
		VPMenuDisplayOptionsMenu.Items[displayOptionsSubItemIndex].Click += delegate
		{
			addPolys(true);
			updateViewport();
		};
		displayOptionsSubItemIndex++;
		VPMenuDisplayOptionsMenu.Items.Add(new ButtonMenuItem { Text = "Toggle Progressive Loading" });
		VPMenuDisplayOptionsMenu.Items[displayOptionsSubItemIndex].Click += delegate
		{
			ovpSettings.setProgressiveLoading(!ovpSettings.progressiveLoadingEnabled());
			updateViewport();
		};
		displayOptionsSubItemIndex++;

		{
			if (Driver.ovpSettings.isLocked())
			{
				vp_menu.Items.Add(new ButtonMenuItem { Text = "Thaw" });
			}
			else
			{
				vp_menu.Items.Add(new ButtonMenuItem { Text = "Freeze" });
			}
			vp_menu.Items[itemIndex].Click += delegate
			{
				Driver.freeze_thaw();
				updateViewport();
			};
			itemIndex++;
			vp_menu.Items.AddSeparator();
			itemIndex++;
			vp_menu.Items.Add(new ButtonMenuItem { Text = "Save bookmark" });
			vp_menu.Items[itemIndex].Click += delegate
			{
				Driver.saveLocation();
			};
			itemIndex++;
			vp_menu.Items.Add(new ButtonMenuItem { Text = "Load bookmark" });
			vp_menu.Items[itemIndex].Click += delegate
			{
				Driver.loadLocation();
			};
			vp_menu.Items[itemIndex].Enabled = Driver.savedLocation_valid switch
			{
				false => false,
				_ => vp_menu.Items[itemIndex].Enabled
			};
			itemIndex++;
		}
		vp_menu.Items.AddSeparator();
		itemIndex++;
		vp_menu.Items.Add(new ButtonMenuItem { Text = "Zoom Extents" });
		vp_menu.Items[itemIndex].Click += delegate
		{
			Driver.zoomExtents(-1);
		};
		itemIndex++;
		vp_menu.Items.AddSeparator();
		itemIndex++;
		vp_menu.Items.Add(new ButtonMenuItem { Text = "Zoom In" });
		vp_menu.Items[itemIndex].Click += delegate
		{
			Driver.zoomIn(-1);
			updateViewport();
		};
		itemIndex++;

		ButtonMenuItem VPMenuZoomInMenu = vp_menu.Items.GetSubmenu("Fast Zoom In");
		itemIndex++;
		int zoomInSubItemIndex = 0;
		VPMenuZoomInMenu.Items.Add(new ButtonMenuItem { Text = "Zoom In (x5)" });
		VPMenuZoomInMenu.Items[zoomInSubItemIndex].Click += delegate
		{
			Driver.zoomIn(-50);
			updateViewport();
		};
		zoomInSubItemIndex++;
		VPMenuZoomInMenu.Items.Add(new ButtonMenuItem { Text = "Zoom In (x10)" });
		VPMenuZoomInMenu.Items[zoomInSubItemIndex].Click += delegate
		{
			Driver.zoomIn(-100);
			updateViewport();
		};
		zoomInSubItemIndex++;
		VPMenuZoomInMenu.Items.Add(new ButtonMenuItem { Text = "Zoom In (x50)" });
		VPMenuZoomInMenu.Items[zoomInSubItemIndex].Click += delegate
		{
			Driver.zoomIn(-500);
			updateViewport();
		};
		zoomInSubItemIndex++;
		VPMenuZoomInMenu.Items.Add(new ButtonMenuItem { Text = "Zoom In (x100)" });
		VPMenuZoomInMenu.Items[zoomInSubItemIndex].Click += delegate
		{
			Driver.zoomIn(-1000);
			updateViewport();
		};
		zoomInSubItemIndex++;

		vp_menu.Items.AddSeparator();
		itemIndex++;

		vp_menu.Items.Add(new ButtonMenuItem { Text = "Zoom Out" });
		vp_menu.Items[itemIndex].Click += delegate
		{
			Driver.zoomOut(-1);
			updateViewport();
		};
		itemIndex++;

		ButtonMenuItem VPMenuZoomOutMenu = vp_menu.Items.GetSubmenu("Fast Zoom Out");
		itemIndex++;
		int zoomOutSubItemIndex = 0;
		VPMenuZoomOutMenu.Items.Add(new ButtonMenuItem { Text = "Zoom Out (x5)" });
		VPMenuZoomOutMenu.Items[zoomOutSubItemIndex].Click += delegate
		{
			Driver.zoomOut(-50);
			updateViewport();
		};
		zoomOutSubItemIndex++;
		VPMenuZoomOutMenu.Items.Add(new ButtonMenuItem { Text = "Zoom Out (x10)" });
		VPMenuZoomOutMenu.Items[zoomOutSubItemIndex].Click += delegate
		{
			Driver.zoomOut(-100);
			updateViewport();
		};
		zoomOutSubItemIndex++;
		VPMenuZoomOutMenu.Items.Add(new ButtonMenuItem { Text = "Zoom Out (x50)" });
		VPMenuZoomOutMenu.Items[zoomOutSubItemIndex].Click += delegate
		{
			Driver.zoomOut(-500);
			updateViewport();
		};
		zoomOutSubItemIndex++;
		VPMenuZoomOutMenu.Items.Add(new ButtonMenuItem { Text = "Zoom Out (x100)" });
		VPMenuZoomOutMenu.Items[zoomOutSubItemIndex].Click += delegate
		{
			Driver.zoomOut(-1000);
			updateViewport();
		};
		zoomOutSubItemIndex++;

		Driver.setContextMenu(ref vp_menu);
	}

	private void updateViewport()
	{
		Application.Instance.Invoke(() =>
		{
			Monitor.Enter(ovpSettings);
			// Note: Removed automatic re-adding of polygons to allow for performance testing
			// addPolys();

			try
			{
				createVPContextMenu();
				Driver.updateViewport();
				// viewPort.Invalidate();
			}
			catch (Exception)
			{
			}
			finally
			{
				Monitor.Exit(ovpSettings);
			}
		});
	}
}