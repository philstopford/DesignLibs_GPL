using System;
using Eto.Drawing;
using Eto.Mac.Forms;
using Eto.Mac.Forms.Controls;
using Eto.Veldrid;
using Eto.Veldrid.Mac64;
using MonoMac.AppKit;
using MonoMac.CoreGraphics;
using Veldrid;

[assembly: Eto.ExportHandler(typeof(VeldridSurface), typeof(Mac64VeldridSurfaceHandler))]

namespace Eto.Veldrid.Mac64
{
	public class MacVeldridView : NSView, IMacControl
	{
		public WeakReference? WeakHandler { get; set; }

		public event EventHandler? DrawRectCalled;
		public event EventHandler? FrameSizeChanged;

		public MacVeldridView()
		{
			WantsLayer = true;
			PostsFrameChangedNotifications = true;
		}

		public MacVeldridView(IntPtr handle) : base(handle)
		{
			WantsLayer = true;
			PostsFrameChangedNotifications = true;
		}

		public override bool AcceptsFirstResponder() => true;

		public override void SetFrameSize(CGSize newSize)
		{
			base.SetFrameSize(newSize);
			FrameSizeChanged?.Invoke(this, EventArgs.Empty);
		}

		public override void DrawRect(CGRect dirtyRect)
		{
			DrawRectCalled?.Invoke(this, EventArgs.Empty);
		}
	}

	public class Mac64VeldridSurfaceHandler : MacView<MacVeldridView, VeldridSurface, VeldridSurface.ICallback>, VeldridSurface.IHandler
	{
		public override NSView ContainerControl => Control;

		public Size RenderSize
		{
			get
			{
				var size = Widget.Size;
				var scale = Scale;
				if (size.Width > 0 && size.Height > 0)
				{
					return Size.Round((SizeF)size * scale);
				}

				var frame = Control.Frame;
				if (frame.Width > 0 && frame.Height > 0)
				{
					return new Size((int)Math.Round(frame.Width * scale), (int)Math.Round(frame.Height * scale));
				}

				return new Size(Math.Max(1, size.Width), Math.Max(1, size.Height));
			}
		}

		private float Scale => Widget.ParentWindow?.LogicalPixelSize ?? Widget.ParentWindow?.Screen?.LogicalPixelSize ?? (float)(Control.Window?.BackingScaleFactor ?? 1.0f);

		public Mac64VeldridSurfaceHandler()
		{
		}

		protected override MacVeldridView CreateControl()
		{
			var view = new MacVeldridView();
			view.FrameSizeChanged += (sender, e) =>
			{
				if (Widget != null && Widget.Loaded)
				{
					Callback.OnResize(Widget, new ResizeEventArgs(RenderSize));
				}
			};
			return view;
		}

		public Swapchain? CreateSwapchain()
		{
			Swapchain? swapchain;

			if (Widget.Backend == GraphicsBackend.OpenGL)
			{
				swapchain = Widget.GraphicsDevice?.MainSwapchain;
			}
			else
			{
				var source = SwapchainSource.CreateNSView(Control.Handle);
				var renderSize = RenderSize;
				var width = (uint)Math.Max(1, renderSize.Width);
				var height = (uint)Math.Max(1, renderSize.Height);
				swapchain = Widget.GraphicsDevice?.ResourceFactory.CreateSwapchain(
					new SwapchainDescription(
						source,
						width,
						height,
						Widget.GraphicsDeviceOptions.SwapchainDepthFormat,
						Widget.GraphicsDeviceOptions.SyncToVerticalBlank,
						Widget.GraphicsDeviceOptions.SwapchainSrgbFormat));
				swapchain?.Resize(width, height);
			}

			return swapchain;
		}

		public override void OnLoad(EventArgs e)
		{
			base.OnLoad(e);
			Callback.OnInitializeBackend(Widget, new InitializeEventArgs(RenderSize));
			Widget.SizeChanged += (sender, args) =>
			{
				if (Widget.Loaded)
				{
					Callback.OnResize(Widget, new ResizeEventArgs(RenderSize));
				}
			};
		}

		public override void OnSizeChanged(EventArgs e)
		{
			base.OnSizeChanged(e);
			if (Widget.Loaded)
			{
				Callback.OnResize(Widget, new ResizeEventArgs(RenderSize));
			}
		}

		public override void Invalidate(bool invalidateChildren)
		{
			base.Invalidate(invalidateChildren);
			if (Widget != null && Widget.Loaded)
			{
				Callback.OnDraw(Widget, EventArgs.Empty);
			}
		}

		public override void Invalidate(Rectangle rect, bool invalidateChildren)
		{
			base.Invalidate(rect, invalidateChildren);
			if (Widget != null && Widget.Loaded)
			{
				Callback.OnDraw(Widget, EventArgs.Empty);
			}
		}

		public override void AttachEvent(string id)
		{
			switch (id)
			{
				case VeldridSurface.DrawEvent:
					Control.DrawRectCalled += (sender, e) => Callback.OnDraw(Widget, EventArgs.Empty);
					break;
				default:
					base.AttachEvent(id);
					break;
			}
		}
	}
}
