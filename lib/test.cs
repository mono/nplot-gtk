#define GTK
#if GTK
using Gtk;
using NPlot.Gtk;
#endif

using NPlot;
using System.Drawing;
using System.Drawing.Imaging;

class X {
	public static float[] makeDaub( int len )
	{
		float[] daub4_h = 
			{ 0.482962913145f, 0.836516303737f, 0.224143868042f, -0.129409522551f };
		
		float[] daub4_g = 
			{ -0.129409522551f, -0.224143868042f, 0.836516303737f, -0.482962913145f };
		
		float[] a = new float[len];
		a[8] = 1.0f;
		float[] t;
		
		int ns = 4;  // number smooth
		while ( ns < len/2 ) {
			t = (float[])a.Clone();
			
			ns *= 2;
			
			for ( int i=0; i<(ns*2); ++i ) 
				a[i] = 0.0f;
			
			// wavelet contribution
			for ( int i=0; i<ns; ++i ) 
				for ( int j=0; j<4; ++j ) 
					a[(2*i+j)%(2*ns)] += daub4_g[j] * t[i+ns];

			// smooth contribution
			for ( int i=0; i<ns; ++i ) 
				for ( int j=0; j<4; ++j ) 
					a[(2*i+j)%(2*ns)] += daub4_h[j]*t[i];
		}
		return a;
	}

	static public void PlotWavelet(IPlotSurface2D plot)
	{	
		plot.Clear();
		
		// Create a new line plot from array data via the ArrayAdapter class.
		LinePlot lp = new LinePlot();
		lp.DataSource = makeDaub(256);
		lp.Color = Color.Green;
		
		Grid myGrid = new Grid();
		myGrid.VerticalGridType = Grid.GridType.Fine;
		myGrid.HorizontalGridType = Grid.GridType.Coarse;
		plot.Add(myGrid);

		// And add it to the plot surface
		plot.Add( lp );
		plot.Title = "Reversed / Upside down Daubechies Wavelet";

		// Ok, the above will produce a decent default plot, but we would like to change
		// some of the Y Axis details. First, we'd like lots of small ticks (10) between 
		// large tick values. Secondly, we'd like to draw a grid for the Y values. To do 
		// this, we create a new LinearAxis (we could also use Label, Log etc). Rather than
		// starting from scratch, we use the constructor that takes an existing axis and
		// clones it (values in the superclass Axis only are cloned). PlotSurface2D
		// automatically determines a suitable axis when we add plots to it (merging
		// current requirements with old requirements), and we use this as our starting
		// point. Because we didn't specify which Y Axis we are using when we added the 
		// above line plot (there is one on the left - YAxis1 and one on the right - YAxis2)
		// PlotSurface2D.Add assumed we were using YAxis1. So, we create a new axis based on
		// YAxis1, update the details we want, then set the YAxis1 to be our updated one.
		LinearAxis myAxis = new LinearAxis( plot.YAxis1 );
		myAxis.NumberOfSmallTicks = 10;
		plot.YAxis1 = myAxis;
	
		// We would also like to modify the way in which the X Axis is printed. This time,
		// we'll just modify the relevant PlotSurface2D Axis directly. 
		plot.XAxis1.WorldMax = 100.0f;
		
		plot.PlotBackColor = Color.Beige;
		plot.XAxis1.Reversed = true;
		plot.YAxis1.Reversed = true;
		
		// Force a re-draw of the control. 
		//plot.Refresh();
	}

	static public void PlotTest(IPlotSurface2D plotSurface)
	{
		plotSurface.Clear();
		
		plotSurface.Title = "My Title";
		
		double[] a = {0, 2, 1, 4, double.NaN, double.NaN, 5, 8, 7, 9};
		LinePlot lp = new LinePlot();
		lp.DataSource = a;
		lp.Label = "My Label";
		
		plotSurface.Add( lp );
		
		// plotSurface.Add( lp );
		
		plotSurface.Legend = new Legend();
		plotSurface.Legend.NeverShiftAxes = true;
		plotSurface.Legend.HorizontalEdgePlacement = Legend.Placement.Inside;
		plotSurface.Legend.VerticalEdgePlacement = Legend.Placement.Inside;
		plotSurface.Legend.XOffset = -10;
		plotSurface.Legend.YOffset = 10;
		//plotSurface.AddAxesConstraint( new AxesConstraint.EqualSpacing() );
		
		((LinearAxis)plotSurface.XAxis1).Offset = 10.0;
		((LinearAxis)plotSurface.XAxis1).Scale = 27.0;
		//((LinearAxis)plotSurface.XAxis1).TicksIndependentOfPhysicalExtent = true;
		//((LinearAxis)plotSurface.YAxis1).TicksIndependentOfPhysicalExtent = true;
		
		//plotSurface.Refresh();
	}
	
	static void Main ()
	{
#if GTK
		Application.Init ();

		Window w = new Window ("Test");
		w.DeleteEvent += delegate {
			Application.Quit ();
		};
		
		NPlot.Gtk.PlotSurface2D plot = new NPlot.Gtk.PlotSurface2D ();

		PlotTest (plot);
		
		plot.Show ();
		w.Add (plot);
		w.ShowAll ();

		Application.Run ();
#else

		NPlot.PlotSurface2D s = new NPlot.PlotSurface2D ();
		Bitmap b = new Bitmap (1000, 1000);
		Graphics g = Graphics.FromImage (b);
		g.FillRectangle  (Brushes.White, 0, 0, 1000, 1000);
		Rectangle bounds = new Rectangle (0, 0, 1000, 1000);
		PlotTest (s);
		s.Draw (g, bounds);
		b.Save ("file.png", ImageFormat.Png);
#endif
	}
}
