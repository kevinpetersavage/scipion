package xmipp.viewer.particlepicker;

import ij.gui.ImageWindow;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.MouseInfo;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Stroke;
import java.awt.Toolkit;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;

import xmipp.ij.commons.XmippImageCanvas;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.Particle;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippResource;
import xmipp.viewer.particlepicker.training.model.ManualParticle;


public abstract class ParticlePickerCanvas extends XmippImageCanvas
{
	public final static BasicStroke dashedst = new BasicStroke(2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] { 10.0f }, 0.0f);
	public final static BasicStroke continuousst = new BasicStroke(2.0f);
	public final static BasicStroke activest = new BasicStroke(3.0f);
	//public final static BasicStroke activedst = new BasicStroke(3.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] { 10.0f }, 0.0f);
        
	protected ParticlePickerJFrame frame;
	protected Micrograph micrograph;
	protected ParticlePicker picker;
	protected XmippImageWindow iw;
	protected Cursor eraserCursor = Toolkit.getDefaultToolkit().createCustomCursor(XmippResource.getIcon("eraser.png").getImage(), new Point(0, 0), "Eraser");
    private Dimension  size;


    public void setMicrograph(Micrograph m)
	{
		micrograph = m;

	}
	


	protected boolean activemoved;

	private static boolean tongleSetSelected = false;
	
	public ParticlePickerCanvas(ParticlePickerJFrame frame)
	{
		this(frame, frame.getMicrograph());
	}

	public ParticlePickerCanvas(ParticlePickerJFrame frame, Micrograph micrograph)
	{
		super(micrograph.getImagePlus(frame.getParticlePicker().getFilters()));
		this.frame = frame;
		this.picker = frame.getParticlePicker();
		this.micrograph = micrograph;
		micrograph.runImageJFilters(picker.getFilters());
		addKeyListener(new KeyListener()
		{

			@Override
			public void keyTyped(KeyEvent arg0)
			{
			}

			@Override
			public void keyReleased(KeyEvent arg0)
			{
			}

			@Override
			public void keyPressed(KeyEvent e)
			{
				ParticlePickerCanvas.this.keyPressed(e);

			}
		});

		addMouseMotionListener(this);

		this.iw = new XmippImageWindow(getImage(), this, null);
		getFrame().displayZoom(getMagnification());
		iw.setTitle(getMicrograph().getName());

	}
	
	public void keyPressed(KeyEvent e)
	{
		PickerParticle active = getActive();
		int step = 1;
		int code = e.getKeyCode();
		int x = screenX(imp.getWidth()/2);
		int y = screenY(imp.getHeight()/2);
		if(active != null)
		{
			if (code == KeyEvent.VK_UP)
			{
				setActiveMoved(true);
				manageActive(active.getX(), active.getY() - step);
			}
			else if (code == KeyEvent.VK_DOWN)
			{
				setActiveMoved(true);
				manageActive(active.getX(), active.getY() + step);
			}
			else if (code == KeyEvent.VK_LEFT)
			{
				setActiveMoved(true);
				manageActive(active.getX() - step, active.getY());
			}
			else if (code == KeyEvent.VK_RIGHT)
			{
				setActiveMoved(true);
				manageActive(active.getX() + step, active.getY());
			}
		}
		if (code == KeyEvent.VK_SPACE)
		{
			getFrame().circlechb.setSelected(tongleSetSelected);
			getFrame().rectanglechb.setSelected(tongleSetSelected);
			tongleSetSelected = !tongleSetSelected;
		}
		
		
		else
			return;// do not repaint if not needed
		repaint();
	}

    protected void captureWindowSize(){
        size = this.imp.getWindow().getSize();
        this.imp.getWindow().setIgnoreRepaint(true);

    }
    protected void restoreWindowSize(){
        if (size != null) {
            this.imp.getWindow().setSize(size);
            this.imp.getWindow().setIgnoreRepaint(false);
        }
    }

	public void zoomIn(int sx, int sy)
	{

        captureWindowSize();

		super.zoomIn(sx, sy);
		getFrame().displayZoom(getMagnification());

        restoreWindowSize();
	}
	
	public void zoomOut(int sx, int sy)
	{
        captureWindowSize();

		super.zoomOut(sx, sy);
		getFrame().displayZoom(getMagnification());

        restoreWindowSize();
	}

	protected abstract Particle getLastParticle();

	protected abstract void manageActive(int x, int y);

	protected void setActiveMoved(boolean b)
	{
		activemoved = b;
	}

	protected void refresh()
	{
		getFrame().updateMicrographsModel();
		getFrame().setChanged(true);
		repaint();
	}

	public ImageWindow getIw()
	{
		return iw;
	}

	

	
	


	public void moveTo(PickerParticle p)
	{
		int width = (int) getSrcRect().getWidth();
		int height = (int) getSrcRect().getHeight();
		int x0 = p.getX() - width / 2;
		if (x0 < 0)
			x0 = 0;
		if (x0 + width > imp.getWidth())
			x0 = imp.getWidth() - width;
		int y0 = p.getY() - height / 2;
		if (y0 < 0)
			y0 = 0;
		if (y0 + height > imp.getHeight())
			y0 = imp.getHeight() - height;
		Rectangle r = new Rectangle(x0, y0, width, height);
		if (!getSrcRect().contains(r))
		{
			setSourceRect(r);
			repaint();
		}
	}

	public void mouseEntered(MouseEvent e)
	{
		super.mouseEntered(e);
		setCustomCursor(e);
	}


	
	public void mouseExited(MouseEvent e)
	{
		super.mouseExited(e);
		setCursor(new Cursor(Cursor.DEFAULT_CURSOR));

	}
	
	
	protected void setCustomCursor(MouseEvent e)
	{
		if (getFrame().isPickingAvailable(e))
		{
			if (!getFrame().isEraserMode())
				setCursor(crosshairCursor);
			else 
				setCursor(eraserCursor);
		}
	}

	public abstract void refreshActive(Particle p);
	

	public abstract PickerParticle getActive();

	
	public abstract ParticlePickerJFrame getFrame();


	protected void drawShape(Graphics2D g2, ManualParticle p, boolean all, Stroke stroke)
	{
		drawShape(g2, p.getX(), p.getY(), picker.getSize(), all, stroke);
	}
	
	protected void drawShape(Graphics2D g2, ManualParticle p, boolean all)
	{
		drawShape(g2, p, all, continuousst);
	}
	
	protected void drawShape(Graphics2D g2, int x, int y, int size, boolean all)
	{
		drawShape(g2, x, y, size, all, continuousst);
				
	}
	protected void drawShape(Graphics2D g2, int x, int y, int size, boolean all, Stroke stroke, Color color)
	{
		g2.setColor(color);
		drawShape(g2, x, y, size, all, stroke);
	}
	protected void drawShape(Graphics2D g2, int x, int y, int size, boolean all, Stroke stroke)
	{

		g2.setStroke(stroke);
		int length = (int) (size * magnification);
		int radius = (int) (size / 2. * magnification);
		x = getXOnImage(x);
		y = getYOnImage(y);
		int distance = Math.min(10, (int) (radius/4. * magnification));

		if (getFrame().isShapeSelected(Shape.Rectangle) || all)
			g2.drawRect(x - radius, y - radius, length, length);
		if (getFrame().isShapeSelected(Shape.Circle) || all)
			g2.drawOval(x - radius, y - radius, length, length);
		if (getFrame().isShapeSelected(Shape.Center) || all)
		{
			g2.setStroke(activest);
			g2.drawLine(x, y - distance, x, y + distance);
			g2.drawLine(x + distance, y, x - distance, y);
		}
	}


	protected void drawLine(double alpha, Graphics2D g2)
	{
		int width = imp.getWidth();
		int height = imp.getHeight();
		double m = 0;
		double x1, y1, x2, y2;
		if (alpha != Math.PI / 2)
		{
			m = Math.tan(alpha - Math.PI / 2);

			double y = height / 2.f;
			double x = y / m;

			if (Math.abs(x) > width / 2.f)// cuts in image sides
			{
				x1 = width;// on image
				y1 = getYCutOnImage(m, width / 2.f);
				x2 = 0;
				y2 = getYCutOnImage(m, -width / 2.f);
			}
			else
			// cuts in image top and bottom
			{
				y1 = 0;
				x1 = getXCutOnImage(m, height / 2.f);
				y2 = height;
				x2 = getXCutOnImage(m, -height / 2.f);
			}
		}
		else
		{
			x1 = 0;
			y1 = y2 = height / 2.f;
			x2 = width;

		}
		Color ccolor = g2.getColor();
		g2.setColor(Color.yellow);
		g2.drawLine((int) (x1 * magnification), (int) (y1 * magnification), (int) (x2 * magnification), (int) (y2 * magnification));
		g2.setColor(ccolor);
	}
	
	
	public int getXOnImage(int x)
	{
		int x0 = (int) getSrcRect().getX();
		return (int) ((x - x0) * magnification);
	}
	
	public int getYOnImage(int y)
	{
		int y0 = (int) getSrcRect().getY();
		return (int) ((y - y0) * magnification);
	}

	private double getYCutOnImage(double m, double x)
	{
		int height = imp.getHeight();
		return height / 2.f - m * x;
	}


	private double getXCutOnImage(double m, double y)
	{
		int width = imp.getWidth();
		return y / m + width / 2.f;
	}
	
	public void updateMicrograph()
	{
		Micrograph m = getFrame().getMicrograph();
		updateMicrograph(m);
	}

	public void updateMicrograph(Micrograph m)
	{
		getMicrograph().releaseImage();
		setMicrograph(m);
		imp = m.getImagePlus(picker.getFilters());
		m.runImageJFilters(picker.getFilters());
		refreshActive(null);
		Rectangle previousRect = (Rectangle)srcRect.clone();
		int previousDstWidth = dstWidth;
		int previousDstHeight = dstHeight;
		if(iw.isClosing())
			iw = new XmippImageWindow(getImage(), this, null);
		iw.setImage(imp);
		double zoom = getParticlePicker().getZoom();
		int imageWidth = imp.getWidth();
		int imageHeight = imp.getHeight();
		if(zoom != -1)
		{
			
			setDrawingSize(previousDstWidth, previousDstHeight);
			if (previousRect.width < imageWidth && previousRect.height < imageHeight)
				setSourceRect(previousRect);
			setMagnification(zoom);
			iw.pack();
		}
		else
			iw.maximize();
		getFrame().displayZoom(getMagnification());
		iw.setTitle(getMicrograph().getName());
		
	}



	public abstract Micrograph getMicrograph();

	protected void moveActiveParticle(int x, int y)
	{
		PickerParticle active = getActive();
		if (active == null)
			return;
		try
		{
			active.setPosition(x, y);
		}
		catch (Exception e)
		{
			XmippDialog.showError(getFrame(), e.getMessage());
		}
		if (getFrame().getParticlesJDialog() != null)
			active.getParticleCanvas(getFrame()).repaint();
	}

	

	public void paint(Graphics g)
	{
		Graphics offgc;//Off screen graphic
		Image offscreen = null;//Off screen image
		Dimension d = getSize();

		// create the offscreen buffer and associated Graphics
		offscreen = createImage(d.width, d.height);
		offgc = offscreen.getGraphics();
		
		super.paint(offgc);//super paint in off screen
		//my paint in offscreen
		Graphics2D g2 = (Graphics2D) offgc;
		doCustomPaint(g2);
		//drawing offscreen image
		g.drawImage(offscreen, 0, 0, this);
	}
	

	protected abstract void doCustomPaint(Graphics2D g2);
        
   
         
    public abstract ParticlePicker getParticlePicker();

}
