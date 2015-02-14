package org.rcsb.structuralSimilarity;

import java.io.Serializable;
import javax.vecmath.Point3d;

/*
 * Class which implements the RogenChainSmoother
 * Will be called by ChainSmootherMapper when needed
 */

public class RogenChainSmoother implements ChainSmoother, Serializable{
	
	private int iterations = 0;
	
	public RogenChainSmoother(int it)
	{
		this.iterations = it;
	}
	
	public Point3d [] getSmoothedPoints(Point3d [] points)
	{
		double[] x = new double[points.length];
		double[] y = new double[points.length];
		double[] z = new double[points.length];
				
		for (int i = 0; i < points.length; i++) {
		//	System.out.println("Hey i = " + i);
		//	System.out.println("Hey x[i] = " + x[i]);
		//	System.out.println("Hey points[i] = " + points[i]);
			
			/*TODO: fix problems with null values*/
			if(points[i] == null)
			{
				x[i] = 0;
				y[i] = 0;
				z[i] = 0;
			}
			else
			{
				x[i] = points[i].x;
				y[i] = points[i].y;
				z[i] = points[i].z;
			}
		}
		
		for (int i = 0; i < iterations; i++) {
			x = smoothRogen(x);
			y = smoothRogen(y);
			z = smoothRogen(z);
		}

		Point3d[] smoothedPoints = new Point3d[x.length];
		for (int i = 0; i < x.length; i++) {
			smoothedPoints[i] = new Point3d(x[i],y[i],z[i]);
		}
		
		return smoothedPoints;
	}
	private double[] smoothRogen(double[] x) {	    	

	/**
	 * P. Rogen, Evaluating protein structure descriptors and tuning Gauss integral
	 * based descriptors, J. Phys.: Condens. Matter (2005) 17, S1523-S1538
	 * @param x
	 * @return
	 */
	
		double[] y = new double[x.length];
        double a = 2.4;
        double b = 2.1;
        
		y[0] = x[0];
		y[1] = (a*x[0] + b*x[1] + a*x[2])/(a+a+b);
		for (int i = 2; i < x.length-2; i++) {
			y[i] = (x[i-2]+ a*x[i-1] + b*x[i] + a*x[i+1] + x[i+2])/(2+a+a+b);
		}
		y[x.length-2] =  (a*x[x.length-3] + b*x[x.length-2] + a*x[x.length-1])/(a+a+b);
		y[x.length-1] = x[x.length-1];

		return y;
	}

}
