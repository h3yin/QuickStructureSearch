package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.UnsupportedEncodingException;

import javax.vecmath.Point3d;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;

public class SavitzkyGolayChainSmoother implements ChainSmoother, Serializable {

	private int iterations = 0;
	
	public SavitzkyGolayChainSmoother(int it)
	{
		this.iterations = it;
	}
	
	
	/**
	 * http://en.wikipedia.org/wiki/Savitzky–Golay_filter
	 * Note, the smoothed curve will have fewer points in its current implementation
	 * @param points
	 * @param iterations
	 * @return
	 */
	@Override
	public Point3d[] getSmoothedPoints(Point3d[] points) {
		int[] coefficients = {1,1,1,1};
		int norm = 4;
		return smoothSavitzkyGolay(points, coefficients, norm, iterations);
	}
	/**
	 * http://en.wikipedia.org/wiki/Savitzky–Golay_filter
	 * Note, the smoothed curve will have fewer points in its current implementation
	 * @param points
	 * @param iterations
	 * @return
	 */
	public Point3d[] smoothSavitzkyGolayCubic7Point(Point3d[] points, int iterations) {	
		int[] coefficients = {-2,3,6,7,6,3,-2};
		int norm = 21;
		return smoothSavitzkyGolay(points, coefficients, norm, iterations);
	}
	
	public Point3d[] smoothSavitzkyGolay(Point3d[] points, int[] coefficients, int norm, int iterations) {	
		double[] x = new double[points.length];
		double[] y = new double[points.length];
		double[] z = new double[points.length];
		
		for (int i = 0; i < points.length; i++) {
			
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
			x = smoothSavitzkyGolay(x, coefficients, norm);
			y = smoothSavitzkyGolay(y, coefficients, norm);
			z = smoothSavitzkyGolay(z, coefficients, norm);
		}

		// there will be fewer points after smoothing
		Point3d[] smoothedPoints = new Point3d[x.length];
		for (int i = 0; i < x.length; i++) {
			smoothedPoints[i] = new Point3d(x[i],y[i],z[i]);
		}
		
		PrintWriter writer = null;
		try {
			writer = new PrintWriter("NumAtoms.txt", "UTF-8");
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		writer.println("The number of atoms is: " + points.length);
		
		writer.close();
		
		return smoothedPoints;
	}
	
    private double[] smoothSavitzkyGolay(double[] x, int[] coefficients, int norm) {
    	int n = coefficients.length - 1;
    	int m = x.length - n;
    	
    	double[] y = new double[m];
    	double[] p = new double[coefficients.length];
    	
        for (int i = 1; i < coefficients.length; i++) {
        	p[i] = x[i-1];	
        }
        for (int i = 0; i < m; i++) {
        	for (int k = 0; k < n; k++) {
        		p[k]= p[k+1];
        	}
        	p[n] = x[i+n];
        	double sum = 0;
        	for (int k = 0; k < coefficients.length; k++) {
        		sum += coefficients[k] * p[k];
        	}
        	y[i] = sum/norm;
        }
        return y;
    }

	private AtomCache initializeCache() {
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();
		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(true);
		params.setLoadChemCompInfo(true);
		cache.setFileParsingParams(params);
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		return cache;
	}

}
