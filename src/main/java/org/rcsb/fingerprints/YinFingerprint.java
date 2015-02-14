package org.rcsb.fingerprints;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.UnsupportedEncodingException;

import javax.vecmath.Point3d;

public class YinFingerprint implements GenericFingerprint, Serializable{
	
	private static final long serialVersionUID = 1L;
	private double MAX_AVE_CA_CA_DISTANCE = 3.3;
    private int length = 9;
    private double binSize = 1.0;
    private int featureCount = (200);
 
    /**
     * Default constructor uses default parameters
     */
    public YinFingerprint() {}
    
    /**
     * Constructor with all parameters
     * @param length fragment length
     */
    public YinFingerprint (int length, double binSize) {
        this.length = length;
        this.binSize = binSize;
        this.featureCount = (int)(this.length * this.MAX_AVE_CA_CA_DISTANCE/this.binSize);
	}

    public String getName() {
    	return this.getClass().getSimpleName() + "_L" + this.length + "B" + this.binSize;
    }

    /*
     * Returns a fingerprint for the given chain. 
     * @param coords coordinates of a macromolecule fragment
     * @return fingerprint
     */
    public double[] getFingerprint(Point3d[] coords) {
    	double scale = 1/this.binSize;
    	
    	double[] features = new double[ (int)(this.featureCount*4*scale)];
    	
    	Point3d center = calculateCenter(coords);
    	
    	if (coords.length-this.length-1 <= 0) {
    		return features;
    	}
    	
    	for (int i = 0; i < coords.length; i++) {
    		if (coords[i] != null)
    		{
    			int temp = (int) center.distance(coords[i]);
    			int slope;
    			
    			String tempB, tempA;
    			
    			if( i-1 >= 0 && coords[i-1] != null &&
    				i+1 < coords.length && coords[i+1]!=null)
    			{
    				if(center.distance(coords[i - 1]) < temp)
    				{
    					tempB = "0";
    				}
    				else
    				{
    					tempB = "1";
    				}
    				
    				
    				if(center.distance(coords[i + 1]) < temp)
    				{
    					tempA = "0";
    				}
    				else
    				{
    					tempA = "1";
    				}
    				
    				slope = Integer.parseInt(tempB+tempA, 2);
    			}
    			else
    			{
    				slope = 0;
    			}
    			
    			try{
    				features[(int) ((slope*this.featureCount + temp)*scale)]++;
    			}
    			catch (Exception e)
    			{
    				e.printStackTrace();
    				System.out.println("Index is: " + (int) ((slope*this.featureCount + temp)*scale));
    			}
    			
    		}
    		
    	}
    	
    	PrintWriter writer;
		try {
			writer = new PrintWriter("featureVector.txt", "UTF-8");
			
			writer.println("This is a feature vector");
		 	writer.println("The first line");
	    	writer.println("The second line");
	    	
	    	for(double elem: features)
	    	{
	    		writer.print( " " + elem + " ");
	    	}
	    	
	    	writer.println("\n\n\n");
	    	writer.println("Size is: " + features.length);
	    	writer.println("\n\n\n");
	    	
	    	writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
   
		return features;
    }
    
    private Point3d calculateCenter (Point3d [] coords)
    {
    	/* Total used to calculate average point (center)*/
    	double sumx = 0.0;
    	double sumy = 0.0;
    	double sumz = 0.0;
    	
    	double count = 0.0;
    	
    	for (int i = 0 ; i < coords.length; i++)
    	{
    		if(coords[i] != null)
    		{
    			sumx += coords[i].x;
    			sumy += coords[i].y;
    			sumz += coords[i].z;
    			
    			count += 1.0;
    		}

    	}
		
		return new Point3d(sumx/count, sumy/count, sumz/count); 	
    }
}
