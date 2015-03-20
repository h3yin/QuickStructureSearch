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
    private double binSize = 40.0;
    private int featureCount = (400);
    private int gapLength = 1;
 
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
    	
    	double[] features = new double[ (int)(this.featureCount*69*scale)];
    	
    	Point3d center = calculateCenter(coords);
    	
    	if (coords.length-this.length-1 <= 0) {
    		return features;
    	}
    	
    	for (int i = 0; i < coords.length; i++) {
    		if (coords[i] != null)
    		{
    			int temp = (int) center.distance(coords[i]);
    			int slope;
    			
    			
				String tempB, tempA, tempC, tempD, tempE, tempF, tempG, tempH, tempI, tempJ;
    			if( i-1*gapLength >= 0 && coords[i-1*gapLength] != null &&
    				i+1*gapLength < coords.length && coords[i+1*gapLength]!=null &&
    				i-2*gapLength >= 0 && coords[i-2*gapLength] != null &&
    				i+2*gapLength < coords.length && coords[i+2*gapLength]!=null &&
    				i-3*gapLength >= 0 && coords[i-3*gapLength] != null &&
     				i+3*gapLength < coords.length && coords[i+3*gapLength]!=null  &&
    				i-4*gapLength >= 0 && coords[i-4*gapLength] != null &&
     				i+4*gapLength < coords.length && coords[i+4*gapLength]!=null  &&
    				i-5*gapLength >= 0 && coords[i-5*gapLength] != null &&
     				i+5*gapLength < coords.length && coords[i+5*gapLength]!=null) 
    			{    				
    				if(center.distance(coords[i - 1*gapLength]) < temp)
    				{
    					tempB = "0";
    				}
    				else
    				{
    					tempB = "1";
    				}
    				
    				
    				if(center.distance(coords[i + 1*gapLength]) < temp)
    				{
    					tempA = "0";
    				}
    				else
    				{
    					tempA = "1";
    				}
    				
    				if(center.distance(coords[i - 2*gapLength]) < center.distance(coords[i - 1*gapLength]))
    				{
    					tempD = "0";
    				}
    				else
    				{
    					tempD = "1";
    				}
    				
    				
    				if(center.distance(coords[i + 2*gapLength]) < center.distance(coords[i + 1*gapLength]))
    				{
    					tempC = "0";
    				}
    				else
    				{
    					tempC = "1";
    				}
    				
    				if(center.distance(coords[i - 3*gapLength]) < center.distance(coords[i - 2*gapLength]))
    				{
    					tempF = "0";
    				}
    				else
    				{
    					tempF = "1";
    				}
    				
    				
    				if(center.distance(coords[i + 3*gapLength]) < center.distance(coords[i + 2*gapLength]))
    				{
    					tempE = "0";
    				}
    				else
    				{
    					tempE = "1";
    				}
    				
    				if(center.distance(coords[i - 4*gapLength]) < center.distance(coords[i - 3*gapLength]))
    				{
    					tempH = "0";
    				}
    				else
    				{
    					tempH = "1";
    				}
    				
    				
    				if(center.distance(coords[i + 4*gapLength]) < center.distance(coords[i + 3*gapLength]))
    				{
    					tempG = "0";
    				}
    				else
    				{
    					tempG = "1";
    				}
    				
    				slope = Integer.parseInt( tempD+tempC+tempB+tempA, 2);
    			}
    			else if (i-1*gapLength >= 0 && coords[i-1*gapLength] != null &&
        				i+1*gapLength < coords.length && coords[i+1*gapLength]!=null &&
        				i-2*gapLength >= 0 && coords[i-2*gapLength] != null &&
        				i+2*gapLength < coords.length && coords[i+2*gapLength]!=null &&
        				i-3*gapLength >= 0 && coords[i-3*gapLength] != null &&
         				i+3*gapLength < coords.length && coords[i+3*gapLength]!=null  &&
        				i-4*gapLength >= 0 && coords[i-4*gapLength] != null &&
         				i+4*gapLength < coords.length && coords[i+4*gapLength]!=null)
    			{
    				if(center.distance(coords[i - 1*gapLength]) < temp)
    				{
    					tempB = "0";
    				}
    				else
    				{
    					tempB = "1";
    				}
				
				
    				if(center.distance(coords[i + 1*gapLength]) < temp)
    				{
    					tempA = "0";
    				}
    				else
    				{
    					tempA = "1";
    				}
				
    				if(center.distance(coords[i - 2*gapLength]) < center.distance(coords[i - 1*gapLength]))
    				{
    					tempD = "0";
    				}
    				else
    				{
    					tempD = "1";
    				}
				
				
    				if(center.distance(coords[i + 2*gapLength]) < center.distance(coords[i + 1*gapLength]))
    				{
    					tempC = "0";
    				}
    				else
    				{
    					tempC = "1";
    				}
				
    				if(center.distance(coords[i - 3*gapLength]) < center.distance(coords[i - 2*gapLength]))
    				{
    					tempF = "0";
    				}
    				else
    				{
    					tempF = "1";
    				}
				
				
    				if(center.distance(coords[i + 3*gapLength]) < center.distance(coords[i + 2*gapLength]))
    				{
    					tempE = "0";
    				}
    				else
    				{
    					tempE = "1";
    				}
    				
    				if(center.distance(coords[i - 4*gapLength]) < center.distance(coords[i - 3*gapLength]))
    				{
    					tempH = "0";
    				}
    				else
    				{
    					tempH = "1";
    				}
    				
    				
    				if(center.distance(coords[i + 4*gapLength]) < center.distance(coords[i + 3*gapLength]))
    				{
    					tempG = "0";
    				}
    				else
    				{
    					tempG = "1";
    				
    				}
    				slope = Integer.parseInt( tempD+tempC+tempB+tempA, 2);
    			}	
    			else if(i-1*gapLength >= 0 && coords[i-1*gapLength] != null &&
        				i+1*gapLength < coords.length && coords[i+1*gapLength]!=null &&
        				i-2*gapLength >= 0 && coords[i-2*gapLength] != null &&
        				i+2*gapLength < coords.length && coords[i+2*gapLength]!=null &&
        				i-3*gapLength >= 0 && coords[i-3*gapLength] != null &&
         				i+3*gapLength < coords.length && coords[i+3*gapLength]!=null)
    			{
    				if(center.distance(coords[i - 1*gapLength]) < temp)
    				{
    					tempB = "0";
    				}
    				else
    				{
    					tempB = "1";
    				}
				
				
    				if(center.distance(coords[i + 1*gapLength]) < temp)
    				{
    					tempA = "0";
    				}
    				else
    				{
    					tempA = "1";
    				}
				
    				if(center.distance(coords[i - 2*gapLength]) < center.distance(coords[i - 1*gapLength]))
    				{
    					tempD = "0";
    				}
    				else
    				{
    					tempD = "1";
    				}
				
				
    				if(center.distance(coords[i + 2*gapLength]) < center.distance(coords[i + 1*gapLength]))
    				{
    					tempC = "0";
    				}
    				else
    				{
    					tempC = "1";
    				}
    				
    				if(center.distance(coords[i - 3*gapLength]) < center.distance(coords[i - 2*gapLength]))
    				{
    					tempF = "0";
    				}
    				else
    				{
    					tempF = "1";
    				}
				
				
    				if(center.distance(coords[i + 3*gapLength]) < center.distance(coords[i + 2*gapLength]))
    				{
    					tempE = "0";
    				}
    				else
    				{
    					tempE = "1";
    				}
    				slope = Integer.parseInt( tempD+tempC+tempB+tempA, 2);
    			}
    			else if (i-1*gapLength >= 0 && coords[i-1*gapLength] != null &&
        				i+1*gapLength < coords.length && coords[i+1*gapLength]!=null &&
        				i-2*gapLength >= 0 && coords[i-2*gapLength] != null &&
        				i+2*gapLength < coords.length && coords[i+2*gapLength]!=null)
    			{
    				if(center.distance(coords[i - 1*gapLength]) < temp)
    				{
    					tempB = "0";
    				}
    				else
    				{
    					tempB = "1";
    				}
				
				
    				if(center.distance(coords[i + 1*gapLength]) < temp)
    				{
    					tempA = "0";
    				}
    				else
    				{
    					tempA = "1";
    				}
				
    				if(center.distance(coords[i - 2*gapLength]) < center.distance(coords[i - 1*gapLength]))
    				{
    					tempD = "0";
    				}
    				else
    				{
    					tempD = "1";
    				}
				
				
    				if(center.distance(coords[i + 2*gapLength]) < center.distance(coords[i + 1*gapLength]))
    				{
    					tempC = "0";
    				}
    				else
    				{
    					tempC = "1";
    				}
    				
    				slope = 66;
    			}
    			else if (i-1*gapLength >= 0 && coords[i-1*gapLength] != null &&
        				i+1*gapLength < coords.length && coords[i+1*gapLength]!=null)
    			{
    				if(center.distance(coords[i - 1*gapLength]) < temp)
    				{
    					tempB = "0";
    				}
    				else
    				{
    					tempB = "1";
    				}
				
				
    				if(center.distance(coords[i + 1*gapLength]) < temp)
    				{
    					tempA = "0";
    				}
    				else
    				{
    					tempA = "1";
    				}
    				slope = 67;
    			}
    			else
    			{
    				slope = 68;
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
