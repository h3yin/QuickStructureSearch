package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

/*
 * Interface for all chain smoothers
 * Contain function header for all getSmoothedPoints
 */
public interface ChainSmoother {
	
	/*Takes in a list of 3d points and returns a list of smoothed 3d points*/
	public Point3d [] getSmoothedPoints(Point3d [] points);

}
