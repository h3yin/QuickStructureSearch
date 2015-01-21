package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;

/**
* This class map a tuple <String, Point3d[]> to another tuple <String, Point3d[]> with
* smoothed points using the ChainSmoother passed into the constructor.
* 
* @author Peter Rose
*/
public class ChainSmootherMapper implements PairFunction<Tuple2<String, Point3d[]>,String, Point3d[]> {
	private static final long serialVersionUID = 1L;
//
//	public ChainSmootherMapper(ChainSmoother smoother) {
//		this.smoother = smoother;
//	}

	/**
	 * Maps a <String, Point3d[]> pair by a <String, Point3d[]> pair with smoothed points.
	 */
	@Override
	public Tuple2<String, Point3d[]> call(Tuple2<String, Point3d[]> tuple) throws Exception {	
//	    Point3d[] smoothPoints = smoother.getSmoothedPoints(tuple._2);
		Point3d[] smoothPoints = new Point3d[0];
		return new Tuple2<String,Point3d[]>(tuple._1, smoothPoints);
	}
}
