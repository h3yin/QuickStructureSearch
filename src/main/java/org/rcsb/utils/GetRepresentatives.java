package org.rcsb.utils;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Arrays;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.align.xml.RepresentativeXMLConverter;

/**
 * TODO Move this to {@link Representatives}.
 */
public class GetRepresentatives {

	private static String clusterUrl = "http://www.rcsb.org/pdb/rest/representatives?cluster=";
	private static String allUrl = "http://www.rcsb.org/pdb/rest/getCurrent/";

	// available sequence clusters
	private static List<Integer> seqIdentities = Arrays.asList(30, 40, 50, 70, 90, 95, 100);

	/**
	 * Returns a representative set of PDB protein chains at 40% sequence 
	 * identity cutoff.
	 * @param sequenceIdentity
	 * @return PdbChainKey set of representatives
	 * @deprecated
	 */
	@Deprecated
	public static SortedSet<StructureName> getRepresentatives() {
        return getRepresentatives(40);
	}

	/**
	 * Returns a representative set of PDB protein chains at the specified sequence 
	 * identity cutoff. See http://www.pdb.org/pdb/statistics/clusterStatistics.do
	 * for more information.
	 * @param sequenceIdentity sequence identity threshold
	 * @return PdbChainKey set of representatives
	 */
	public static SortedSet<StructureName> getRepresentatives(int sequenceIdentity) {
		SortedSet<StructureName> representatives = new TreeSet<StructureName>();

		if (!seqIdentities.contains(sequenceIdentity)) {
			System.err.println("Error: representative chains are not available for %sequence identity: "
							+ sequenceIdentity);
			return representatives;
		}


		try {

			URL u = new URL(clusterUrl + sequenceIdentity);

			InputStream stream = HTTPConnectionTools.getInputStream(u, 60000);

			String xml = null;

			if (stream != null) {
				xml = JFatCatClient.convertStreamToString(stream);

				SortedSet<String> reps = RepresentativeXMLConverter.fromXML(xml);

				for (String s : reps) {
					StructureName k = new StructureName(s);
					representatives.add(k);
				}

			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		return representatives;
	}
	
	/**
	 * Returns the current list of all PDB IDs.
	 * @return PdbChainKey set of all PDB IDs.
	 */
	public static SortedSet<String> getAll() {
		SortedSet<String> representatives = new TreeSet<String>();

		try {

			URL u = new URL(allUrl);

			InputStream stream = HTTPConnectionTools.getInputStream(u, 60000);

			if (stream != null) {
				BufferedReader reader = new BufferedReader(
						new InputStreamReader(stream));

				String line = null;

				while ((line = reader.readLine()) != null) {
					int index = line.lastIndexOf("structureId=");
					if (index > 0) {
						representatives.add(line.substring(index + 13, index + 17));
					}
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		return representatives;
	}
}
