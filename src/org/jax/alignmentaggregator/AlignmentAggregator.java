package org.jax.alignmentaggregator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;

public class AlignmentAggregator {
	
	private String _barcodeattribute;
	private int _forwardcorrection, _reversecorrection, _insertsizeoffset;
	
	public static void main(String[] args) {
		String[] parsedargs = new String[4];
		String barcodeattribute = "CB";
		String forwardcorrection = "4";
		String reversecorrection = "-5";
		
		int argidx = 0;
		for(int i = 0; i < args.length; i++) {
			if(i < args.length-1) {
				switch(args[i]) {
					case "--startbases":  
						forwardcorrection = args[i+1];
						i++;
						break;
					case "--endbases":  
						reversecorrection = args[i+1];
						i++;
					default:
						if(args[i].startsWith("-")) {
							argidx = 5;
						}
						else {
							parsedargs[argidx++] = args[i];
						}
						break;
				}
			}
			else {
				parsedargs[argidx++] = args[i];
			}
			
			if(argidx > 3) {
				break;
			}
		}
		
		if(argidx != 4) {
			System.out.println("Usage: bamfile clusterbarcodes outputdirectory");
			System.out.println("Options: --bambc     Bamfile attribute used for the barcode. (Default: \"CB\")");
			System.out.println("         --startbases The amount of bases add to the start position (must be 0 or positive). (Default: 4)");
			System.out.println("         --endbases The amount of bases to add to the end position (must be 0 or negative). (Default: -5)");

			System.exit(0);
		}
		
		
		long timestart = System.currentTimeMillis();
		
		String bamfile = parsedargs[0];
		String cellcluster = parsedargs[1];
		String outdir = parsedargs[2];
		
		AlignmentAggregator bcs = new AlignmentAggregator(barcodeattribute, Integer.parseInt(forwardcorrection), Integer.parseInt(reversecorrection));
		try {
			bcs.splitByCluster(bamfile, cellcluster, outdir);
		} catch (IOException e) {
			e.printStackTrace();
		}
		long timeend = System.currentTimeMillis();

		double duration = (double)(timeend-timestart)/1000;
		
		System.out.println("Completed in: "+Double.toString(duration)+" seconds.");
	}
	
	public AlignmentAggregator(String barcodeattr, int forwardoffset, int reverseoffset) {
		_barcodeattribute = barcodeattr;
		_forwardcorrection = forwardoffset;
		_reversecorrection = reverseoffset;
		_insertsizeoffset = forwardoffset-reverseoffset;
		
		if(_forwardcorrection < 0) {
			System.out.println("--startbases must be 0 or positive");
			System.exit(0);
		}
		
		if(_reversecorrection > 0) {
			System.out.println("--endbases must be 0 or negative");
			System.exit(0);
		}

	}
	
	public void splitByCluster(String bamfile, String cellclusterfile, String outdir) throws IOException{

		SamReaderFactory factory = SamReaderFactory.makeDefault()
	              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
	              .validationStringency(ValidationStringency.SILENT);
		
		final SamReader reader = factory.open(new File(bamfile));

		SAMFileHeader header = reader.getFileHeader();

		SAMRecordIterator it = reader.iterator();
		ClusterInfo ci = readCellClusters(cellclusterfile);
		TreeMap<String, String> cellclusters = ci.cellclustermap;
		TreeSet<String> clusters = ci.clusters;
		
		
		TreeMap<String, SAMFileWriter> writers = new TreeMap<String, SAMFileWriter>();
		SAMFileWriterFactory wfactory = new SAMFileWriterFactory();
		for(Iterator<String> clusterit = clusters.iterator(); clusterit.hasNext();) {
			String curcluster = clusterit.next();
			File f = new File(outdir+"/cluster_"+curcluster+".bam");
			if(f.exists()) {
				System.out.println("Error:"+f.getAbsolutePath()+" exists.  Exiting.");
				System.exit(1);
			}
			SAMFileWriter curwriter = wfactory.makeBAMWriter(header, true, new File(outdir+"/clusterbam_"+curcluster+".bam"));
			writers.put(curcluster, curwriter);
		}
		System.out.println("Splitting BAM files into clusters.");

		while(it.hasNext()){
			SAMRecord next = it.next();
			
			String barcode = getBarcode(next, _barcodeattribute);

			if(barcode == null) {
				continue;
			}
			
			if(cellclusters.containsKey(barcode)) {
				String cluster = cellclusters.get(barcode);
				
				if(cluster == null) {
					System.out.println("Missing cluster for barcode "+barcode);
				}
				else {
					SAMFileWriter curwriter = writers.get(cluster);
					boolean isnegative = next.getReadNegativeStrandFlag();
					
					if (isnegative) {
						//Reverse Strand
						curwriter.addAlignment(updateReverseRead(next));
					}
					else {
						//Forward Strand
						curwriter.addAlignment(updateForwardRead(next));
					}
					

				}
			}
		}

		it.close();
		reader.close();
		
		for(Iterator<SAMFileWriter> wit = writers.values().iterator(); wit.hasNext();) {
			SAMFileWriter next = wit.next();
			next.close();
		}
		
	}
	
	private SAMRecord updateForwardRead(SAMRecord r) {
		int alignmentstart = r.getAlignmentStart();
		
		if(alignmentstart == 0) {
			return r;
		}
		
		r.setAlignmentStart(alignmentstart+_forwardcorrection);
		r.setInferredInsertSize(r.getInferredInsertSize()-_insertsizeoffset);
		
		String readbases = r.getReadString();
		r.setReadString(readbases.substring(_forwardcorrection));
		
		String basequals = r.getBaseQualityString();
		r.setBaseQualityString(basequals.substring(_forwardcorrection));
		
		r.setCigar(getForwardTrimmedCigar(r.getCigar(), _forwardcorrection));
		
		
		return r;
	}
	
	private SAMRecord updateReverseRead(SAMRecord r) {
		int alignmentstart = r.getAlignmentStart();
		
		if(alignmentstart == 0) {
			return r;
		}
		
		r.setMateAlignmentStart(alignmentstart+_forwardcorrection);
		r.setInferredInsertSize(r.getInferredInsertSize()+_insertsizeoffset);
		
		String readbases = r.getReadString();
		r.setReadString(readbases.substring(0,readbases.length()+_reversecorrection));
		
		String basequals = r.getBaseQualityString();
		r.setBaseQualityString(basequals.substring(0, basequals.length()+_reversecorrection));
		
		r.setCigar(getReverseTrimmedCigar(r.getCigar(), -_reversecorrection));

		return r;
	}
	
	private Cigar getReverseTrimmedCigar(Cigar curcigar, int correction) {
		return getReverseCigar(getForwardTrimmedCigar(getReverseCigar(curcigar), correction));
	}
	
	private Cigar getReverseCigar(Cigar c) {
		List<CigarElement> newcigarlist = new LinkedList<CigarElement>();
		for(Iterator<CigarElement> cigit = c.iterator(); cigit.hasNext();) {
			((LinkedList<CigarElement>) newcigarlist).addFirst(cigit.next());
		}
		return new Cigar(newcigarlist);
	}
	
	private Cigar getForwardTrimmedCigar(Cigar curcigar, int correction) {
		Cigar newcigar = new Cigar();
		Iterator<CigarElement> cigit = curcigar.iterator();
		
		int remainder = correction;
		boolean first = true;
		while(cigit.hasNext() && remainder > 0) {
			CigarElement next = cigit.next();
			CigarOperator curop = next.getOperator();
			if(curop.consumesReadBases()) {
				int length = next.getLength();
				int newlength = length-remainder;
				
				if (newlength > 0) {
					newcigar.add(new CigarElement(length, curop));
					break;
				}
				else {
					remainder = remainder+newlength;
				}
				first = false;
			}
			else {
				if (first) {
					newcigar.add(next); //Keep only the first non-consumed elements
				}
			}
		}
		
		//Add the remainder of the cigar
		while(cigit.hasNext()) {
			newcigar.add(cigit.next());
		}
		
		return newcigar;
	}
	
	private ClusterInfo readCellClusters(String cellclusterfile) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(cellclusterfile));
		TreeMap<String, String> clustermap = new TreeMap<String, String>();
		TreeSet<String> clusters = new TreeSet<String>();
		while(br.ready()) {
			String line = br.readLine();
			try {

				String[] split = line.split("\t|,");
				String cellid = split[0];
				String clusterid = split[1];
				clustermap.put(cellid, clusterid);
				clusters.add(clusterid);
			}
			catch(Exception e) {
				System.out.println("Error:" +line);
			}
		}
		br.close();
		ClusterInfo rv = new ClusterInfo();
		rv.cellclustermap = clustermap;
		rv.clusters = clusters;
		return rv;
	}
	
	private class ClusterInfo {
		
		TreeMap<String, String> cellclustermap;
		TreeSet<String> clusters;
		
	}
	
	private String getBarcode(SAMRecord read, String attributeid) {
		List<SAMTagAndValue> attributes = read.getAttributes();
		for(Iterator<SAMTagAndValue> it = attributes.iterator(); it.hasNext();) {
			SAMTagAndValue next = it.next();
			if(next.tag.equals(attributeid)){
				return (String) next.value;
			}
		}
		return null;
	}
	
}
