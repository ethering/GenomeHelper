/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.vcf;

import edu.unc.genomics.VCFEntry;
import edu.unc.genomics.io.VCFFileReader;
import edu.unc.genomics.io.VCFFileWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.stat.descriptive.*;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

/**
 *
 * @author ethering
 */
public class VCFParser
{

    //test class to print some vcf info
    public void printVars(File in) throws IOException
    {
        VCFFileReader reader = new VCFFileReader(in.toPath());
        for (VCFEntry vcf : reader)
        {

            System.out.println(vcf.getChr() + "\t" + vcf.getStart() + "\t" + vcf.getStop() + "\t" + vcf.getRef() + "\t");

            String[] alts = vcf.getAlt();
            Map<String, String> info = (Map<String, String>) vcf.getInfo();
            for (Map.Entry pairs : info.entrySet())
            {
                System.out.println(pairs.getKey() + " = " + pairs.getValue());
            }
            List<String[]> genotypes = vcf.getGenotypes();
            for (String alt : alts)
            {
                System.out.println(alt + " ");
            }
            for (String[] gts : genotypes)
            {
                System.out.println(Arrays.toString(gts));
            }
            System.out.println("");
        }
    }

    public HashMap calculateGATKParams(File vcfFile) throws IOException
    {
        DescriptiveStatistics snpQualStats = new DescriptiveStatistics();
        DescriptiveStatistics depthStats = new DescriptiveStatistics();
        DescriptiveStatistics mappingStats = new DescriptiveStatistics();
        Percentile snpQualP = new Percentile();
        Percentile depthP = new Percentile();
        Percentile mappingP = new Percentile();

        try (VCFFileReader reader = new VCFFileReader(vcfFile.toPath()))
        {
            for (VCFEntry vcf : reader)
            {
                Map<String, String> info = (Map<String, String>) vcf.getInfo();
                snpQualStats.addValue(vcf.getQual());

                String stringCov = info.get("DP");
                depthStats.addValue(Double.parseDouble(stringCov));

                String stringMappingQual = info.get("MQ");
                if (!Double.isNaN(Double.parseDouble(stringMappingQual)))
                {
                    mappingStats.addValue(Double.parseDouble(stringMappingQual));
                }
            }
        }
        double[] snpQualArray = snpQualStats.getValues();
        double[] depthArray = depthStats.getValues();
        double[] mappingArray = mappingStats.getValues();
        snpQualP.setData(snpQualArray);
        depthP.setData(depthArray);
        mappingP.setData(mappingArray);

        //System.out.println("\nSNP Quality (QUAL)");
        //System.out.println("Mean = " + snpQualStats.getMean());
        //System.out.println("SD = " + snpQualStats.getStandardDeviation());
        //System.out.println("Lower 5 percentile = " + snpQualP.evaluate(5));
        Double lowQUAL = snpQualP.evaluate(1);
        System.out.println("lowQUAL\t"+lowQUAL);

        //System.out.println("\nCoverage (DP)");
        //System.out.println("Mean = " + depthStats.getMean());
        //System.out.println("SD = " + depthStats.getStandardDeviation());
        //System.out.println("Lower 5 percentile = " + depthP.evaluate(5));
        Double lowDP = depthP.evaluate(1);
        System.out.println("lowDP\t"+lowDP);

        //System.out.println("\nMapping quality (MQ)");
        //System.out.println("Mean = " + mappingStats.getMean());
        //System.out.println("SD = " + mappingStats.getStandardDeviation());
        //System.out.println("Lower 5 percentile = " + mappingP.evaluate(5));
        Double lowMQ = mappingP.evaluate(1);
        System.out.println("lowMQ\t"+lowMQ);
        //System.out.println("Returning HashMap for 'lowQual, lowDP and lowMQ'");
        
        HashMap<String, Double> params = new HashMap(3);
        params.put("lowQUAL", lowQUAL);
        params.put("lowDP", lowDP);
        params.put("lowMQ", lowMQ);
        return params;
    }

    public HashMap calculateGATKParams(File vcfFile, int maxRecords) throws IOException
    {
        DescriptiveStatistics snpQualStats = new DescriptiveStatistics();
        DescriptiveStatistics depthStats = new DescriptiveStatistics();
        DescriptiveStatistics mappingStats = new DescriptiveStatistics();
        Percentile snpQualP = new Percentile();
        Percentile depthP = new Percentile();
        Percentile mappingP = new Percentile();
        int recordCount = 0;
        try (VCFFileReader reader = new VCFFileReader(vcfFile.toPath()))
        {
            Iterator it = reader.iterator();
            while (it.hasNext() && recordCount < maxRecords)
            {
                VCFEntry vcf = (VCFEntry) it.next();

                Map<String, String> info = (Map<String, String>) vcf.getInfo();
                snpQualStats.addValue(vcf.getQual());

                String stringCov = info.get("DP");
                depthStats.addValue(Double.parseDouble(stringCov));

                String stringMappingQual = info.get("MQ");
                if (!Double.isNaN(Double.parseDouble(stringMappingQual)))
                {
                    mappingStats.addValue(Double.parseDouble(stringMappingQual));
                }

                recordCount++;
            }
        }

        double[] snpQualArray = snpQualStats.getValues();
        double[] depthArray = depthStats.getValues();
        double[] mappingArray = mappingStats.getValues();
        snpQualP.setData(snpQualArray);
        depthP.setData(depthArray);
        mappingP.setData(mappingArray);

        System.out.println("\nSNP Quality (QUAL)");
        System.out.println("Mean = " + snpQualStats.getMean());
        System.out.println("SD = " + snpQualStats.getStandardDeviation());
        //System.out.println("Lower 5 percentile = " + snpQualP.evaluate(5));
        Double lowQUAL = snpQualP.evaluate(1);
        System.out.println("Lower 1 percentile = " + lowQUAL);

        System.out.println("\nCoverage (DP)");
        System.out.println("Mean = " + depthStats.getMean());
        System.out.println("SD = " + depthStats.getStandardDeviation());
        //System.out.println("Lower 5 percentile = " + depthP.evaluate(5));
        Double lowDP = depthP.evaluate(1);
        System.out.println("Lower 1 percentile = " + lowDP);

        System.out.println("\nMapping quality (MQ)");
        System.out.println("Mean = " + mappingStats.getMean());
        System.out.println("SD = " + mappingStats.getStandardDeviation());
        //System.out.println("Lower 5 percentile = " + mappingP.evaluate(5));
        Double lowMQ = mappingP.evaluate(1);
        System.out.println("Lower 1 percentile = " + lowMQ);

        System.out.println("Returning HashMap for 'lowQual, lowDP and lowMQ'");
        HashMap<String, Double> params = new HashMap(3);
        params.put("lowQUAL", lowQUAL);
        params.put("lowDP", lowDP);
        params.put("lowMQ", lowMQ);
        return params;

    }

    public void printHeterozygotes(File vcfFile) throws IOException
    {
        VCFFileReader reader = new VCFFileReader(vcfFile.toPath());
        //loop thru VCFEntrys
        for (VCFEntry vcf : reader)
        {
            String[] alts = vcf.getAlt();
            if (alts.length > 1)
            {
                boolean allsnps = true;
                for (String alt : alts)
                {
                    if (alt.length() > 1)
                    {
                        allsnps = false;
                    }
                }

                if (allsnps)
                {
//                    System.out.print(vcf.getChr() + "\t" + vcf.getStart());
//                    for (String alt : alts)
//                    {
//                        System.out.print("\t" + alt);
//                    }
//                    System.out.println();
                    System.out.println(vcf.toOutput());
                }

            }
        }
    }

    public double CountHeterozygousSites(ArrayList<File> vcfs) throws IOException
    {
        HashMap<String, ArrayList<Integer>> hets = new HashMap();
        double heterozygousSites = 0;
        //parser.printVars(vcf1);
        for (File f : vcfs)
        {
            //open each file as a VCFReader
            VCFFileReader reader = new VCFFileReader(f.toPath());
            //loop thru VCFEntrys
            for (VCFEntry vcf : reader)
            {
                SnpLoci snp;
                String ref = vcf.getRef();
                String[] alts = vcf.getAlt();

                int depth = Integer.parseInt(vcf.getInfo().get("DP"));

                //we only want single nt snps with a depth > 10
                if (ref.length() == 1 && depth >= 10)
                {
                    //look at the alts - are they all singles?
                    boolean allsnps = true;
                    for (String alt : alts)
                    {
                        if (alt.length() > 1)
                        {
                            allsnps = false;
                        }
                    }
                    if (allsnps)
                    {
                        String chr = vcf.getChr();
                        int start = vcf.getStart();
                        List<String[]> genotypes = vcf.getGenotypes();
                        String[] genotypeArray = genotypes.get(0);
                        String genotype = genotypeArray[0];

                        //System.out.println(genotype);
                        if ((genotype.equals("0/1") && alts.length == 1) || (genotype.equals("1/2") && alts.length == 2))
                        {
                            if (hets.containsKey(chr))
                            {
                                ArrayList<Integer> al = hets.get(chr);
                                if (!al.contains(start))
                                {
                                    al.add(start);
                                    hets.put(chr, al);
                                }
                            }
                            else
                            {
                                ArrayList<Integer> al = new ArrayList();
                                al.add(start);
                                hets.put(chr, al);
                            }
                        }
                    }
                }
            }
        }
        for (Map.Entry pairs : hets.entrySet())
        {
            ArrayList<Integer> ints = (ArrayList<Integer>) pairs.getValue();
            heterozygousSites += ints.size();
        }
        return heterozygousSites;
    }

}
