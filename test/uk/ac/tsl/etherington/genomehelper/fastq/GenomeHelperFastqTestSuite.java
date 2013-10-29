/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fastq;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.runner.RunWith;
import org.junit.runners.Suite;

/**
 *
 * @author ethering
 */
@RunWith(Suite.class)
@Suite.SuiteClasses(
{
    uk.ac.tsl.etherington.genomehelper.fastq.FastqCompressionTest.class, uk.ac.tsl.etherington.genomehelper.fastq.FastqInterlacerTest.class, uk.ac.tsl.etherington.genomehelper.fastq.FastqQCTest.class, uk.ac.tsl.etherington.genomehelper.fastq.FastqMotifFinderTest.class, uk.ac.tsl.etherington.genomehelper.fastq.FastqParserTest.class, uk.ac.tsl.etherington.genomehelper.fastq.FastqJoinerTest.class
})
public class GenomeHelperFastqTestSuite
{

    @BeforeClass
    public static void setUpClass() throws Exception
    {
    }

    @AfterClass
    public static void tearDownClass() throws Exception
    {
    }
    
}