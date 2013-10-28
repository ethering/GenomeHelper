/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fastq;

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
    uk.ac.tsl.etherington.piculus.fastq.FastqCompressionTest.class, uk.ac.tsl.etherington.piculus.fastq.FastqInterlacerTest.class, uk.ac.tsl.etherington.piculus.fastq.FastqQCTest.class, uk.ac.tsl.etherington.piculus.fastq.FastqMotifFinderTest.class, uk.ac.tsl.etherington.piculus.fastq.FastqParserTest.class, uk.ac.tsl.etherington.piculus.fastq.FastqJoinerTest.class
})
public class PiculusFastqTestSuite
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