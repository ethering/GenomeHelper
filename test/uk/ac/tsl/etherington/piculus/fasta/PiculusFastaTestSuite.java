/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fasta;

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
    uk.ac.tsl.etherington.piculus.fasta.RandomFastaTest.class, uk.ac.tsl.etherington.piculus.fasta.FastaTranslatorTest.class, uk.ac.tsl.etherington.piculus.fasta.FastaFeaturesTest.class, uk.ac.tsl.etherington.piculus.fasta.FastaSubstringsTest.class, uk.ac.tsl.etherington.piculus.fasta.ContigPositionTest.class, uk.ac.tsl.etherington.piculus.fasta.FastaMotifFinderTest.class
})
public class PiculusFastaTestSuite
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