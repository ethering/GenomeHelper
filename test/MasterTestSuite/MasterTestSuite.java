/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package MasterTestSuite;

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
    uk.ac.tsl.etherington.genomehelper.bam.MappedSamRecordsTest.class,
    uk.ac.tsl.etherington.genomehelper.fasta.RandomFastaTest.class,
    uk.ac.tsl.etherington.genomehelper.fasta.FastaTranslatorTest.class,
    uk.ac.tsl.etherington.genomehelper.fasta.FastaFeaturesTest.class,
    uk.ac.tsl.etherington.genomehelper.fasta.FastaSubstringsTest.class,
    uk.ac.tsl.etherington.genomehelper.fasta.ContigPositionTest.class,
    uk.ac.tsl.etherington.genomehelper.fasta.FastaMotifFinderTest.class,
    uk.ac.tsl.etherington.genomehelper.fastq.FastqCompressionTest.class,
    uk.ac.tsl.etherington.genomehelper.fastq.FastqInterlacerTest.class,
    uk.ac.tsl.etherington.genomehelper.fastq.FastqQCTest.class,
    uk.ac.tsl.etherington.genomehelper.fastq.FastqMotifFinderTest.class,
    uk.ac.tsl.etherington.genomehelper.fastq.FastqParserTest.class,
    uk.ac.tsl.etherington.genomehelper.fastq.FastqJoinerTest.class,
    uk.ac.tsl.etherington.genomehelper.gff.GFFFeatureStatsTest.class
})
public class MasterTestSuite
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