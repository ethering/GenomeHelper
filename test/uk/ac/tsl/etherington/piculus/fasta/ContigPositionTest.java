/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fasta;

import java.util.ArrayList;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author ethering
 */
public class ContigPositionTest
{
    
     ContigPosition cp1 = new ContigPosition("cp1", 1, 10, '+');
     ContigPosition cp2 = new ContigPosition("cp2", 1, 10);
     ContigPosition cp3 = new ContigPosition("cp2", 10, 20);
     ContigPosition cp4 = new ContigPosition("cp4", 21, 30);
     ContigPosition cp5 = new ContigPosition("cp4", 21, 30);
     ContigPosition cp6 = new ContigPosition("cp6", 21, 30, '-');
     ContigPosition cp7 = new ContigPosition("cp2", 1, 10);
     

    public ContigPositionTest()
    {
    }

    @BeforeClass
    public static void setUpClass()
    {
    }

    @AfterClass
    public static void tearDownClass()
    {
    }

    /**
     * Test of validateNewStrand method, of class ContigPosition.
     */
    @Test
    public void testValidateNewStrand()
    {
        System.out.println("validateNewStrand");
        char strand = 'x';
        //try setting an invalid strand
        boolean result = cp2.validateNewStrand(strand);
        assertEquals(result, false);
        //now try setting a valid strand
        strand = '-';
        result = cp2.validateNewStrand(strand);
        assertEquals(result, true);

    }

    /**
     * Test of validateCurrentStrand method, of class ContigPosition.
     */
    @Test
    public void testValidateCurrentStrand()
    {
        System.out.println("validateStrand");
        //strand for cp1 has been set manually, so should be good (true)
        boolean result = cp1.validateCurrentStrand();
        assertEquals(result, true);
        //strand for cp2 has not be set, so should be '+' by default
        result = cp2.validateCurrentStrand();
        //make sure it's valid
        assertEquals(result, true);
        //and then make sure it's a '+'
        assertEquals(cp2.getStrand(), '+');
    }

    /**
     * Test of setStrand method, of class ContigPosition.
     */
    @Test
    public void testSetStrand()
    {
        System.out.println("setStrand");
        char strand = '-';
        cp2.setStrand(strand);
        assertEquals(cp2.getStrand(), '-');
    }

    /**
     * Test of getStrand method, of class ContigPosition.
     */
    @Test
    public void testGetStrand()
    {
        System.out.println("getStrand");
        char result = cp2.getStrand();
        assertEquals(result, '+');
        result = cp6.getStrand();
        assertEquals(result, '-');
    }

    /**
     * Test of getContigid method, of class ContigPosition.
     */
    @Test
    public void testGetContigid()
    {
        System.out.println("getContigid");
        String result = cp1.getContigid();
        assertEquals(result, "cp1");

    }

    /**
     * Test of getStart method, of class ContigPosition.
     */
    @Test
    public void testGetStart()
    {
        System.out.println("getStart");
        int result = cp1.getStart();
        assertEquals(result, 1);
    }

    /**
     * Test of getEnd method, of class ContigPosition.
     */
    @Test
    public void testGetEnd()
    {
        System.out.println("getEnd");
        int result = cp1.getEnd();
        assertEquals(result, 10);
    }

    /**
     * Test of setContigid method, of class ContigPosition.
     */
    @Test
    public void testSetContigid()
    {
        System.out.println("setContigid");
        String contigid = "cp22";
        cp1.setContigid(contigid);
        assertEquals(cp1.contigid, "cp22");
    }

    /**
     * Test of setStart method, of class ContigPosition.
     */
    @Test
    public void testSetStart()
    {
        System.out.println("setStart");
        cp1.setStart(2);
        assertEquals(cp1.start, 2);
    }

    /**
     * Test of setEnd method, of class ContigPosition.
     */
    @Test
    public void testSetEnd()
    {
        System.out.println("setEnd");
        cp1.setEnd(20);
        assertEquals(cp1.end, 20);
    }

    /**
     * Test of equals method, of class ContigPosition.
     */
    @Test
    public void testEquals()
    {
        System.out.println("equals");
        boolean result = cp4.equals(cp5);
        assertEquals(result, true);
        result = cp4.equals(cp6);
        assertEquals(result, false);
    }

    /**
     * Test of overlaps method, of class ContigPosition.
     */
    @Test
    public void testOverlaps()
    {
        System.out.println("overlaps");
        boolean result = cp2.overlaps(cp3);
        assertEquals(result, true);
        result = cp3.overlaps(cp4);
        assertEquals(result, false);
        
        result = cp2.overlaps(cp7);
        assertEquals(result, true);
        
        result = cp7.overlaps(cp2);
        assertEquals(result, true);
    }
    /**
     * Test of exists method, of class ContigPosition.
     */
    @Test
    public void testExists()
    {
        System.out.println("exists");

        //add some CPs to an ArrayList
        ArrayList<ContigPosition> cpArray = new ArrayList<>();
        cpArray.add(cp2);
        cpArray.add(cp3);
        cpArray.add(cp4);
        //get a CP that exists and check that it exist
        boolean result = cp5.exists(cpArray);
        assertEquals(result, true);
        //change the end postion
        cp5.setEnd(31);
        System.out.println("cp5 " +cp5.toString());
        //and test to see if it's still found
        result = cp5.exists(cpArray);
        assertEquals(result, false);
    }

    /**
     * Test of toString method, of class ContigPosition.
     */
    @Test
    public void testToString()
    {
        System.out.println("toString");
        String result = cp2.toString();
        //get the class of the returned object (should be a java.lang.String)
        Class cl = result.getClass();
        //create a String object directly
        String strClass = "I'm a String";
        //and test that both are Strings
        assertEquals(strClass.getClass(), cl);
    }
}