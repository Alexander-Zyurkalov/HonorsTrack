package org.stepic.bioinformatics;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class GenomeListTest {

    @Test
    void distanceBetweenPatternAndThis() {
        var list = new GenomeList(
                "TTACCTTAAC",
                "GATATCTGTC",
                "ACGGCGTTCG",
                "CCCTAAAGAG",
                "CGTCAGAGGT"
        );
        assertEquals(5,list.distanceBetweenPatternAndThis(new Sequence("AAA")));

        list = new GenomeList(
                "TTTATTT CCTACAC GGTAGAG".split(" ")
        );
        assertEquals(3,list.distanceBetweenPatternAndThis(new Sequence("TAA")));

        list = new GenomeList(
                "AAACT AAAC AAAG".split(" ")
        );
        assertEquals(0,list.distanceBetweenPatternAndThis(new Sequence("AAA")));

        list = new GenomeList(
                "TTTTAAA CCCCAAA GGGGAAA".split(" ")
        );
        assertEquals(0,list.distanceBetweenPatternAndThis(new Sequence("AAA")));

        list = new GenomeList(
                "AAATTTT AAACCCC AAAGGGG".split(" ")
        );
        assertEquals(0,list.distanceBetweenPatternAndThis(new Sequence("AAA")));

    }

    @Test
    void medianString() {
        var list = new GenomeList(
                "AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTACGGGACAG".split(" ")
        );
        var result = list.medianString(3);
        assertTrue(
                result.equals(new Sequence("ACG")) ||
                result.equals(new Sequence("GAC"))
        );

        list = new GenomeList(
                "ACGT ACGT ACGT".split(" ")
        );
        result = list.medianString(3);
        assertTrue(
                result.equals(new Sequence("ACG")) ||
                result.equals(new Sequence("CGT"))
        );

        list = new GenomeList(
                "ATA ACA AGA AAT AAC".split(" ")
        );
        result = list.medianString(3);
        assertEquals(new Sequence("AAA"), result);


        list = new GenomeList(
                "AAG AAT".split(" ")
        );
        result = list.medianString(3);
        assertTrue(
                result.equals(new Sequence("AAG")) ||
                        result.equals(new Sequence("AAT"))
        );

    }
}