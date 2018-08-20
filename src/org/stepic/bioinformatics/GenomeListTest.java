package org.stepic.bioinformatics;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;

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

    @Test
    void motifEnumeration() {

        var list = new GenomeList();
        list.add(new FuzzyGenome("ATTTGGC"));
        list.add(new FuzzyGenome("TGCCTTA"));
        list.add(new FuzzyGenome("CGGTATC"));
        list.add(new FuzzyGenome("GAAAATT"));
        final var result = list.motifEnumeration(3,1);
        assertTrue(result.size() == 4 &&
                Arrays
                        .asList("ATA ATT GTT TTT".split(" "))
                        .stream()
                        .allMatch(str -> result.contains(new Sequence(str)))
        );

        list = new GenomeList();
        list.add(new FuzzyGenome("ACGT"));
        list.add(new FuzzyGenome("ACGT"));
        list.add(new FuzzyGenome("ACGT"));
        final var result2 = list.motifEnumeration(3,0);
        assertTrue(result2.size() == 2 &&
                        Arrays
                                .asList("ACG CGT".split(" "))
                                .stream()
                                .allMatch(str -> result2.contains(new Sequence(str)))
                ,"This dataset checks for off\u00ADby\u00ADone errors, both at the beginning and at the end. The\n" +
                        "3\u00ADmers “ACG” and “CGT” both appear perfectly in all 3 strings in Dna. Thus, if your output\n" +
                        "doesn’t contain “ACG”, you are most likely not counting the first k\u00ADmer of every string.\n" +
                        "Similarly, if your output doesn’t contain “CGT”, you are most likely not counting the last k\u00ADmer\n" +
                        "of every string"
        );

    }

    @Test
    void greedyMotifSearch() {
        var list = new GenomeList(
                "GGCGTTCAGGCA\n" +
                        "AAGAATCAGTCA\n" +
                        "CAAGGAGTTCGC\n" +
                        "CACGTCAATCAC\n" +
                        "CAATAATATTCG".split("\n")
        );
        var result = list.greedyMotifSearch(3).stream().map(FuzzyGenome::getText).sorted().collect(Collectors.toList());
        var expected = Arrays.asList("CAG\n" +
                "CAG\n" +
                "CAA\n" +
                "CAA\n" +
                "CAA".split("\n")).stream().sorted().collect(Collectors.toList());
        assertIterableEquals(expected,result);

    }

//    @Test
//    void score() {
//        var list = new GenomeList(
//                "TCGGGGGTTTTT",
//                "CCGGTGACTTAC",
//                "ACGGGGATTTTC",
//                "TTGGGGACTTTT",
//                "AAGGGGACTTCC",
//                "TTGGGGACTTCC",
//                "TCGGGGATTCAT",
//                "TCGGGGATTCCT",
//                "TAGGGGAACTAC",
//                "TCGGGTATAACC"
//        );
//        assertEquals(30,list.score());
//    }
}