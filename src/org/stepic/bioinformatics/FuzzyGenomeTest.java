package org.stepic.bioinformatics;

import org.junit.jupiter.api.Test;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

class FuzzyGenomeTest {

    /*@Test
    void findAllPositionsOfTheApprPattern() {
        var genome = new FuzzyGenome("CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT");
        assertIterableEquals(Arrays.asList(6, 7, 26, 27),
                genome.findAllPositionsOfTheApproximatePattern(new Sequence("ATTCTGGA"), 3));

        genome = new FuzzyGenome("TTTTTTAAATTTTAAATTTTTT");
        assertIterableEquals(Arrays.asList(4,5,6, 7, 8, 11, 12, 13, 14, 15),
                genome.findAllPositionsOfTheApproximatePattern(new Sequence("AAA"), 2),
                "This dataset checks if you are only counting instances where the number of mismatches is\n" +
                        "exactly equal to d (i.e. ignoring instances where mismatch < d)");

        genome = new FuzzyGenome("GAGCGCTGGGTTAACTCGCTACTTCCCGACGAGCGCTGTGGCGCAAATTGGCGATGAAACTGCAGAGAGAACTGGTCATCCAACTGAATTCTCCCCGCTATCGCATTTTGATGCGCGCCGCGTCGATT");
        assertIterableEquals(Arrays.asList(0, 30, 66),
                genome.findAllPositionsOfTheApproximatePattern(new Sequence("GAGCGCTGG"), 2),
                "This dataset checks if your code has an off-by-one error at the beginning of Text (i.e. your\n" +
                        "code is not checking the the leftmost substring of Text).");

    }*/

    @Test
    void approximatePatternCount() {
        var genome = new FuzzyGenome("TTTAGAGCCTTCAGAGG");
        assertEquals(4,
                genome.approximatePatternCount(new Sequence("GAGG"),2));

        genome = new FuzzyGenome("AAA");
        assertEquals(2,
                genome.approximatePatternCount(new Sequence("AA"),0),
                "This dataset first checks that your code is returning a value larger than 0 (returning a\n" +
                        "value of 0 would imply that your code is probably flipping the two lines). Then, it checks that\n" +
                        "your code is correctly handling overlapping occurrences (i.e. returning 2 instead of 1).");

    }


    @Test
    void fastFrequentWordsWithFreqFunctionally() {
        var genome = new FuzzyGenome("ACGTTGCATGTCGCATGATGCATGAGAGCT");
        ArrayList<Sequence> output = new ArrayList<>(
                genome.frequentWordsWithMismatch(4, 1));
        var expected = new ArrayList<Sequence>();
        assertEquals(2,output.size());
        assertTrue(output.contains(new Sequence("ACAT")));
        assertTrue(output.contains(new Sequence("ATGT")));


        genome = new FuzzyGenome("AAAAAAAAAA");
        output = new ArrayList<>(
                genome.frequentWordsWithMismatch(2, 1));

        assertEquals(2,output.size());
        assertTrue(output.contains(new Sequence("AT")));
        assertTrue(output.contains(new Sequence("TA")),
                "This dataset checks that your code includes k-mers that do not actually appear in Text.\n" +
                        "Notice here that, although AT nor TA actually appear in Text, they are valid because they appear\n" +
                        "in Text with up to 1 mismatch (i.e. 0 or 1 mismatch).");



        genome = new FuzzyGenome("AGTCAGTC");
        output = new ArrayList<>(
                genome.frequentWordsWithMismatch(4, 2));
        assertEquals(2,output.size());
        assertTrue(output.contains(new Sequence("AATT")));
        assertTrue(output.contains(new Sequence("GGCC")),
                "This dataset makes sure that your code is not accidentally swapping k and d");


        genome = new FuzzyGenome("AAT");
        output = new ArrayList<>(
                genome.frequentWordsWithMismatch(3, 0));
        assertEquals(2,output.size());
        assertTrue(output.contains(new Sequence("AAT")));
        assertTrue(output.contains(new Sequence("ATT")),
                "This dataset checks that your code is looking at BOTH Text and its Reverse Complement\n" +
                        "(i.e. not just looking at Text, and not just looking at the Reverse Complement of Text, but looking\n" +
                        "at both).");




        genome = new FuzzyGenome("AATTAATTGGTAGGTAGGTA");
        output = new ArrayList<>(
                genome.frequentWordsWithMismatch(4, 0));
        assertEquals(1,output.size());
        assertTrue(output.contains(new Sequence("AATT")),
                "This dataset makes sure you are finding k-mers in both Text and the Reverse Complement\n" +
                        "of Text.");


        genome = new FuzzyGenome("ATA");
        final ArrayList<Sequence> output2 = new ArrayList<>(
                genome.frequentWordsWithMismatch(3, 1));
        assertEquals(20,output2.size());
        assertTrue(
            Arrays
                .asList("AAA AAT ACA AGA ATA ATC ATG ATT CAT CTA GAT GTA TAA TAC TAG TAT TCT TGT TTA TTT".split(" "))
                .stream()
                .allMatch(str -> output2.contains(new Sequence(str)))
        );


        genome = new FuzzyGenome("TAGCG");
        final ArrayList<Sequence> output3 = new ArrayList<>(
                genome.frequentWordsWithMismatch(2, 1));
        assertEquals(4,output3.size());
        assertTrue(
            Arrays
                .asList("CA CC GG TG".split(" "))
                .stream()
                .allMatch(str -> output3.contains(new Sequence(str)))
        );


    }
}