package org.stepic.bioinformatics;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class GenomeTest {
    @Test
    void patternCount() {

        assertEquals(2, new Genome("GCGCGC").patternCount("GCG")  );
        assertEquals(3, new Genome("GCGCGC").patternCount("GC")  );
        assertEquals(2, new Genome("GGACTTACTGACGTACG").patternCount("ACT")  );
        assertEquals(3, new Genome("ACGTACGTACGT").patternCount("CG")  );
        assertEquals(5, new Genome("ATCCGATCCCATGCCCATG").patternCount("CC")  );
        assertEquals(
                4,
                new Genome("AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAATAATTACAGAGTACACAACATCCAT")
                .patternCount("AAA"));
        assertEquals(4,
                new Genome( "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT")
                .patternCount("TTT"));
        assertEquals(
                9,
                new Genome("CTGTTTTTGATCCATGATATGTTATCTCTCCGTCATCAGAAGAACAGTGACGGATCGCCCTCTCTCTTGGTCAGGCGACCGTTTGCCATAATGCCCATGCTTTCCAGCCAGCTCTCAAACTCCGGTGACTCGCGCAGGTTGAGTA")
                .patternCount("CTC"));

    }

    @Test
    void frequentWords() {

        var genome = new Genome("ACGTTGCATGTCGCATGATGCATGAGAGCT");
        var output = genome.frequentWords(4);
        var expected = new ArrayList<Sequence>();
        expected.add(new Sequence("GCAT"));
        expected.add(new Sequence("CATG"));
        assertIterableEquals(expected,output);

        genome = new Genome("TGGTAGCGACGTTGGTCCCGCCGCTTGAGAATCTGGATGAACATAAGCTCCCACTTGGCTTATTCAGAGAACTGGTCAACACTTGTCTCTCCCAGCCAGGTCTGACCACCGGGCAACTTTTAGAGCACTATCGTGGTACAAATAATGCTGCCAC");
        output = genome.frequentWords(3);
        expected = new ArrayList<>();
        expected.add(new Sequence("TGG"));
        assertIterableEquals(expected,output);

        genome = new Genome("CAGTGGCAGATGACATTTTGCTGGTCGACTGGTTACAACAACGCCTGGGGCTTTTGAGCAACGAGACTTTTCAATGTTGCACCGTTTGCTGCATGATATTGAAAACAATATCACCAAATAAATAACGCCTTAGTAAGTAGCTTTT");
        output = genome.frequentWords(4);
        expected = new ArrayList<>();
        expected.add(new Sequence("TTTT"));
        assertIterableEquals(expected,output);

        genome = new Genome("ATACAATTACAGTCTGGAACCGGATGAACTGGCCGCAGGTTAACAACAGAGTTGCCAGGCACTGCCGCTGACCAGCAACAACAACAATGACTTTGACGCGAAGGGGATGGCATGAGCGAACTGATCGTCAGCCGTCAGCAACGAGTATTGTTGCTGACCCTTAACAATCCCGCCGCACGTAATGCGCTAACTAATGCCCTGCTG");
        output = genome.frequentWords(5);
        expected = new ArrayList<>();
        expected.add(new Sequence("AACAA"));
        assertIterableEquals(expected,output);

        genome = new Genome("CCAGCGGGGGTTGATGCTCTGGGGGTCACAAGATTGCATTTTTATGGGGTTGCAAAAATGTTTTTTACGGCAGATTCATTTAAAATGCCCACTGGCTGGAGACATAGCCCGGATGCGCGTCTTTTACAACGTATTGCGGGGTAAAATCGTAGATGTTTTAAAATAGGCGTAAC");
        output = genome.frequentWords(5);
        expected = new ArrayList<>();
        expected.add(new Sequence("TTTTA"));
        expected.add(new Sequence("AAAAT"));
        expected.add(new Sequence("GGGGT"));

        assertIterableEquals(expected,output);




    }

    @Test
    void findAllPositionOfThe() {
        var seq = new Genome("GATATATGCATATACTT");
        var found = seq.findAllPositionsOfThe(new Sequence("ATAT"));
        ArrayList<Integer> expected = new ArrayList<>();
        expected.add(1); expected.add(3); expected.add(9);
        assertIterableEquals(expected,found);
    }

    @Test
    void clumpFinding() {
//        ======================================================
        var expected = new HashSet<Sequence>();
        expected.add(new Sequence("CGACA"));
        expected.add(new Sequence("GAAGA"));
        var genome = new Genome("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA");
//        assertIterableEquals(expected,genome.clumpFinding(5,50,4));
//        ======================================================
        expected = new HashSet<>();
        expected.add(new Sequence("A"));
        expected.add(new Sequence("C"));
        expected.add(new Sequence("G"));
        expected.add(new Sequence("T"));
        genome = new Genome("ACGTACGT");
        assertIterableEquals(expected,genome.clumpFinding(1,5,2));
//        ======================================================
        expected = new HashSet<>();
        genome = new Genome("CCACGCGGTGTACGCTGCAAAAAGCCTTGCTGAATCAAATAAGGTTCCAGCACATCCTCAATGG" +
                "TTTCACGTTCTTCGCCAATGGCTGCCGCCAGGTTATCCAGACCTACAGGTCCACCAAAGAACTT" +
                "ATCGATTACCGCCAGCAACAATTTGCGGTCCATATAATCGAAACCTTCAGCATCGACATTCAAC" +
                "ATATCCAGCG");
        Set<Sequence> clumps = genome.clumpFinding(3,25,3);

        assertTrue(clumps.contains(new Sequence("AAA")));
        assertTrue(clumps.contains(new Sequence("CAG")));
        assertTrue(clumps.contains(new Sequence("CAT")));
        assertTrue(clumps.contains(new Sequence("CCA")));
        assertTrue(clumps.contains(new Sequence("GCC")));
        assertTrue(clumps.contains(new Sequence("TTC")));
        assertTrue(clumps.size() == 6);

    }

    @Test
    void findAllPositionsOfTheApprPattern() {
        var genome = new Genome("CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT");
        assertIterableEquals(Arrays.asList(6, 7, 26, 27),
                genome.findAllPositionsOfTheApprPattern(new Sequence("ATTCTGGA"), 3));

        genome = new Genome("TTTTTTAAATTTTAAATTTTTT");
        assertIterableEquals(Arrays.asList(4,5,6, 7, 8, 11, 12, 13, 14, 15),
                genome.findAllPositionsOfTheApprPattern(new Sequence("AAA"), 2),
                "This dataset checks if you are only counting instances where the number of mismatches is\n" +
                        "exactly equal to d (i.e. ignoring instances where mismatch < d)");

        genome = new Genome("GAGCGCTGGGTTAACTCGCTACTTCCCGACGAGCGCTGTGGCGCAAATTGGCGATGAAACTGCAGAGAGAACTGGTCATCCAACTGAATTCTCCCCGCTATCGCATTTTGATGCGCGCCGCGTCGATT");
        assertIterableEquals(Arrays.asList(0, 30, 66),
                genome.findAllPositionsOfTheApprPattern(new Sequence("GAGCGCTGG"), 2),
                "This dataset checks if your code has an off-by-one error at the beginning of Text (i.e. your\n" +
                        "code is not checking the the leftmost substring of Text).");

    }




}