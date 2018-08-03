package org.stepic.bioinformatics;

import static org.junit.jupiter.api.Assertions.*;


class PatternCountTest {

    @org.junit.jupiter.api.Test
    void patternCount() {
        assertEquals(2, PatternCount.patternCount("GCGCG","GCG")  );
        assertEquals(2, PatternCount.patternCount("GGACTTACTGACGTACG","ACT")  );
        assertEquals(3, PatternCount.patternCount("ACGTACGTACGT","CG")  );
        assertEquals(5, PatternCount.patternCount("ATCCGATCCCATGCCCATG","CC")  );
        assertEquals(4,PatternCount.patternCount(
                "AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAATTTTATTGAC" +
                "TTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAATAATTACAGAGTA" +
                "CACAACATCCAT","AAA"));
        assertEquals(4,PatternCount.patternCount(
                "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTG" +
                     "CATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGC" +
                     "CGACTTT","TTT"));
        assertEquals(9,PatternCount.patternCount(
                "CTGTTTTTGATCCATGATATGTTATCTCTCCGTCATCAGAAGAACAGTGACGGATCGCCCTCTC" +
                        "TCTTGGTCAGGCGACCGTTTGCCATAATGCCCATGCTTTCCAGCCAGCTCTCAAACTCCGGTGA" +
                        "CTCGCGCAGGTTGAGTA","CTC"));

    }
}