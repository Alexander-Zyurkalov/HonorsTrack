package org.stepic.bioinformatics;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.*;

class GenomeListWithPseudocountsTest {


    @Test
    void greedyMotifSearchFunctional() {
        var list = new GenomeListWithPseudocounts(
                ("GGCGTTCAGGCA\n" +
                "AAGAATCAGTCA\n" +
                "CAAGGAGTTCGC\n" +
                "CACGTCAATCAC\n" +
                "CAATAATATTCG").split("\n")
        );
        var result = list.greedyMotifSearchFunctional(3).stream().map(FuzzyGenome::getText).sorted().collect(Collectors.toList());
        var expected = Arrays.asList(
                ("TTC\n" +
                "ATC\n" +
                "TTC\n" +
                "ATC\n" +
                "TTC").split("\n")).stream().sorted().collect(Collectors.toList());
        assertIterableEquals(expected,result);



        list = new GenomeListWithPseudocounts(
                ("AGGCGGCACATCATTATCGATAACGATTCGCCGCATTGCC\n" +
                        "ATCCGTCATCGAATAACTGACACCTGCTCTGGCACCGCTC\n" +
                        "AAGCGTCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT\n" +
                        "AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG\n" +
                        "AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT\n" +
                        "AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT\n" +
                        "AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG\n" +
                        "AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA").split("\n")
        );
        result = list.greedyMotifSearchFunctional(5).stream().map(FuzzyGenome::getText).sorted().collect(Collectors.toList());
        expected = Arrays.asList(
                ("AGGCG\n" +
                        "ATCCG\n" +
                        "AAGCG\n" +
                        "AGTCG\n" +
                        "AACCG\n" +
                        "AGGCG\n" +
                        "AGGCG\n" +
                        "AGGCG").split("\n")).stream().sorted().collect(Collectors.toList());
        assertIterableEquals(expected,result);

        list = new GenomeListWithPseudocounts(
                ("GCACATCATTAAACGATTCGCCGCATTGCCTCGATAGGCG\n" +
                        "TCATAACTGACACCTGCTCTGGCACCGCTCATCCGTCGAA\n" +
                        "AAGCGGGTATAGCCAGATAGTGCCAATAATTTCCTTCGGC\n" +
                        "AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG\n" +
                        "AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT\n" +
                        "AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT\n" +
                        "AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTGAGGCG\n" +
                        "AATTGAACATCTTACTCTTTTCGCTTTCAAAAAAAAGGCG\n").split("\n")
        );
        result = list.greedyMotifSearchFunctional(5).stream().map(FuzzyGenome::getText).sorted().collect(Collectors.toList());
        expected = Arrays.asList(
                ("AGGCG\n" +
                        "TGGCA\n" +
                        "AAGCG\n" +
                        "AGGCA\n" +
                        "CGGCA\n" +
                        "AGGCG\n" +
                        "AGGCG\n" +
                        "AGGCG").split("\n")).stream().sorted().collect(Collectors.toList());
        assertIterableEquals(expected,result);


        list = new GenomeListWithPseudocounts(
                ("GCACATCATTATCGATAACGATTCATTGCCAGGCGGCCGC\n" +
                        "TCATCGAATAACTGACACCTGCTCTGGCTCATCCGACCGC\n" +
                        "TCGGCGGTATAGCCAGATAGTGCCAATAATTTCCTAAGCG\n" +
                        "GTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTGAGTCG\n" +
                        "GACGGCAACTACGGTTACAACGCAGCAAGAATATTAACCG\n" +
                        "TCTGTTGTTGCTAACACCGTTAAGCGACGGCAACTAGGCG\n" +
                        "GCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTGAAGCG\n" +
                        "AAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAAAATTG").split("\n")
        );
        result = list.greedyMotifSearchFunctional(5).stream().map(FuzzyGenome::getText).sorted().collect(Collectors.toList());
        expected = Arrays.asList(
                ("GGCGG\n" +
                        "GGCTC\n" +
                        "GGCGG\n" +
                        "GGCAG\n" +
                        "GACGG\n" +
                        "GACGG\n" +
                        "GGCGC\n" +
                        "GGCGC").split("\n")).stream().sorted().collect(Collectors.toList());
        assertIterableEquals(expected,result);


        list = new GenomeListWithPseudocounts(
                ("ACGAGAACTGTAATGGAGACCAATCGGGTCGTATGTACGACACGGATCTTCTGTATCGATCATCGCTTAACTTATACGATCTCATTCTCACGACGATCCTCAACCCCGGATACCCGCACTGCCTCAATCCGAAGTACTGCGTAGTACTTTACCCTT\n" +
                        "ACCTGTGATCACTCAGAGAACAGAGGATCCGGGTTGGATGTCAGTGTTATGCCAAGAAACGAGACCTAAGGTGCCGTCCCCGGCGGAATGCTTCTCGCTTCCCCTTCTAAAGGGTCCTGGCAAAATGCTTGTGACTTTGAATGCCCTCACTGAACT\n" +
                        "AAATGTGAAACCTATATCAGTCATTATACCGGCGCGATGTTAACTCGCACCTGATTGCAAGGTCACTGATCGCGTCACTACACTAGAGTTTATTCATACCTGCATGGGGGGCATGATGGATGAATTTTAACTAGGTGATCGTGACCAATGTTCACC\n" +
                        "TTTAACCGAATGAGAAGGGTTTCGTTTTAGGCCGTCGATCCGCCGCTTCTTCTCTCGACTGAATCGGAGGTTTTAATGTGTGTGCTAAAGTGAACTGCCTAATAACCCGCTAGGGGGATGATTTATTGTATCTGCAATGTACGGTGGTACATCCCG\n" +
                        "GCTTTCCCGTGTTCTTCGTGCGCTCGGGACCAGAGATCAGAACTAGTGCTAAAGCCATGGAAGCATTTAGCCGGCCGATGTAGTAAAGTTGCCCATATTTCTCCCTAAACCAGCGTATACGTCGAAAACTTTCTTCACTCCACTATGTATAGATGG\n" +
                        "CGCTTAGCGCATCCTCGCTAAACTAATTTGTAGGAGCAATCGCATGTCGACACCGAGGAAGACAACAAGTACATTAGTTACACGCTTCTTACTGGGGATAAAAATCTAGGATCGCGTGATCGGCGTGCTCCGGCGTCAGGTGACTTGATGGCCCCA\n" +
                        "AGTCATGACAATTCGCACGAGGTGACTTTCAATCGACTTACGTACTGCGCGTCGATGTCGCTTCACTCCACTAATACCCATATTCCTAACACGCCAGTACGGTTCAGAGTCGGCGGCTGAGGGGCCCTAGAAACGAGACACCCTAGAGGCTCTGGG\n" +
                        "ATAATAAACGAGATTTAGTGCTCCAGGTCCCTCTAACTTCGCTAGACTGTTCACGGTACTTAGAGTAGCTCAGAAATCGCCTGTCTTCGGGTTGGTTTTTTCAGGAGGTGCTCTGTGCGTTAACATACCAAGCCATAGCTGCTTTTCCTCTACAAA\n" +
                        "CCTACCGGAGGACTTCACTGAACTAAGTACTGGGGTGCTAAGGTCGAACGAGATGACTGGACCCTACTCTCTACGGGACCGACGCCCCAGGGCTTAATTCATATTGACTAGATTTATGATAATAATAGACTCGGCGGTTTGTAGCTTTCCCCTGAA\n" +
                        "CCGTTACTCCGCCGCTGCTGATCCTACTCCACGTGGGGCCCCCCAATATGCACATCTTATCGTCCCTGGACTGGCAGTTGATGCGAAATAATATTCGTGGGTATGATAACGCGCTATACTGATAGAACCACGGGGACTCCTGTATTCGTCTCGCCA\n" +
                        "ATTCGTGGGTTGAAGCCTTTAAACGGGATGGCCAAGTTGATTGGGTCTAATTGATATTAATTCTGGTGTACTGTACGAGACCGGTGCAGCACGGACGGGCGGTTTCAACAATACTCGTGCGCCTGGAGGGTACCTCGCTGAACTCGCATATCAGGT\n" +
                        "TGGTTACCCTCTCTTCACTTAACTCTTGTATCAAGACGTTTCTGTGAGACAAAGCAATGGCCGGACTTTGGGGCGCGTGCTCTGGAGATCCAAGCACTCGAGGTCAGCGGTATAATATAACGCATCCAACATGCAGACTGTGCGTGGGGGCCCAAA\n" +
                        "CCTTAGCGGTGTCGGCCATCATTTATCGAGCTGGAACTCCGGTGGGTAACGTGGATCCGCAGCAGGCCTTACCGTCACTTAACTTGTTACTAGAACATACGTGGAACCTATGCATGATCGAGATAGAGTGCGCTCCGCGGTACACCGCGCTCTATA\n" +
                        "AGTTACAACACGACAGGAACTATTTGCTAGGCGTTACATCTCTTTACTTTTAGCCCCTGGATTTTGAACGCATGTCAACACGTTCCACCATGGGTATAAGAATGCATGGACAGGGTTAATGAATGTGTCTCGGTCCGTTAGCTGTTACAATATACC\n" +
                        "CTTAGCAGCCCACGTCGCTGGACTGTGTTATATTACGAGGTCGAATAGTGAGGTTAACAGTCTCCGTTGTAACTTAATCCCGATATCACGCAGTGTATATGGTCGCTGTAGCTTTCTGGGCAGCTCGCACACCGCCAATTCGCAGAGGCGACCAGA\n" +
                        "TGGTAATAGGCTTCGAAATAACTCTTGGATTGCAACGAAGGTCCGAGCCTTCTCTGCACTTGGATACATTTTGGACATATGAAAGGATGGGTGCTCGGGATGGGACTTTGGTTGCCTGCAAGACGGCGAGACCACCTTACTGAAACCAACATCTTA\n" +
                        "TTACCGGTAATCTCTGATCGCCCATGCCGTCAGGTGCCTTAATTTAAGCGAGAGCTAAATAGAAAGCTGCGCGGTTTTAGCAAAATGAAGTATCAGGAATAACATGGGTTAATGTCACAAACCTGAGGGTTACCTCTCTGCACTCCAGGTCCAGCA\n" +
                        "AATTCAGCTGGTCAGTCACCACAACGTCTCTAGACTCGCACACCCAACTATTATATCACGTACAAGCCGCCCCACAACCGGCATGATAATGTCTGCACGGCCCAACTAACACGCCAGATGACGTACTTTTCGCGGCAGAGCAGTATTCGAACTCAA\n" +
                        "CGTACCCGTTATACAAGCACCATCTACAAAACGTTAGGTGTCAACGATCGTGGGCCGGACTTAGGGGTGAGACCTTAAAGCACACATTCCCTCCCACATCACTTCACTTGTCAAAAATAAAGTCGAATGATGACTACCTCAATTCCTCGCGAAAGC\n" +
                        "GCTGACACGTATTACGAACCGAGACCAGGCGCGACCCCAACCTGGACTAACAGCTATCCCTTTGTTACTTGGCACGTACGGTCAAATCCCGGTGGAGTTATTTACCGAGGGTGCGCCATGCATCGCTCAACTCCTAATGGTCGCCTGTACTTGGTC\n" +
                        "CGGCCTTGGTCCAATTCATCGTAACATTCTGTGAATCTACGGGTACTACTCGACCACGCTTGCTTCGCTGAGCTGTACCGCAAAATCGAATGGACCCAATAATCTGAATCCTTCGGTATACATCACTAGACTCACGATTCAGTATGCCCTCAATCA\n" +
                        "ATCCGAAACATTCAGTTCCGATGGCAACGACGACCACCCAGCACGCATCATCACTCGACTGACCCTGCTCGAATACAGGCGTATCTAACAACCAGCGGATCCAGGGCCCATGCTGAGGCTATTGTCACTCCCGCCCACCCTGTATGTATTCGGATT\n" +
                        "CCCGGATCGCCCGGTCAAGCACTCCAGGGCTTCAGCGCTGTGTACCTCTCTGCACACTCCGGTGGTGGGCCCCGTCCCTACACTGCGTATTCAAGTGATAGCTCGTACGATCCGCATCGAAGGCTGACTGCCCCTCTACAACAGTGCGCTCGCACT\n" +
                        "CCCGTATGTAGGTGATTAGACCCACCAAGCAATCCGCATGTTGCGTCACCGTAGATATATGCAGCGGTCATTCTTCGCTTGACTCTCTGACTTGCCCATTTAGAAACTATCCACTTAATTCCCCTTGGACTTGGGGCCTGTTTTCGGTCGCTTTGT\n" +
                        "CATTTCTATAAAGCTACAATAATAATCCGCGCTGTCGGCAGACGTGGTACCGACCCTACTCCTACCGTTTGAGAGATGGAGGGTCTTCCCTGAACTAACGGCATGCATGAGAGGGGTACGACCCTGGTACTTCTGAAACCAGCATCCGCGGCGACG").split("\n")
        );
        result = list.greedyMotifSearchFunctional(12).stream().map(FuzzyGenome::getText).sorted().collect(Collectors.toList());
        expected = Arrays.asList(
                ("CATCGCTTAACT\n" +
                        "CCTCACTGAACT\n" +
                        "CGTCACTACACT\n" +
                        "CTTCTCTCGACT\n" +
                        "CTTCACTCCACT\n" +
                        "CCTCGCTAAACT\n" +
                        "CTTCACTCCACT\n" +
                        "CTTCGCTAGACT\n" +
                        "CTTCACTGAACT\n" +
                        "CGTCCCTGGACT\n" +
                        "CCTCGCTGAACT\n" +
                        "CTTCACTTAACT\n" +
                        "CGTCACTTAACT\n" +
                        "CATCTCTTTACT\n" +
                        "CGTCGCTGGACT\n" +
                        "CTTCTCTGCACT\n" +
                        "CCTCTCTGCACT\n" +
                        "CGTCTCTAGACT\n" +
                        "CATCACTTCACT\n" +
                        "CATCGCTCAACT\n" +
                        "CATCACTAGACT\n" +
                        "CATCACTCGACT\n" +
                        "CGTCCCTACACT\n" +
                        "CTTCGCTTGACT\n" +
                        "CTTCCCTGAACT").split("\n")).stream().sorted().collect(Collectors.toList());
        assertIterableEquals(expected,result);
    }

    @Test
    void getMotifsFromProfile() {
        var profile = new ProfileHash();
        profile.put('A', Arrays.asList(4.0/5, 0.0, 0.0, 0.1 ));
        profile.put('C', Arrays.asList(0.0, 3.0/5, 1.0/5, 0.0));
        profile.put('G', Arrays.asList(1.0/5, 1.0/5, 4.0/5.0, 0.0));
        profile.put('T', Arrays.asList(0.0, 1.0/5, 0.0, 4.0/5));
        var profileList = profile.toProfileList();

        var genomeList = new GenomeListWithPseudocounts(
                "TTACCTTAAC",
                "GATGTCTGTC",
                "ACGGCGTTAG",
                "CCCTAACGAG",
                "CGTCAGAGGT"
        );
        var motifs = genomeList.getMotifsFromProfile(profileList);
        var extected = new GenomeListWithPseudocounts(
                "ACCT",
                "ATGT",
                "GCGT",
                "ACGA",
                "AGGT"
        );
        assertIterableEquals(extected,motifs);

    }

    @Test
    void randomizedMotifSearch() {
        var list = new GenomeListWithPseudocounts(
                ("TTACCTTAAC\n" +
                "GATGTCTGTC\n" +
                "ACGGCGTTAG\n" +
                "CCCTAACGAG\n" +
                "CGTCAGAGGT").split("\n")
        );
        var expectedMotifs = new GenomeListWithPseudocounts(
                ("ACCT\n" +
                        "ATGT\n" +
                        "GCGT\n" +
                        "ACGA\n" +
                        "AGGT").split("\n"));

        var randomized = list.randomizedMotifSearch(4,10);
        assertEquals(expectedMotifs.score(),randomized.score(),3);


    }

    @Test
    void getRandomMotifs() {

        var list = new GenomeListWithPseudocounts("ATGCCGTA" );
        List<String> listofList =
            Stream.iterate(
                    list.getRandomMotifs(3),
                    l -> list.getRandomMotifs(3))
                .limit(24)
                .map(l->l.get(0).getText())
                .distinct()
                .sorted()
                .collect(Collectors.toList()
            );
        var expected = Arrays.asList(
                ("ATG\n" +
                        "TGC\n" +
                        "GCC\n" +
                        "CCG\n" +
                        "CGT\n"+
                        "GTA"
                ).split("\n")).stream().sorted().collect(Collectors.toList());
        assertIterableEquals(expected,listofList);
    }

    @Test
    void gibblsSampler() {
        var list = new GenomeListWithPseudocounts(
                ("CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA\n" +
                "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG\n" +
                "TAGTACCGAGACCGAAAGAAGTATACAGGCGT\n" +
                "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC\n" +
                "AATCCACCAGCTCCACGTGCAATGTTGGCCTA").split("\n")
        );
        var expectedMotifs = new GenomeListWithPseudocounts(
                ("TCTCGGGG\n" +
                "CCAAGGTG\n" +
                "TACAGGCG\n" +
                "TTCAGGTG\n" +
                "TCCACGTG").split("\n"));

        var randomized = list.gibblsSampler(8,2000 ,583);
        System.out.println(expectedMotifs.score());
        System.out.println(randomized.score());
//        assertEquals(expectedMotifs.score(),randomized.score(),3);
        assertEquals(expectedMotifs,randomized);
    }

    @Test
    void getMotifWithoutStr() {
        var motifs = new GenomeListWithPseudocounts(
                ("TCTCGGGG\n" +
                 "CCAAGGTG\n" +
                 "TACAGGCG\n" +
                 "TTCAGGTG\n" +
                 "TCCACGTG").split("\n"));
        var result = motifs.getMotifWithoutStr(0);
        var expected = new GenomeListWithPseudocounts(
            (
            "CCAAGGTG\n" +
            "TACAGGCG\n" +
            "TTCAGGTG\n" +
            "TCCACGTG").split("\n"));

        assertIterableEquals(expected,result);


        result = motifs.getMotifWithoutStr(1);
        expected = new GenomeListWithPseudocounts(
                (
                "TCTCGGGG\n"  +
                "TACAGGCG\n" +
                "TTCAGGTG\n" +
                "TCCACGTG"
        ).split("\n"));
        assertIterableEquals(expected,result);


        result = motifs.getMotifWithoutStr(2);
        expected = new GenomeListWithPseudocounts(
                (
                "TCTCGGGG\n"  +
                "CCAAGGTG\n" +
                "TTCAGGTG\n" +
                "TCCACGTG"
        ).split("\n"));
        assertIterableEquals(expected,result);



        result = motifs.getMotifWithoutStr(3);
        expected = new GenomeListWithPseudocounts(
                ("TCTCGGGG\n" +
                "CCAAGGTG\n" +
                "TACAGGCG\n" +
                "TCCACGTG").split("\n"));
        assertIterableEquals(expected,result);

        result = motifs.getMotifWithoutStr(4);
        expected = new GenomeListWithPseudocounts(
                ("TCTCGGGG\n" +
                "CCAAGGTG\n" +
                "TACAGGCG\n" +
                "TTCAGGTG"
                ).split("\n"));
        assertIterableEquals(expected,result);


    }
}