package org.stepic.bioinformatics;

import org.junit.jupiter.api.Test;

import java.util.function.Function;
import java.util.stream.IntStream;

import static org.junit.jupiter.api.Assertions.*;

class SequenceTest {

    @Test
    void equalsAt() {
        assertTrue(new Sequence("GCGCG").equalsAt("GCG",0));
        assertFalse(new Sequence("GCGCG").equalsAt("GCG",1));
        assertTrue(new Sequence("GCGCG").equalsAt("GCG",2));
        assertFalse(new Sequence("GCGCG").equalsAt("GCG",3));

    }


    @Test
    void reverseComplement() {
        var seq = new Sequence("AAAACCCGGT");
        var expected = new Sequence("ACCGGGTTTT");
        assertEquals(expected,seq.reverseComplement());
    }

    @Test
    void patternToNumber() {
        int p = 4;
        var pow = new Sequence("ATGCAA").pow;
        assertEquals(1,pow.apply(0).intValue());
        assertEquals(p,pow.apply(1).intValue());
        assertEquals(p*p,pow.apply(2).intValue());
        assertEquals(p*p*p,pow.apply(3).intValue());
        assertEquals(p*p*p*p,pow.apply(4).intValue());
        assertEquals(p*p*p*p*p,pow.apply(5).intValue());

        assertEquals(0,new Sequence("AA").patternToNumber());
        assertEquals(1,new Sequence("AC").patternToNumber());
        assertEquals(2,new Sequence("AG").patternToNumber());
        assertEquals(3,new Sequence("AT").patternToNumber());
        assertEquals(4,new Sequence("CA").patternToNumber());
        assertEquals(5,new Sequence("CC").patternToNumber());
        assertEquals(6,new Sequence("CG").patternToNumber());
        assertEquals(7,new Sequence("CT").patternToNumber());
        assertEquals(8,new Sequence("GA").patternToNumber());
        assertEquals(9,new Sequence("GC").patternToNumber());
        assertEquals(10,new Sequence("GG").patternToNumber());
        assertEquals(11,new Sequence("GT").patternToNumber());
        assertEquals(12,new Sequence("TA").patternToNumber());
        assertEquals(13,new Sequence("TC").patternToNumber());
        assertEquals(14,new Sequence("TG").patternToNumber());
        assertEquals(15,new Sequence("TT").patternToNumber());
    }
    @Test
    void patternToNumberRecursively() {
        assertEquals(0,new Sequence("AA").patternToNumberRecursively());
        assertEquals(1,new Sequence("AC").patternToNumberRecursively());
        assertEquals(2,new Sequence("AG").patternToNumberRecursively());
        assertEquals(3,new Sequence("AT").patternToNumberRecursively());
        assertEquals(4,new Sequence("CA").patternToNumberRecursively());
        assertEquals(5,new Sequence("CC").patternToNumberRecursively());
        assertEquals(6,new Sequence("CG").patternToNumberRecursively());
        assertEquals(7,new Sequence("CT").patternToNumberRecursively());
        assertEquals(8,new Sequence("GA").patternToNumberRecursively());
        assertEquals(9,new Sequence("GC").patternToNumberRecursively());
        assertEquals(10,new Sequence("GG").patternToNumberRecursively());
        assertEquals(11,new Sequence("GT").patternToNumberRecursively());
        assertEquals(12,new Sequence("TA").patternToNumberRecursively());
        assertEquals(13,new Sequence("TC").patternToNumberRecursively());
        assertEquals(14,new Sequence("TG").patternToNumberRecursively());
        assertEquals(15,new Sequence("TT").patternToNumberRecursively());
    }
}