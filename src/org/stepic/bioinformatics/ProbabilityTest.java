package org.stepic.bioinformatics;

import org.junit.jupiter.api.Test;

import java.math.BigDecimal;

import static org.junit.jupiter.api.Assertions.*;

class ProbabilityTest {

    @Test
    void pow() {
        assertEquals(2*2,Probability.pow(2,2));
        assertEquals(2*2*2,Probability.pow(2,3));
        assertEquals(4*4*4*4,Probability.pow(4,4));
        assertEquals(5*5*5*5*5,Probability.pow(5,5));
        assertEquals(1,Probability.pow(1,5));
        assertEquals(0,Probability.pow(0,5));
        assertEquals(6*6*6*6*6*6,Probability.pow(6,6));
    }

    @Test
    void binomialCoefficient() {
        assertEquals(10,Probability.binomialCoefficient(5,2),0.0001);
        assertEquals(816,Probability.binomialCoefficient(18,3),0.0001);
    }

    @Test
    void factorial() {
        assertEquals(BigDecimal.valueOf(1),Probability.factorial(0));
        assertEquals(BigDecimal.valueOf(1),Probability.factorial(1));
        assertEquals(BigDecimal.valueOf(2),Probability.factorial(2));
        assertEquals(BigDecimal.valueOf(2*3),Probability.factorial(3));
        assertEquals(BigDecimal.valueOf(2*3*4),Probability.factorial(4));
        assertEquals(BigDecimal.valueOf(2*3*4*5),Probability.factorial(5));
        assertEquals(BigDecimal.valueOf(2*3*4*5*6),Probability.factorial(6));
        assertEquals(BigDecimal.valueOf(2*3*4*5*6*7),Probability.factorial(7));
        assertEquals(new BigDecimal("2432902008176640000"),Probability.factorial(20));

    }

    @Test
    void pr() {
        assertEquals(7.599592208862305E-7,Probability.pr(30,4,5,3),0.001);
        assertEquals(270.0/2187.0, Probability.pr(7,3,2,2),0.00001);
    }
}