package org.stepic.bioinformatics;

import java.math.BigDecimal;
import java.util.List;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;


public class Probability {

    public static BigDecimal factorial(int n) {
        return Stream.iterate(1, i->i+1).limit(n)
                .map(BigDecimal::valueOf)
                .reduce(BigDecimal.valueOf(1), BigDecimal::multiply);
    }
    public static double binomialCoefficient (int n,int k){
        if (k>n) return 0;
        BigDecimal up = factorial(n);
        BigDecimal down = factorial(n-k).multiply(factorial(k));
        return up.divide(down).doubleValue();
    }

    public static long pow (int a, int b)
    {
        if ( b == 0)        return 1;
        if ( b == 1)        return a;
        if ( b % 2 == 0)    return     pow ( a * a, b/2); //even a=(a^2)^b/2
        else                return a * pow ( a * a, b/2); //odd  a=a*(a^2)^b/2

    }
    public static double pr(final int N, int A, final int k, final int t) {
        final double up = binomialCoefficient(N-t*(k-1),t);
        final long down = pow(A,t*k);
        return up/(double)down;
    }

    private static BigDecimal sumNonMutuallyExclusive(BigDecimal p1, BigDecimal p2){
        return p1.add(p2).add( p1.multiply(p2).negate());
//        return p1.multiply(p2).negate();
    }
    public static BigDecimal  sumNonMutuallyExclusive(List<Double> probabilities) {
        return sumNonMutuallyExclusive(probabilities.stream().mapToDouble(n->n));
    }
    public static BigDecimal sumNonMutuallyExclusive(DoubleStream stream) {
        return stream.mapToObj(BigDecimal::valueOf)
                .reduce(Probability::sumNonMutuallyExclusive)
                .orElse(BigDecimal.ZERO);
    }

    public static void main(String[] args) {
//        System.out.println();
        final var one_strand = pr(1000,4,9,1) * 1000;
        System.out.println(one_strand);
//        final var one_mer = BigDecimal.valueOf(0.25);
//        System.out.println(one_mer);
//        final var nine_mer = Stream
//                .iterate(one_mer,i->one_mer).limit(9)
//                .reduce(BigDecimal::multiply).orElse(one_mer);
//        System.out.println(nine_mer);

//        System.out.println(
//            sumNonMutuallyExclusive(
//                DoubleStream.iterate(one_strand,i->one_strand)
//                        .limit(500)
//            ).doubleValue() * 1000
//        );
    }
}
