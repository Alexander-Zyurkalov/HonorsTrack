package org.stepic.bioinformatics;

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;


public class Probability {

    public static double logb( double a, double b ) {
        return Math.log(a) / Math.log(b);
    }
    public static double log2( double a ) {
        if (a == 0.0)  return 0;
        return logb(a,2);
    }

    public static double computeEntropy(Stream<Double> nums) {
        return - nums
            .map(p -> p * log2(p))
            .reduce(0.0, (s, n) -> s + n);
    }

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

    public static double probabilityByProfile(final String str,
                                              final HashMap<Character,List<Double>> profile){
        final List<Double> defaultArray =
                Stream.iterate(0, i -> i + 1).limit(str.length())
                .map(i -> 0.0)
                .collect(Collectors.toList());
        return Stream.iterate(0,i->i+1).limit(str.length())
                .map(i->profile.getOrDefault(str.charAt(i),defaultArray).get(i))
                .reduce(1.0,(m,n) -> m*n);
    }

    public static void main(String[] args) {

        var profile = new HashMap<Character, List<Double>>();
        profile.put('A', Arrays.asList(0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0));
        profile.put('C', Arrays.asList(0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6));
        profile.put('G', Arrays.asList(0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0));
        profile.put('T', Arrays.asList(0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4));
        System.out.println(
            probabilityByProfile("TCGTGGATTTCC",profile)
        );
    }
}
