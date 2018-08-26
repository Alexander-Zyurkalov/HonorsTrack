package org.stepic.bioinformatics;

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
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
    }
    public static BigDecimal  sumNonMutuallyExclusive(List<Double> probabilities) {
        return sumNonMutuallyExclusive(probabilities.stream().mapToDouble(n->n));
    }
    public static BigDecimal sumNonMutuallyExclusive(DoubleStream stream) {
        return stream.mapToObj(BigDecimal::valueOf)
                .reduce(Probability::sumNonMutuallyExclusive)
                .orElse(BigDecimal.ZERO);
    }


    public static ProfileList hashToList(final ProfileHash profile){
        final var profileList = new ProfileList();
        for (Character key :
                profile.keySet()) {
            var list = profile.get(key);
            for (int i = 0; i < list.size(); i++) {
                if (profileList.size() <= i) {
                    var hash = new HashMap<Character,Double>();
                    hash.put(key,list.get(i));
                    profileList.add(hash);
                }
                else {
                    profileList.get(i).put(key,list.get(i));
                }
            }
        }
        return profileList;
    }

    public static double probabilityByProfile(final String str,
                                              final ProfileHash profile) {
        return probabilityByProfile(str,hashToList(profile));
    }
    public static double probabilityByProfile(final String str,
                                              final ProfileList profile) {
        return Stream.iterate(0,i->i+1).limit(str.length())
                .map(i->profile.get(i).getOrDefault(str.charAt(i),0.0))
                .reduce(1.0,(m,n) -> m*n);
    }

    public static void main(String[] args) {

        var profile = new ProfileHash();
        profile.put('A', Arrays.asList(0.4, 0.3, 0.0, 0.1, 0.0, 0.9));
        profile.put('C', Arrays.asList(0.2, 0.3, 0.0, 0.4, 0.0, 0.1));
        profile.put('G', Arrays.asList(0.1, 0.3, 1.0, 0.1, 0.5, 0.0));
        profile.put('T', Arrays.asList(0.3, 0.1, 0.0, 0.4, 0.5, 0.0));
        System.out.println("AGGTGA = " +
            probabilityByProfile("AGGTGA",profile)
        );
        System.out.println("ACGTTA = " +
                probabilityByProfile("ACGTTA",profile)
        );
        System.out.println("ATGCTA = " +
                probabilityByProfile("ATGCTA",profile)
        );
        System.out.println("AAGTGA = " +
                probabilityByProfile("AAGTGA",profile)
        );
        System.out.println("ACGTTT = " +
                probabilityByProfile("ACGTTT",profile)
        );
        System.out.println("TCGCGA = " +
                probabilityByProfile("TCGCGA",profile)
        );
        System.out.println();
        System.out.println(probabilityByProfile("CAGTGA",profile));
    }
}
