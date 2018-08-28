package org.stepic.bioinformatics;

import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
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

    public static int getRandomIndexOld(final double ... probabilities){
        final int MAX_N = 100;
        final var rnd = new Random();
        final int chosen = rnd.nextInt(MAX_N)+1;
        final double sum = Arrays.stream(probabilities).sum();
        return Stream.iterate(0,i->i+1).limit(probabilities.length).sequential()
            .map( i -> new int[]{i,(int)Math.floor(probabilities[i]/sum * MAX_N)})
            .flatMap(n ->
                    Stream.iterate(0,i->i+1).limit(n[1])
                            .map(i->n[0]))
            .limit(chosen)
            .reduce((r,n) -> n).orElse(0);
    }
    public static int getRandomIndex(final double ... probabilities){
        final var rnd = new Random();
        final double sum = Arrays.stream(probabilities).sum();
        final double chosen = rnd.nextDouble()*sum;
        final double[] sums = new double[probabilities.length+1];
        sums[0] = 0.0;
        for (int i = 1; i < sums.length; i++) {
            sums[i] = (probabilities[i-1]/sum) + sums[i-1];
            if (chosen < sums[i])
                return i-1;
        }
        return 0;
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
        System.out.println("result = " + getRandomIndex(1.0/6,2.0/6,3.0/6));
        double p1 = (600.0 - 15) / (600-15+1);
        double p2 = 1 - p1;
        double counter = binomialCoefficient(10,2);
        var formatter = new DecimalFormat(".#######");
        System.out.println(formatter.format(Math.pow(p2,2)*Math.pow(p1,8) * counter));

    }
}
