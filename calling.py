import sys
import math
import pprint
import random
import scipy.misc as misc

data = [(1, 33),
        (1, 32),
        (1, 34),
        (1, 33),
        (1, 35),
        (1, 32),
        (1, 33),
        (0, 39),
        (0, 32),
        (0, 33),]

def generate_data(coverage, freq, mean_qual = 35, stdev = 2.0):
    data = []
    for i in range(coverage):
        if i < coverage*freq:
            base = 0
        else:
            base = 1

        qual = random.gauss(mean_qual, stdev)
        data.append((base, qual))

    return data

def qual_to_p(quality):
    return 10**(float(quality)/-10.0)

def p_to_qual(p_value):
    return math.log10(p_value) * -10.0

def binomial_distribution(total, success, p):
    return float(misc.comb(total, success, exact=True)) * p**float(success) * (1.0-p)**float(total-success)

def initialize_likelihood(data):
    genotype = 2
    likelihood = []
    for n in data:
        w = []
        for g in range(genotype):
            w.append(0.0)
        likelihood.append(w)

    return likelihood

def genotype_likelihood(data, likelihood):
    ''' Haploid-based genotype likelihood '''
    genotype = 2
    for n, d in enumerate(data):
        base = d[0]
        qual = d[1]
        for g in range(genotype):
            if base == 0:
                likelihood[n][g] = ((1.0-g) * (1.0-qual_to_p(qual))) + (g * qual_to_p(qual))
            elif base == 1:
                likelihood[n][g] = ((1.0-g) * qual_to_p(qual)) + (g * (1.0-qual_to_p(qual)))
    return likelihood

def sj_gt_config(n, data):
    counter = 0
    for d in data[0:n+1]:
        counter += d[0]
    return counter

def likelihood_of_allele_count(data, likelihood):
    ploidy = 1
    genotype = 2
    M = len(data) * ploidy
    lac = []

    for K in range(1, M+1):
        z = 0.00
        for m in range(M):
            total_ref_allele = sj_gt_config(m, data)
            for g in range(genotype):
                l = 0.00
                if K == total_ref_allele:
                    for m in range(M):
                        l += math.log10(misc.comb(ploidy, g, exact=True) * likelihood[m][g])

                z += l
        ll = 10**(z) / float(misc.comb(M, K))
        lac.append(ll)
        # print "[K=%s] Likelihood: %s" % (K, ll)

    return lac


def saf_em(data, likelihood):
    ploidy = 1
    genotype = 2
    saf = 0.3
    M = len(data) * ploidy
    counter = 0
    precision = 1.00

    max_iterations = 10
    threshold = 0.00001
    
    while counter < max_iterations and precision >= threshold:
        u = 0.00
        for n, d in enumerate(data):
            v = 0.00
            # calculate numerator
            for g in range(genotype):
                v += g * likelihood[n][g] * binomial_distribution(1, g, saf)

            w = 0.00
            # calculate denominator
            for g in range(genotype):
                w += likelihood[n][g] * binomial_distribution(1, g, saf)

            u += v/w
        new_saf = u/M
        precision = abs(new_saf - saf)
        saf = new_saf
        counter += 1
        
        # print "Iteration %s: site-allele-frequency = %s" % (counter, saf)
        
    return saf

def phi(n, k, saf):
    return binomial_distribution(n, k, saf)

def variant_quality(saf, lac, M):
    total = 0.00

    # Prior: Wright-Fisher Model
    for k in range(1, M+1):
        total+= phi(M, k, saf) * lac[k-1]

    return -10.0 * math.log10((phi(M, M, saf) * lac[M-1]) / total)
    
    
if __name__ == "__main__":
    data = generate_data(10, 0.3)
    ploidy = 1
    M = len(data)*ploidy
    
    likelihood = initialize_likelihood(data)
    # pprint.pprint(likelihood)
    
    likelihood = genotype_likelihood(data, likelihood)
    # pprint.pprint(likelihood)
    
    lac = likelihood_of_allele_count(data, likelihood)

    saf = saf_em(data, likelihood)

    print "Variant Quality: %s" % (variant_quality(saf, lac, M), )
