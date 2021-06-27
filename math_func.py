#!/usr/bin/env python
# coding=utf-8

import math
from math import log10
import numpy as np

#sys.float_info.dig to show 
_MAX_CONFIDENCE = 1.0 - 1e-15

# rounded to, for numerical stability.
_GL_PRECISION = 50 
_QUAL_PRECISION = 1
LOG_10 = math.log(10.0)

prior = float(1/3)
Genotype = ["0/0", "0/1", "1/1", '0/2', '1/2', '2/2']

def perror_to_phred(perror):
    return -10 * math.log10(perror)

def ptrue_to_bounded_phred_dv(ptrue, max_prob = _MAX_CONFIDENCE):
    if not 0 <= ptrue <= 1:
        raise ValueError("ptrue must be between zero and one: {}".format(ptrue))

    return perror_to_phred(1.0 - min(ptrue, max_prob))


def computer_quals(GL_P, n_alleles):
    """Computes GQ and QUAL values from a set of prediction probabilities."""
    ## GQ is prob(GQ) / prob(all genotypes)
    ## GQ is rounded to the nearest integer to comply with the VCF spec
    #index, gt = get_max_index(predictions, n_alleles)
    #gq = int(np.around(ptrue_to_bounded_phred_dv(predictions[index])))
    ##print(min(sum(predictions[1:]), 1.0))
    #qual = ptrue_to_bounded_phred_dv(min(sum(predictions[1:]), 1.0))
    #rounded_qual = round(qual, _QUAL_PRECISION)

    #return gq, rounded_qual, gt
    index_of_max = np.argmax(GL_P)
    PL = [int(np.around(-10*log10(i))) for i in GL_P]
    #GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
    GQ = int(-10*log10(sum(GL_P[:index_of_max]+GL_P[index_of_max+1:])))
    QUAL = abs(np.around(-10*log10(GL_P[0]), 1))

    return Genotype[GL_P.index(max(GL_P))], GQ, QUAL


def get_max_index(predictions, n_alleles = 2):
    index_of_max = np.argmax(predictions)

    index = 0
    for h1 in range(0, n_alleles):
        for h2 in range(0, h1 + 1):
            if index == index_of_max:
                return index, [h2, h1]
            index += 1
    raise ValueError("No corresponding genotype for predictions", predictions)


def round_gls(gls, precision = _GL_PRECISION):
    """Returns genotype likelihoods rounded to the desired precision level.
    
    Returns genotype likelihoods rounded to the desired precision level.
    
    Args:
        gls: A list of floats. The input genotype likelihoods at any precision.
        precision: Positive int. The number of places past the decimal point to
                round to. If None, no rounding is performed.
    """
    if abs(sum(gls) - 1) > 1e-6:
        raise ValueError('Invalid genotype likelihoods do not sum to one: sum({}) = {}'.format(gls, sum(gls)))

    min_ix = 0
    min_gl = gls[0]
    for ix, gl in enumerate(gls):
        if gl < min_gl:
            min_gl = gl
            min_ix = ix

    rounded_gls = [round(gl, precision) for gl in gls]

    rounded_gls[min_ix] = max(
            0.0,
            round(1 - sum(rounded_gls[:min_ix] + rounded_gls[min_ix + 1:]), 
            precision))

    return rounded_gls

def toRealSpace(log10_probs):
    """
    sum(10^log10_probs) close to one
    """
    return [math.pow(10.0, x) for x in log10_probs]

def log10sumexp(log10_probs):
    """Returns log10(sum(10^log10_probs)) computed in a numerically-stable way"""

    m = max(log10_probs)

    return m + math.log10(sum(pow(10.0, x-m) for x in log10_probs))

def normalize_log10_probs(log10_probs):
    """
    Approximately normalizes log10 probabilities
    """
    #log10_probs = [math.log10(p) for p in probs]
    log10_probs = np.array(log10_probs)
    if np.max(log10_probs) > 0.0:
        raise ValueError("log10_probs all must be <=0", log10_probs)

    lse = log10sumexp(log10_probs)

    return np.minimum(log10_probs - lse, 0.0)

def calc_reference_confidence(n_ref, n_total, p_error):
    if n_ref < 0:
        raise ValueError("n_ref={} must be >=0".format(n_ref))
    if n_total < n_ref:
        raise ValueError("n_total={} must be >= n_ref={}".format(n_total, n_ref))
    n_alts = n_total - n_ref
    logp = math.log(p_error) / LOG_10
    log1p = math.log1p(-p_error) / LOG_10

    p_ref = n_ref * log1p + n_alts * logp
    p_het = -n_total * math.log(2) / LOG_10
    p_hom_alt= n_ref *logp + n_alts * log1p

    return normalize_log10_probs([p_ref, p_het, p_hom_alt])

def cal_GL(n_ref, n_total, err):
    # Approximate adjustment of events with larger read depth
    # original genotype likelihood
    n_alts = n_total - n_ref
    c0, c1 = n_ref, n_alts
    #ori_GL00 = np.float64(pow((1-err), c0)*pow(err, c1)*(1-prior)/2)
    #ori_GL11 = np.float64(pow(err, c0)*pow((1-err), c1)*(1-prior)/2)
    #ori_GL01 = np.float64(pow(0.5, c0+c1)*prior)

    #ori_GL00 = np.float64(pow((1-err), c0)*pow(err, c1))
    #ori_GL11 = np.float64(pow(err, c0)*pow((1-err), c1))
    #ori_GL01 = np.float64(pow(0.5, c0+c1))
    #log10_probs = list(normalize_log10_probs([log10(ori_GL00), log10(ori_GL01), log10(ori_GL11)]))

    logp = math.log(err) / LOG_10
    log1p = math.log1p(-err) / LOG_10
    ori_GL00 = n_ref * log1p + n_alts * logp
    ori_GL01 = -n_total * math.log(2) / LOG_10
    ori_GL11 = n_ref *logp + n_alts * log1p
    log10_probs = list(normalize_log10_probs([ori_GL00, ori_GL01, ori_GL11]))

    return toRealSpace(log10_probs)

def rescale_read_counts(n_ref, n_total, max_allowed_reads=100):
    """Ensures that n_total <= max_allowed_reads, rescaling if necessary."""
    if n_total > max_allowed_reads:
        ratio = n_ref / (1.0 * n_total)
        n_ref = int(math.ceil(ratio * max_allowed_reads))
        n_total = max_allowed_reads

    return n_ref, n_total
