import math
from collections import Counter

# ==================================================
# Statistical tools
# ==================================================

# efficient combinatorial function to handle extremely large numbers
def log_choose(n, k):
    r = 0.0
    # swap for efficiency if k is more than half of n
    if k * 2 > n:
        k = n - k

    for  d in xrange(1,k+1):
        r += math.log(n, 10)
        r -= math.log(d, 10)
        n -= 1

    return r

# return the genotype and log10 p-value
def bayes_gt(ref, alt, is_dup):
    # probability of seeing an alt read with true genotype of of hom_ref, het, hom_alt respectively
    if is_dup: # specialized logic to handle non-destructive events such as duplications
        p_alt = [1e-2, 0.2, 1/3.0]
    else:
        p_alt = [1e-3, 0.5, 0.9]

    total = ref + alt
    log_combo = log_choose(total, alt)

    lp_homref = log_combo + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
    lp_het = log_combo + alt * math.log(p_alt[1], 10) + ref * math.log(1 - p_alt[1], 10)
    lp_homalt = log_combo + alt * math.log(p_alt[2], 10) + ref * math.log(1 - p_alt[2], 10)

    return (lp_homref, lp_het, lp_homalt)

# get the number of entries in the set
def countRecords(myCounter):
    numRecords = sum(myCounter.values())
    return numRecords

# median is approx 50th percentile, except when it is between
# two values in which case it's the mean of them.
def median(myCounter):
    #length is the number of bases we're looking at
    numEntries = countRecords(myCounter)

    # the ordinal value of the middle element
    # if 2 middle elements, then non-integer
    limit = 0.5 * numEntries

    # a list of the values, sorted smallest to largest
    # note that this list contains unique elements only
    valueList = list(myCounter)
    valueList.sort()

    # number of entries we've gone through
    runEntries = 0
    # index of the current value in valueList
    i = 0
    # initiate v, in case list only has one element
    v = valueList[i]

    # move through the value list, iterating by number of
    # entries for each value
    while runEntries < limit:
        v = valueList[i]
        runEntries += myCounter[v]
        i += 1
    if runEntries == limit:
        return (v + valueList[i]) / 2.0
    else:
        return v

# calculate upper median absolute deviation
def upper_mad(myCounter, myMedian):
    residCounter = Counter()
    for x in myCounter:
        if x > myMedian:
            residCounter[abs(x - myMedian)] += myCounter[x]
    return median(residCounter)

# sum of the entries
def sumRecords(myCounter):
    mySum = 0.0
    for c in myCounter:
        mySum += c * float(myCounter[c])
    return mySum

# calculate the arithmetic mean, given a counter and the
# length of the feature (chromosome or genome)
# for x percentile, x% of the elements in the set are
# <= the output value
def mean(myCounter):
    # the number of total entries in the set is the
    # sum of the occurrences for each value
    numRecords = countRecords(myCounter)

    # u holds the mean
    u = float()

    u = sumRecords(myCounter) / numRecords
    return u

def stdev(myCounter):
    # the number of total entries in the set is the
    # sum of the occurrences for each value
    numRecords = countRecords(myCounter)

    # u holds the mean
    u = mean(myCounter)
    sumVar = 0.0

    # stdev is sqrt(sum((x-u)^2)/#elements)
    for c in myCounter:
        sumVar += myCounter[c] * (c - u)**2
    myVariance = float(sumVar) / numRecords
    stdev = myVariance**(0.5)
    return stdev

