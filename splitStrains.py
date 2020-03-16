import numpy as np
import matplotlib.pyplot as plt
import math
from sklearn import mixture
import pysam
import os
import argparse
import seaborn as sns
import sys
import logging
from mixem.distribution.distribution import Distribution
import mixem
from scipy.stats import binom
from scipy.stats import chi2
from scipy.stats import multinomial
from scipy.stats import norm
from scipy.stats import entropy
from scipy.special import factorial, comb, xlogy

from scipy.optimize import linprog
from scipy.optimize import minimize, Bounds


TITLE_FONT_SIZE = 8
TICK_FONT_SIZE = 8
AXES_FONT_SIZE = 8
LABEL_FONT_SIZE = 8
DPI = 300
PLOT_ENTROPY = False

sns.set(style="darkgrid")
sns.set_context("paper", rc={"lines.linewidth": 1, "patch.linewidth" : 0.5 })

np.set_printoptions(suppress=True, precision=2)

plt.rc('xtick', labelsize=TICK_FONT_SIZE)
plt.rc('ytick', labelsize=TICK_FONT_SIZE)
plt.rc('axes', labelsize=AXES_FONT_SIZE)

def roundUP(vec):
    """ This function rounds up means to integers so that the result sums to 100.
        If the means are not in proportion i.e. 70 and 10, then LP solver yields an error. """
    compons = len(vec)
    c = np.array([1 for i in range(0,compons)] + [-100])
    A = -c
    b = np.zeros((1,1))
    r = 1
    intVec = vec.astype(int)
    bound = [(val-r, val+r) for val in intVec] + [(1,1)]
    res = linprog(c, A_ub=A, b_ub=b, bounds=bound)
    return res['x'][:-1]


def plotConvergence(estimated_p, estimated_err):
    """ This function produces a plot which displays the estimated values of error and p for alternative hypothesis """

    plt.figure(figsize=(15,10))
    plt.subplot(211)
    plt.title('Converges of p')
    plt.plot(np.arange(len(estimated_p)),estimated_p)

    plt.subplot(212)
    plt.title('Converges of error')
    plt.plot(np.arange(len(estimated_err)), estimated_err)
    # plt.show()


# objective function for error
def obj_func_err(er, p, M, m, c_1, c_2):
    """ Objective function for optimization of error """
    res = -(M*np.log(p*(1 - 3*er) + (1-p)*er) + m*np.log((1-p)*(1-3*er)+p*er) + (c_1 + c_2 )*np.log(er))
    return np.sum(res)

# objective function for p
def obj_func_p(p, er, M, m, c_1, c_2):
    """ Objective function for optimization of proportion of a major """
    res = -(M*np.log(p*(1 - 3*er) + (1-p)*er) + m*np.log((1-p)*(1-3*er)+p*er) + (c_1 + c_2)*np.log(er))
    return np.sum(res)


# iterative optimization
def optimize(num_iter, init_p, init_err, M, m, c_1, c_2):
    """ This method is used to iteratively optimize obj_func_p and obj_func_err. The estimated parameters are used for alternative hypothesis test"""
    estimated_err = [init_err]
    estimated_p = [init_p]

    p_bound = Bounds(0.5, 0.9)
    e_bound = Bounds(1e-100, 0.3)
    opts = {'disp':False, 'gtol':1e-100}
    for i in range(num_iter):
        # estimate proportion of a major strain given estimated error
        res_p = minimize(obj_func_p, [init_p], args = (estimated_err[i], M, m, c_1, c_2), method='TNC',bounds=p_bound, options=opts)
        estimated_p.append(res_p.x[0])
        # estimate error given proportion
        res_err = minimize(obj_func_err, [init_err], args = (estimated_p[i], M, m, c_1, c_2), method='TNC', bounds=e_bound, options=opts)
        estimated_err.append(res_err.x[0])
    return [estimated_p, estimated_err]


def likelyhood_ratio_test(freqVec, upperLimit, num_iter=30, init_p=0.7, init_err=0.001):
    """ This function runs likelyhood ratio test: Under null hypothesis we assume a single strain, under alternative we assume 2 strains are present """
    # convert proportion to count
    depthVec = freqVec[:,-1]
    # make all proportions positive and sort them
    countVec = np.absolute(freqVec[:,0:-2])
    countVec = np.sort(countVec)

    # filter out noise
    depthVec = depthVec[countVec[:,-1] < upperLimit]
    countVec = countVec[countVec[:,-1] < upperLimit]

    countVec = countVec / 100 * depthVec[:, np.newaxis]
    major = countVec[:,3].astype(int)
    minor = countVec[:,2].astype(int)
    error1 = countVec[:,1].astype(int)
    error2 = countVec[:,0].astype(int)
    trials = major + minor + error1 + error2

    # estimate error for null hypo
    null_err = np.sum(trials - major) / (3*np.sum(trials))
    logging.debug(f'null hypo error estimate: {null_err}')

    # likelihood null hypo
    errors = error1 + error2 + minor
    null_hypo = np.sum(np.log(comb(trials, major)) + xlogy(errors,null_err) + xlogy(major,1 - 3*null_err))
    # estimate error and probability of major for alternative hypo
    estimated_p, estimated_err = optimize(num_iter, init_p, init_err, major, minor, error1, error2)
    # plotConvergence(estimated_p, estimated_err)
    alt_p = estimated_p[-1]
    alt_err = estimated_err[-1]
    logging.debug(f'alt hypo error estimate: {alt_err}')
    logging.debug(f'alt hypo p estimate: {alt_p}')
    # compute alternative hypothesis likelihood
    p_major = alt_p*(1-3*alt_err) + (1-alt_p)*alt_err
    p_minor = (1-alt_p)*(1-3*alt_err) + alt_p*alt_err
    errors = (error1 + error2)
    quantiles = np.array([major, minor, errors]).T
    category_probab_vec = [p_major, p_minor, alt_err]

    alt_hypo = np.sum(multinomial.logpmf(quantiles, n=trials, p=category_probab_vec))

    logging.debug(f'null hypothesis log likelihood: {null_hypo}')
    logging.debug(f'alter hypothesis log likelihood: {alt_hypo}')

    # Run Pearson's chi2 test on likelihood ratio statistic
    df = 2
    alpha = 0.95
    scale = 100
    critical = chi2.ppf(alpha, df, loc=0, scale=scale)
    log_ratio = -2*(null_hypo - alt_hypo)
    tresh = critical

    logging.info(f'Likelihood Ratio Statistic: -2*log(LR) = {log_ratio}, treshold: {int(tresh)}')
    # logging.info(f'Likelihood Ratio Statistic: -2*log(LR) = {int(log_ratio)}, treshold: {int(tresh)}')

    strainType = 'mixed'

    if  log_ratio < tresh:
        strainType = 'single'

    return strainType


def fitDataGMM(data, components=2):
    """ Fit provided data and return the Gaussian Mixture Model """
    data = data.reshape(-1,1)
    gmm = mixture.GaussianMixture(components, max_iter=150, n_init=20)
    gmm.fit(data)
    return gmm


def fitDataBMM(data, depth, lowerLimit, upperLimit, init_proportions, components=2):
    """ Fit data and return Binomial Mixture Model"""

    if components > 9:
        logging.error('Too many components specified. Max components 9.')
        exit()

    distros = []

    # Init distros
    for i in range(components):
        distros.append(BinomialDistribution(init_proportions[i], depth))

    # Format data as pairs of success and trials
    data_points = []

    for x,y in zip(data[:,:-2].flatten(), np.repeat(data[:,-1],4)):

        # do filtering on each proportion
        if x > lowerLimit and x < upperLimit:
            charCount = x*0.01*y    # convert proportion to count
            data_points.append([charCount,y])

        else:
            continue

    data = np.array(data_points)
    weights, distros, log_like = mixem.em(data, distros, initial_weights=None, progress_callback=None, max_iterations=500, tol_iters=200, tol=0.1)
    return BinomialMixture(weights, distros, log_like)


class BinomialDistribution(Distribution):
    """Binomial distribution with parameters p and n."""

    def __init__(self, p, n):
        self.p = p
        self.n = n

    def log_density(self, data):
        """ This function is called by EM algorithm. """
        x = data[:,0].astype(int)
        n = data[:,1].astype(int)
        # return = np.log(binom.pmf(x, n, self.p))
        # when computing probability we ignore the binomial coeff np.log(comb(n,x))
        return x*np.log(self.p) + (n - x)*np.log(1 - self.p)


    def estimate_parameters(self, data, weights):
        """ This method is used by EM algorithm to compute weights and estimate p for a mixture model."""
        self.p = np.sum(data[:,0]*weights)/(np.sum(weights)*self.n)

    def log_likelihood(self, x, n):
        """ Call this method for non EM use. """
        # return np.log(binom.pmf(x, n, self.p))
        # when computing probability we ignore the binomial coeff np.log(comb(n,x))
        return x*np.log(self.p) + (n - x)*np.log(1 - self.p)

    def __repr__(self):
        return "Binom[p={p:.4g}, trials={trials:.4g}]".format(p=self.p, trials=self.n)

    def get_proportion(self):
        return self.p


class BinomialMixture():
    """ Binomial Mixture Model class"""
    def __init__(self, weights, distros, log_like):
        self.weights = weights
        self.distros = distros
        self.log_like = log_like
        self.proportions = []

        for distro in self.distros:
            self.proportions.append(distro.get_proportion())

    def get_proportions(self):
        return np.array(self.proportions)

    def get_weights(self):
        return np.array(self.weights)

    def get_distros(self):
        return self.distros

    def get_num_components(self):
        return len(self.distros)

    def set_prob(self, list_p):
        for i in range(len(self.distros)):
            self.distros[i].p = list_p[i]


class Model():
    """ A general mixture model. The class provides a single interface for Binomial and Gaussian Mixtures."""
    def __init__(self, model):
        self.model = model
        self.means = None
        self.stdDivs = None
        self.distros = None
        self.modelName = ''

        if type(self.model) == BinomialMixture:
            self.modelName = 'BMM'
            self.distros = self.model.get_distros()
            self.means = self.model.get_proportions()*100

        elif type(self.model) == mixture.GaussianMixture:
            self.modelName = 'GMM'
            self.means = self.model.means_.flatten()
            self.means = self.means/np.sum(self.means)*100
            self.stdDivs = np.sqrt(gmm.covariances_.flatten())

    def get_weights(self):
        """ Get weights of each component."""

        if type(self.model) == BinomialMixture:
            return self.model.get_weights()

        elif type(self.model) == mixture.GaussianMixture:
            return [weight for weight in self.model.weights_]

    def component_log_density(self, x, component_index, depth=None):
        """ Compute log density of a component at an input data."""

        if type(self.model) == BinomialMixture:
            x = int(x*0.01*depth)
            return self.distros[component_index].log_likelihood(x, depth)

        elif type(self.model) == mixture.GaussianMixture:
            return norm.logpdf(x, self.means[component_index], self.stdDivs[component_index])

    def get_num_components(self):
        """ Get number of components."""
        if type(self.model) == BinomialMixture:
            return self.model.get_num_components()

        elif type(self.model) == mixture.GaussianMixture:
            return self.model.n_components

    def get_strain_proportions(self):
        """ Returns means for GaussianMixture or Binomial Mixture"""
        return np.array(self.means)

    def set_prob(self, p_list):
        if type(self.model) == BinomialMixture:
            model.set_prob(p_list)

    def __repr__(self):
        return self.modelName


class Interval():
    """ Simple class that contains genome interval start end info """
    def __init__(self, i_start, i_end):
        self.start = i_start
        self.end = i_end

    def isInside(self, pos):

        if pos > start and pos < end:
            return True

        return False


def getIntervals(gffFilePath, regionStart, regionEnd):
    """ Given a gff file this function produces a list of Interval objects """
    try:
        f = open(gffFilePath, 'r')
    except IOError as e:
        logging.error(e)
        return 1

    # the list of intervals
    intervals = []

    for line in f:
        lineParams = line.split()
        lineStart = lineParams[0]
        # skip headers
        if lineStart[0] == '#':
            continue

        # parse only gene regions
        feature = lineParams[2]
        if feature != 'gene':
            continue

        # try getting region
        try:
            start = int(lineParams[3]) - 1      # since sam file start indexing from 0 do the adjustment
            end = int(lineParams[4]) - 1
        except:
            logging.debug('failed to get a gff region. Skipping.')
            continue
        # stop iterating gff as we reached the end
        if start > regionEnd:
            break
        # skip if start is not in the range
        if start < regionStart:
            continue

        intervals.append(Interval(start, end))

    f.close()
    return intervals


def bayesClassifyReads(outputDir, freqVec, chromosome, samfile, refFile, model, components, lowerLimit, upperLimit,  step=150):
    """ Run the classification using Naive Bayes on a genome position of interest and write reads belonging to different strains in according files """
    outputSuffix = '_strain.reads'      # set output files suffix
    allele = {'a':0, 'c': 1, 't': 2, 'g':3}

    strain_proportions = model.get_strain_proportions()/np.sum(model.get_strain_proportions())*100

    for strain_proportion in strain_proportions:
        try:
            strainFile = open(f'{outputDir}/{str(int(strain_proportion))}{outputSuffix}', 'w')
            strainFile.close()
        except OSError as e:
            logging.error(e)
            return 1

    components = model.get_num_components()

    k = 0
    # start classification
    while k < len(freqVec):
        # locate SNP cluster
        snpList = [freqVec[k,4]]    # create a list which will hold a cluster of SNPs positions

        while k+1 < len(freqVec):

            if abs(freqVec[k,4]-freqVec[k+1,4]) < 1:
                snpList.append(freqVec[k+1,4])
                k += 1
            else:
                break

        readBuffer = dict()

        # iterate through each snp index in a cluster
        for snpPos in snpList:

            # get the base of the reference fasta at snpPos
            refBase = refFile.fetch(reference=chromosome, start=snpPos, end=snpPos + 1).lower()

            # get all reads at the snpPos positions belonging to a cluster
            for pileupcolumn in samfile.pileup(chromosome, snpPos, snpPos + 1, truncate=True, min_base_quality=baseQuality, min_mapping_quality=mapQuality, fastafile=refFile, stepper='samtools'):

                # for each read in the pileup column do the classification
                for read in pileupcolumn.pileups:

                    if type(read.query_position) is not int:
                        continue

                    # get the char at the query position
                    readBase = read.alignment.query_sequence[read.query_position].lower()

                    # ignore reads that match with reference file, this is valid only for more than 2 strains situation
                    if readBase == refBase and components > 2:
                        # print('skipping pos: ', snpPos, 'readBase: ', readBase)
                        continue

                    # get proportion of this char in freqVec
                    index = np.searchsorted(freqVec[:,4], snpPos)
                    proportion = freqVec[index, allele[readBase]]

                    # filter proportion
                    if proportion < 0:
                        # if proportion belongs to reference char then we need to flip the sign
                        proportion = -proportion

                    if proportion < 3 or proportion > 97:
                        continue

                    depth = freqVec[index, -1]

                    # if read has been seen before, find it in the buffer and update the probabilities
                    if read.alignment.query_name in readBuffer:
                        # compute likelyhoods for each component
                        for i in range(components):
                            readBuffer[read.alignment.query_name][0][i] += model.component_log_density(proportion, i, depth)

                        readBuffer[read.alignment.query_name][1] += readBase
                        readBuffer[read.alignment.query_name][2].append(snpPos)
                    # else, add the read to the buffer and initiate its probabilities at current snpPos
                    else:
                        readBasePos = [snpPos]

                        strainProbab = np.log(model.get_weights())

                        for i in range(components):
                            strainProbab[i] += model.component_log_density(proportion, i, depth)
                        readBuffer[read.alignment.query_name] = [strainProbab, readBase, readBasePos]


        # write each read pulled out during cluster processing to a file
        for read in readBuffer:

            strainName = strain_proportions[np.argmax(readBuffer[read][0])]

            strainFile = open(f'{outputDir}/{str(int(strainName))}{outputSuffix}', 'a')
            strainFile.write(
                read + ':   '
                + readBuffer[read][1] + ':    ' + str(int(proportion)) + ':    '
                + ''.join([' ' + str(int(x)) for x in readBuffer[read][2]]) + ':    '     # write each base posiion
                + ''.join([str(x) + ' ' for x in readBuffer[read][0]])     # write probabilities of each strain
                + '\n')
            strainFile.close()

        k += 1

    return 0


def computeDataFromSam(freqVec, samfile, refFile, baseQuality, mapQuality, regionStart=None, regionEnd=None):
    """ Use samtools mpileup engine to compute pileups and create an allele frequency vector for each position in the genome """
    columnCount = 1
    totalDepth = 0
    chromosome = samfile.references[0]

    # flag_filter=3 means we get reads that are mapped in a proper pair, flag_filter=8 ignore reads that have unmapped mate
    for pileupcolumn in samfile.pileup(chromosome, regionStart, regionEnd, truncate=True, min_base_quality=baseQuality, min_mapping_quality=mapQuality, fastafile=refFile, stepper='samtools', flag_require=3, flag_filter=3852):

        i = pileupcolumn.pos
        pile = ''.join(pileupcolumn.get_query_sequences(mark_ends=False)).lower()
        depth = len(pile)

        if depth == 0:
            continue

        totalDepth += depth
        columnCount += 1

        vec = []
        doNotAppend = 0
        refBase = refFile.fetch(chromosome, start=i, end=i + 1).lower()

        # Do filtering of each proportion vector derived from the column
        for char in ['a', 'c', 't', 'g']:

            charProp = 0
            # calculate char proportions at snp

            if char == refBase:
                # so if char equals reference char we set it to negative value so it is ignored in model fitting and historgram
                # however, we flip the sign to positive when doing classification in bayesClassifyReads()
                charProp = -pile.count(char) / depth * 100
            else:
                charProp = pile.count(char) / depth * 100

            # if no variation in the position do not consider it
            if abs(charProp) >= 98:
                doNotAppend = 1
                break

            vec.append(charProp)

        if doNotAppend == 1:
            continue

        vec.append(i)           # append position of the current column
        vec.append(depth)       # append depth of the current column
        freqVec.append(vec)     # append vector of proportions if it has variation

    # logging.info(f'avg depth {totalDepth/columnCount}')
    return freqVec


def plotScatter(outputDir, freqVec, originalFrecVec, figureFileName, entropyVec, regionStart, regionEnd, lowerLimit, upperLimit):
    """ Plot 2d scatter plot of SNP proportions on a genome"""
    # plt.figure(figsize = (7.5,4.25), dpi = DPI)
    fig, ax = plt.subplots()
    fig.set_size_inches(7.5, 6.2)

    ax.set_title('Base proportions at SNPs over genome', fontsize = TITLE_FONT_SIZE)
    ax.set_xlabel('genome positions')
    ax.set_ylabel('proportions')

    freqVec = np.absolute(freqVec[freqVec[:,4] > 0])
    originalFrecVec = originalFrecVec[originalFrecVec[:,4] > 0]
    fpoints = []
    points = []

    for i in range(0,4):
        fpoints.append(np.stack((freqVec[:,4], freqVec[:,i]), axis=1))
        points.append(np.stack((originalFrecVec[:,4], originalFrecVec[:,i]), axis=1))

    fpoints = np.concatenate(fpoints, axis=0)
    points = np.concatenate(points, axis=0)
    lowerLimit = 5
    upperLimit = 95
    fpoints = fpoints[fpoints[:,1] > lowerLimit]
    fpoints = fpoints[fpoints[:,1] < upperLimit]

    points = points[points[:,1] > lowerLimit]
    points = points[points[:,1] < upperLimit]

    if len(entropyVec) != 0:
        ax.step(entropyVec[:,0], entropyVec[:,1]*10, where='post', c='g')
        ax.scatter(points[:,0], points[:,1], c='r', marker='+')
        legend = ['10x entropy', 'original', 'survived']
    else:
        legend = ['proportions at SNPs']

    ax.scatter(fpoints[:,0], fpoints[:,1], c='b', s=9)

    ax.set_yticks(np.arange(0,101,10))
    ax.legend(legend, prop={'size': AXES_FONT_SIZE})

    try:
        plt.savefig(outputDir + '/' + figureFileName + '-scatter.png', dpi=DPI)
        # plt.show()
        plt.close()
        return 0
    except IOError as e:
        logging.error(e)
        return 1


def plotHist(outputDir, originalFreqVecFlat,  freqVecFlat, gmm, figureFileName):
    """ Plot the histogram and computed cdfs used for classification"""
    # Get estimated MLE and weights to build component pdfs
    means = gmm.means_.flatten()
    stdDiv = np.sqrt(gmm.covariances_.flatten())
    weights = gmm.weights_.flatten()

    title1 = f'Strain means {means.astype(int)} (log scale) '
    title2 = 'GMM pdf'

    fig, axs = plt.subplots(2,1)
    fig.set_size_inches(7.5, 6.25)
    fig.subplots_adjust(hspace=0.5)

    numBins = 100
    axs[0].hist(originalFreqVecFlat, bins=numBins, range=(0,numBins), alpha=0.5, facecolor='r')
    axs[0].hist(freqVecFlat, bins=numBins, range=(0,numBins), facecolor='b')
    axs[0].set_yscale('log', basey=2, nonposy='clip')
    axs[0].set_title(title1, fontsize=TITLE_FONT_SIZE)
    axs[0].set_ylabel("proportion frequency", fontsize=LABEL_FONT_SIZE)
    axs[0].set_xticks(range(0, numBins+1, 10))
    axs[0].legend(['original', 'processed'], prop={'size': AXES_FONT_SIZE})

    # Plot GMM pdf
    axs[1].set_title(title2, fontsize=TITLE_FONT_SIZE)
    axs[1].set_xticks(range(0, numBins+1, 10))
    x_axis = np.arange(0,numBins+1, 0.5).reshape(-1,1)

    # axs[1].plot(x_axis, np.exp(gmm.score_samples(x_axis))) # compute log likelyhood of each sample

    for params in zip(means,stdDiv):
        axs[1].plot(x_axis, norm.pdf(x_axis, params[0], params[1])) # compute log likelyhood of each sample

    axs[1].set_xlabel('proportions')

    try:
        plt.savefig(outputDir + '/' + figureFileName + '-hist.png', dpi=DPI)
        # plt.show()
        plt.close()
        return 0
    except Exception as e:
        logging.error(e)
        return 1


def computeData(filename):
    """
    Depricated method.
    Build matrix of proportions and a vector of features out of a pileup file. This method requires some preprocessing for pileup.
    """
    datafile = np.genfromtxt(filename, dtype='str')

    vecLen = np.vectorize(len)
    datafile[:,1] = np.char.lower(datafile[:,1])
    posIndexCol = datafile[:,0].reshape(len(datafile[:,0]),1)
    depthAtPos = vecLen(datafile[:,1])
    countMatrix = np.dstack((
        np.char.count(datafile[:,1],'a'),
        np.char.count(datafile[:,1],'c'),
        np.char.count(datafile[:,1],'t'),
        np.char.count(datafile[:,1],'g')))[0]
    data = np.divide(countMatrix[:], depthAtPos[:,None])*100     # percentages of all alternative alleles for each position in the pile
    freqVec = data.flatten()
    freqVec = freqVec[freqVec > TRESHOLD]
    data = np.append(posIndexCol, data, axis=1).astype(np.float_)
    return [freqVec, data]


def filterVec(freqVec, depthThreshold, ethreshold, entropy_step, lowerLimit, upperLimit):
    """ This function encapsulats 2 filtering steps: depth filtering and entropy filtering """
    # do depth filtering.
    freqVec =  freqVec[freqVec[:,-1] > depthThreshold]

    entropyVec = []

    if len(freqVec) < 2:
        logging.warning('No SNPs remained after depth filtering.')
        return (freqVec, entropyVec)

    # do entropy filtering
    if ethreshold != 0:
        sizeBeforeEntropy = len(freqVec)
        freqVec, entropyVec = entropyFilter(freqVec, ethreshold, lowerLimit, upperLimit, entropy_step)
        logging.info(f'max entropy: {max(entropyVec[:,1])}')
        logging.info(f'mean entropy: {entropyVec[:,1].mean()}')

    return (freqVec, entropyVec)


def entropyFilter(freqVec, threshold, lowerLimit, upperLimit, step=200):
    """ Given a freqVec, take intervals of length step an compute entropy for this interval. If entropy is higher than threshold remove the interval """
    bins = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 20, 30, 40, 50, 60, 70, 80, 84, 86, 88, 90, 92, 94,96, 98, 100])
    entropyVec = []
    sizeBeforeEntropy = len(freqVec)
    # make a copy of original freqVec in case filtering is too aggressive
    copyFreqVec = freqVec.copy()

    for i in range(0, len(copyFreqVec), step):
        freqVecFlat = copyFreqVec[i:i+step,:-2].flatten()      # take the interval of freqVec excluding position and depth and flatten
        freqVecFlat = freqVecFlat[freqVecFlat > lowerLimit]      # filter out values lower than lowerLimit
        freqVecFlat = freqVecFlat[freqVecFlat < upperLimit]      # filter out values higher than upperLimit

        pos = freqVec[i,4]

        if freqVecFlat.size == 0:
            entropyVec.append([pos, 0])
            continue

        distribution, bins = np.histogram(freqVecFlat, bins, density=True)      # compute distribution
        entropyVal = entropy(distribution)      # compute entropy
        entropyVec.append([pos, entropyVal])

        # if entropy high on the interval i to i + step subsitute all vectors with neg vectors
        if entropyVal > threshold:
            copyFreqVec[i:i+step,:-2] = np.array([-200,-200,-200,-200])

    copyFreqVec  = copyFreqVec[copyFreqVec[:,:-2].min(1) >= -100]  # filter out proportion entries which have [-200 -200 -200 -200]

    survivedProportion = int((len(copyFreqVec) / sizeBeforeEntropy)*100)
    entropyVec = np.array(entropyVec)

    if survivedProportion < 60:
        logging.warning(f'Entropy or Depth filtering is too aggressive, skipping. Consider lowering entropy step (-fes option) or depth threshold.')
        return (freqVec, entropyVec)

    logging.info(f'survived after entropy filtering: {survivedProportion}%')
    return (copyFreqVec, entropyVec)


def convolveVec(freqVecFlat, proprtionCountThresh=2, boxPoints=4):
    """ Convolve histogram to get rid of the noise """
    # setup a filter box
    box = np.ones(boxPoints)/boxPoints
    # bin the flattened vector
    hist, bins = np.histogram(freqVecFlat, np.arange(0,101, 1))
    # do convolution on the hist values
    histSmooth = np.convolve(hist, box, mode='same')

    # parse processed hist vector
    vec = []
    # for each index of hist
    for i in range(0, len(histSmooth)):
        count = int(histSmooth[i])
        # if index passes the count threshold then append this index to flat vector for gmm
        if count > proprtionCountThresh:
            for j in range(0, count):
                vec.append(i)

    return np.array(vec)


if __name__ == "__main__":

    description = """ SplitStrains detects minor/major strains and classify reads. """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-c','--classify', action='store_true', help='if this option is specified then the program will run reads classification, otherwise it will detect means and produce histogram png')
    parser.add_argument('-z','--reuse', action='store_true', help='if this flag is specified the program will reuse the csv file from the previous run')
    parser.add_argument('-mo', metavar='gmm/bmm', dest='model', type=str, help='specify clustering model: GMM or BMM. Default GMM.', default='gmm')
    parser.add_argument('-f', metavar="plotName", dest='plotName', default='plot', help='name for the histogram figure')
    parser.add_argument('-s', metavar='n', required=True, dest='regionStart', type=int, help='specify the start position on the genome')
    parser.add_argument('-e', metavar='n', required=True, dest='regionEnd', type=int, help='specify the end position on the genome')
    parser.add_argument('-r', metavar='ref', required=True, dest='ref', help='genome reference')
    parser.add_argument('-b', metavar='gff', dest='gff', help='use gff file to process only gff regions', default='')
    parser.add_argument('-o', metavar='dir', required=True, dest='outputDir', help='output directory')
    parser.add_argument('-i', metavar='n', default=150, dest='step', type=int, help='step for snp cluster detection. Default=150')
    parser.add_argument('-g', metavar='n', default=2, type=int, dest='components', help='gmm model components. Default=2')
    parser.add_argument('-ft', metavar='n', default=1, dest='proportion_count_threshold', help='Filter out proportions which have count less than n. Default=1')
    parser.add_argument('-fe', metavar='n', default=1, dest='entropy_thresh', help='Entropy filtering threshold. Set to 0 to turn off entropy filtering. Default=1.0')
    parser.add_argument('-fes', metavar='n', type=int, default=70, dest='entropy_step', help='Entropy filtering step. Defines the step length on freqVec.csv for entropy filtering computation. Default=200')
    parser.add_argument('-fd', metavar='n', required=True, default=100, dest='depthThreshold', type=int, help='Do not consider pileup columns with the depth less than n. Higher values help to reduce noise for gmm. Good values are avg depth of a bam file. Default=100')
    parser.add_argument('-u', metavar='n', type=int, default=90, dest='upperLimit', help='Do not consider proportion of bases beyond n value. Default=90')
    parser.add_argument('-l', metavar='n', type=int, default=10, dest='lowerLimit', help='Do not consider proportion of bases below n value. Default=10')
    parser.add_argument('-m', metavar='n', type=int, default=40, dest='mapQuality', help='Do not consider reads below n map quality. Default=40')
    parser.add_argument('-q', metavar='n', type=int, default=15, dest='baseQuality', help='Do not consider bases below n quality. Default=15')
    parser.add_argument(dest='bamFilePath', metavar='bamFilePath', help='input bam file')

    args = parser.parse_args()

    components = args.components    # gmm components. For 2 strains 2 components.
    proprtionCountThresh = args.proportion_count_threshold
    depthThreshold = args.depthThreshold      # pileup columns with depth less than filter value are skipped. Helps to reduce noise for gmm fitting
    lowerLimit = args.lowerLimit
    upperLimit = args.upperLimit
    regionStart = args.regionStart
    regionEnd = args.regionEnd
    step = args.step
    baseQuality = args.baseQuality    # samtools default mpileup quality filter is 13
    mapQuality = args.mapQuality
    outputDir = args.outputDir
    plotName = args.plotName
    refFastaPath = args.ref     # path to a ref fasta file
    bamFilePath = args.bamFilePath  # path to bam file
    gffFilePath = args.gff
    entropy_step = args.entropy_step
    ethreshold = float(args.entropy_thresh)
    useModel = args.model
    reuseFreqVec = args.reuse

    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO, stream=sys.stdout)

    try:
        samfile = pysam.AlignmentFile(bamFilePath, "rb" )     # read bam file
        refFile = pysam.FastaFile(refFastaPath)     # read reference fasta file
    except FileNotFoundError:
        logging.error(f'{bamFilePath} or {refFastaPath} is not found.')
        exit()

    refName = samfile.references[0]
    refLength = samfile.lengths[0]

    logging.info('splitStrain.py has started.')

    if (regionEnd > refLength):
        logging.warning('regionEnd > reference length.')

    interval = regionEnd - regionStart

    if interval < 1000000:
        logging.warning(f'the interval length {interval} is too small.')

    logging.info(f'sample name: {bamFilePath}')
    logging.info(f'reference name: {refName}, reference length: {refLength}')
    logging.info(f'regionStart: {regionStart}, regionEnd: {regionEnd}')
    logging.info(f'depth threshold: {depthThreshold}')
    logging.info(f'entropy threshold: {ethreshold}')

    intervals = []  # list of Interval objects. This will be populated if gff file is provided
    freqVec = []    # vector format [a prop, c prop, t prop, g prop, position, depth]
    freqVecCSV = 'freqVec.csv'

    # compute freqVec
    if reuseFreqVec == False:

        # If gff file is provided, compute on regions specified in a gff file
        if gffFilePath != '':
            logging.info(f'using gff: {gffFilePath}')
            intervals = getIntervals(gffFilePath, regionStart, regionEnd)
            for interval in intervals:
                freqVec = computeDataFromSam(freqVec, samfile, refFile, baseQuality, mapQuality, interval.start, interval.end)


        else:
            freqVec = computeDataFromSam(freqVec, samfile, refFile, baseQuality, mapQuality, regionStart, regionEnd)

        freqVec = np.array(freqVec)

        # terminate if freqVec has less than 2 entries
        if freqVec.size < 2:
            logging.warning('No SNPs found on the given interval.')
            exit()
        # write freqVec to a file
        try:
            np.savetxt(f'{outputDir}/{freqVecCSV}', freqVec, delimiter=',')
            # np.savetxt(f'{outputDir}/{freqVecCSV}', freqVec, delimiter=',', fmt='%i')
        except IOError:
            logging.error(f'failed to save the csv {outputDir}/{freqVecCSV}.')
            exit()
    # if reuse is set then load freqVec
    else:

        try:
            logging.info(f'loading csv {outputDir}/{freqVecCSV} from the previous run')
            freqVec = np.loadtxt(open(f'{outputDir}/{freqVecCSV}', 'rb'), delimiter=',', dtype=float)
            assert len(freqVec) != 0, f'{freqVecCSV} is empty.'

        except IOError:
            logging.error(f'failed to load the csv {outputDir}/{freqVecCSV}. Please check if the file exists.')
            exit()

        except AssertionError as error:
            logging.error(error)
            exit()


    logging.debug('Starting filterVec()')
    originalFreqVec = freqVec.copy()
    freqVec, entropyVec = filterVec(freqVec, depthThreshold, ethreshold, entropy_step, lowerLimit, upperLimit)
    plotScatter(outputDir, freqVec, originalFreqVec, plotName, entropyVec, regionStart, regionEnd, lowerLimit, upperLimit)


    num_iter = 20
    init_p = 0.7
    init_err = 0.001

    freqVec = freqVec[np.max(freqVec[:,:4], axis=1) < upperLimit]

    # call single strain if not enough variation is found
    if len(freqVec) < 5:
        logging.info(f'Not enough data points. result: {bamFilePath} Single strain.')
        exit()

    # test null and alt hypthesis
    strainType = likelyhood_ratio_test(freqVec, upperLimit, num_iter, init_p, init_err)

    # if test calls single strain exit
    if strainType == 'single':
        logging.info(f'LR test result: {bamFilePath} Single strain.')
        exit()


    if components == 2:
        # consider reference base frequencies in the histogram and fitting
        freqVecFlat = np.absolute(freqVec[:,:-2].flatten())
    else:
        # do not consider base frequencies. Ref bases frequencies will be filtered out since they are negative
        freqVecFlat = freqVec[:,:-2].flatten()

    freqVecFlat = freqVecFlat[freqVecFlat > lowerLimit]
    freqVecFlat = freqVecFlat[freqVecFlat < upperLimit]

    # TODO change box size to a parameter
    freqVecFlat = convolveVec(freqVecFlat, proprtionCountThresh, 3)

    if freqVecFlat.size < components:
        logging.info(f'Not enough SNP frequencies. result: {bamFilePath} Single strain.')
        exit()

    # Fit data with Gaussian Mixture
    gmm = fitDataGMM(freqVecFlat, components)
    init_proportions = gmm.means_.flatten()/100
    # Fit data with Binomial Mixture
    avgDepth = int(freqVec[:,-1].mean())
    bmm = fitDataBMM(freqVec, avgDepth, lowerLimit, upperLimit, init_proportions, components)
    bmm.set_prob(bmm.get_proportions()/np.sum(bmm.get_proportions()))
    # specify which model to use
    if useModel == 'bmm':
        model = Model(bmm)
    elif useModel == 'gmm':
        model = Model(gmm)
    else:
        logging.error('Wrong model name: Use either gmm or bmm.')
        exit()

    logging.info(f'using the model:{model}')


    means = model.get_strain_proportions()

    means = roundUP(means)

    if components == 2:
        if (means[0] > 50 and means[1] > 50) or (means[0] < 50 and means[1] < 50):
            logging.warning(f'result: Could not fit the data {bamFilePath}. Incorrect means:{means[0]}, {means[1]}. Possibly 50:50 split.')
            exit()

    logging.info(f'result: {bamFilePath} {means/np.sum(means)}')

    originalFrecVecFlat = originalFreqVec[:,:-2].flatten()
    originalFrecVecFlat = originalFrecVecFlat[originalFrecVecFlat > 2]
    originalFrecVecFlat = originalFrecVecFlat[originalFrecVecFlat < 98]
    plotHist(outputDir, originalFrecVecFlat, freqVecFlat, gmm, plotName)

    if args.classify == True:
        logging.info('starting strain separation')
        result = bayesClassifyReads(outputDir, originalFreqVec, refName, samfile, refFile, model, components, lowerLimit, upperLimit, step)

        if result == 0:
            logging.info('separation is complete.')
        else:
            logging.error('separation was not completed.')
