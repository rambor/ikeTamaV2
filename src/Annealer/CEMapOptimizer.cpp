//
// Copyright 2020 Robert P. Rambo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//
#include "../Anneal.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * Cross-entropy refinement of the lattice probability map (non-symmetric).
 *
 * Same algorithm as before (sample bead sets from the per-point probabilities, keep the topN by D_KL,
 * update probabilities with smoothing updateAlpha, same termination tests). Changes are for speed only:
 *
 *  1. Pair-bin table: the lattice is fixed for the whole run, so the P(r) bin of every pair of lattice points
 *     is computed once (M x M bytes). Scoring a trial is then N^2/2 table reads instead of N^2/2 square roots
 *     and bin conversions.
 *  2. Trials of a round are independent, so they run in parallel (OpenMP), each thread with its own random
 *     generator and bin counts. The topN of the round are chosen after all trials (std::partial_sort); this is
 *     the same elite set the serial insert-and-sort produced.
 *  3. Bookkeeping by lattice position (vectors) instead of std::map lookups by bead index.
 *  4. Sampling makes one shuffled pass per sweep over the points not yet chosen, as before, but without
 *     repeated std::iter_swap bookkeeping over the full lattice when N is reached early.
 */
void Anneal::ceMapOptimization(PointSetModel *pModel, PofRData *pData, unsigned int topN, const unsigned int minN, const unsigned int maxN, std::vector<ProbabilityBead> & lattice){

    const double invTotal = 1.0/(double)topN;
    const double updateAlpha = 0.67;

    maxbin = pModel->getMaxBin(pData);
    char score[9];
    (pData->getIsPr()) ? std::sprintf(score, "    D_KL") : std::sprintf(score, "CHI_FREE");
    maxbin += 1;
    totalBins = pData->getShannonBins();
    maxbin = (maxbin > totalBins) ? maxbin : totalBins;
    contactCutOff = interconnectivityCutOff;
    pData->setScoringFunction(maxbin);

    std::sort(lattice.begin(), lattice.end()); // sort lattice based on probabilities (as before)
    const auto M = (unsigned int)lattice.size();
    const unsigned int last = topN - 1;

    /*
     * pair-bin table over lattice positions (symmetric, full matrix for simple indexing)
     */
    if (maxbin > 255) throw std::runtime_error("ceMapOptimization: more than 255 P(r) bins, widen the pair-bin table type");
    std::vector<unsigned char> pairBin((size_t)M*M, 0);
    #pragma omp parallel for schedule(static)
    for (int a = 0; a < (int)M; a++){
        const vector3 & va = pModel->getBead(lattice[a].index)->getVec();
        for (unsigned int b = a + 1; b < M; b++){
            const unsigned int bin = pData->convertToBin((va - pModel->getBead(lattice[b].index)->getVec()).length());
            pairBin[(size_t)a*M + b] = pairBin[(size_t)b*M + a] = (unsigned char)bin;
        }
    }

    std::cout << "*******************             MAP REFINEMENT             *******************" << std::endl;
    logger("MAX BIN", std::to_string(maxbin));
    logger("LATTICE POINTS", std::to_string(M));
    logger("BOUNDS LOWER", std::to_string(minN));
    logger("BOUNDS UPPER", std::to_string(maxN));
    logger("TOTAL ROUNDS", std::to_string(highTempRounds));
    logger("TRIALS PER ROUND", std::to_string(ccmultiple));
    logger("TOPN", std::to_string(topN));
#ifdef _OPENMP
    logger("THREADS", std::to_string(omp_get_max_threads()));
#endif
    std::cout << "*******************       CROSS ENTROPY OPTIMIZATION       *******************" << std::endl;

    const unsigned int trials = ccmultiple;
    std::vector<double> trialEnergy(trials);
    std::vector<std::vector<unsigned int>> trialPositions(trials);   // lattice positions (not bead indices)
    std::vector<unsigned int> order(trials);
    std::vector<unsigned int> useCount(M);
    std::vector<double> prob(M);

    std::vector<double> running_average(3, 0.0);
    const double currentmaximum = 0.50;
    const unsigned int seed0 = seedForRun(0);

    for (unsigned int round = 0; round < 10*highTempRounds; round++){

        for (unsigned int i = 0; i < M; i++) prob[i] = lattice[i].prob;
        const double t0 =
#ifdef _OPENMP
            omp_get_wtime();
#else
            (double)std::clock()/CLOCKS_PER_SEC;
#endif

        #pragma omp parallel
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            std::mt19937 gen(seed0 + 7919u*round + 104729u*(unsigned int)tid);
            std::uniform_real_distribution<double> uniform(0.0, 1.0);
            std::uniform_int_distribution<unsigned int> randomN(minN, maxN - 1);
            std::vector<unsigned int> pool(M), chosen; chosen.reserve(maxN);
            std::vector<unsigned int> binCount(maxbin);

            #pragma omp for schedule(dynamic, 16)
            for (int trial = 0; trial < (int)trials; trial++){
                for (unsigned int i = 0; i < M; i++) pool[i] = i;
                std::shuffle(pool.begin(), pool.end(), gen);

                const unsigned int targetN = randomN(gen);
                unsigned int remaining = M;
                chosen.clear();
                // sweep the shuffled pool, accepting each point with its probability, until targetN are chosen
                while (chosen.size() < targetN && remaining > 0){
                    unsigned int before = (unsigned int)chosen.size();
                    for (unsigned int rn = 0; rn < remaining && chosen.size() < targetN; ){
                        const unsigned int pos = pool[rn];
                        if (uniform(gen) < prob[pos]){
                            chosen.push_back(pos);
                            pool[rn] = pool[--remaining];      // remove from the pool, test the swapped-in point next
                        } else {
                            rn++;
                        }
                    }
                    if (chosen.size() == before && remaining > 0){   // all probabilities ~0: avoid spinning
                        bool any = false;
                        for (unsigned int rn = 0; rn < remaining; rn++) if (prob[pool[rn]] > 1e-12) { any = true; break; }
                        if (!any) break;
                    }
                }

                // P(r) of the trial from the pair-bin table (sorted positions: rows read in increasing order, cache friendly)
                std::sort(chosen.begin(), chosen.end());
                std::fill(binCount.begin(), binCount.end(), 0);
                const unsigned int n = (unsigned int)chosen.size();
                for (unsigned int a = 0; a < n; a++){
                    const unsigned char * row = &pairBin[(size_t)chosen[a]*M];
                    for (unsigned int b = a + 1; b < n; b++) binCount[row[chosen[b]]]++;
                }
                trialEnergy[trial] = pData->getScore(binCount);     // strictly Kullback-Leibler divergence
                trialPositions[trial] = chosen;
            }
        }

        // elite set: topN lowest energies of this round
        for (unsigned int t = 0; t < trials; t++) order[t] = t;
        const unsigned int keep = std::min(topN, trials);
        std::partial_sort(order.begin(), order.begin() + keep, order.end(),
                          [&](unsigned int x, unsigned int y){ return trialEnergy[x] < trialEnergy[y]; });
        const double firstValue = trialEnergy[order[0]], lastValue = trialEnergy[order[keep - 1]];

        const double t1 =
#ifdef _OPENMP
            omp_get_wtime();
#else
            (double)std::clock()/CLOCKS_PER_SEC;
#endif
        logger("TRIAL TIME", formatNumber((float)(t1 - t0), 6));

        double diff = lastValue - firstValue;
        double average = (running_average[0] + running_average[1] + running_average[2])/3.0;
        double epsilon = (average > 0) ? (std::abs(running_average[0] - lastValue)/3)/average : 1.0;
        std::printf(" RND: %2i 1st: %.1E LAST: %.1E GAP: %.1E EPSLON: %.2f \n", round, firstValue, lastValue, diff, epsilon);
        running_average[0] = running_average[1];
        running_average[1] = running_average[2];
        running_average[2] = lastValue;

        // count how often each lattice point appears among the elite
        std::fill(useCount.begin(), useCount.end(), 0);
        for (unsigned int e = 0; e < keep; e++){
            for (unsigned int pos : trialPositions[order[e]]) useCount[pos]++;
        }

        // smoothed probability update
        double maxprob = 0, minprob = 1;
        unsigned int highcount = 0;
        for (unsigned int i = 0; i < M; i++){
            auto & point = lattice[i];
            point.prob = updateAlpha*useCount[i]*invTotal + (1.0 - updateAlpha)*point.prob;
            if (point.prob > maxprob) maxprob = point.prob;
            if (point.prob < minprob) minprob = point.prob;
            if (point.prob > currentmaximum) highcount++;
        }

        logger("MAX PROBABILITY", formatNumber(maxprob,3));
        logger("MIN PROBABILITY", formatNumber(minprob,3));
        logger("MIN COUNT", std::to_string(minN));
        logger(">0.60 COUNT", formatNumber(highcount));

        if (round > 5 && (highcount > minN || (maxprob > 0.6 && epsilon < 0.001 && lastValue < 0.000091))){
            std::string line = "EPSILON " + formatNumber((float)epsilon,3) + " < 0.01";
            logger("TERMINATING", line);
            break;
        }
    }
}


/*
 * Cross-entropy refinement of the lattice probability map with symmetry.
 *
 * Same algorithm as before (sample bead sets from the per-point probabilities in one shuffled pass, keep the topN by
 * D_KL + beta*violations, smoothed update with alpha 0.61, same termination test). Changes:
 *
 *  1. Clash count fixed. calculateModelPrDistributionSymCE stores the coordinates interleaved per bead
 *     (bead, its mates, next bead, ...) but counted clashes as if the first workingLimit entries were the subunit
 *     and the rest the mates, i.e. the first half of the beads (with their mates) against the second half. Clashes are
 *     now counted between every bead and every symmetry mate (ClashTracker, same definition as getViolations).
 *  2. Pair-bin table: the bins of all S x S mate pairs of every pair of lattice points are computed once
 *     (if they fit in memory, otherwise computed per trial from the cached mates). Scoring a trial is table reads.
 *  3. Trials of a round run in parallel (OpenMP), each thread with its own random generator; the elite is chosen
 *     after the round (std::partial_sort) instead of sorting topTrials after every insertion.
 *  4. topN larger than the number of trials no longer reads past the end of the elite list.
 *
 * pData->getScore must be thread-safe (it only reads the target distribution); build with -DCE_SERIAL otherwise.
 */
void Anneal::ceMapOptimizationSym(PointSetModel *pModel, PofRData *pData, unsigned int topN, const unsigned int minN, const unsigned int maxN, std::vector<ProbabilityBead> & lattice){

    contactCutOff = pModel->getNeighborCutOffLimit();
    logger("Contact cutoff", formatNumber(contactCutOff, 3));
    logger("Interconnectivity", formatNumber(interconnectivityCutOff, 3));
    violation_limit = (float)(sqrt(3)*pModel->getBeadRadius());

    ensureSymmetryCache(pModel);
    ensureClashTable(pModel);
    const unsigned int S = symCacheSubUnits, N = symCacheN;

    maxbin = pModel->getMaxBin(pData);
    maxbin += 1;  // maximum bin index, so vector size must be +1
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);
    if (maxbin > 255) throw std::runtime_error("ceMapOptimizationSym: more than 255 P(r) bins, widen the pair-bin table type");

    std::sort(lattice.begin(), lattice.end());
    const auto M = (unsigned int)lattice.size();
    std::vector<unsigned int> latticeIndex(M);
    for (unsigned int i=0; i<M; i++) latticeIndex[i] = lattice[i].index;

    // bins of the pairs between the mates of one bead (counted for every bead in a trial)
    std::vector<unsigned char> selfBins;   // (S*(S-1)/2) per lattice position
    selfBins.reserve((size_t)M*S*(S-1)/2);
    for (unsigned int a=0; a<M; a++)
        for (unsigned int s=0; s<S; s++)
            for (unsigned int t=s+1; t<S; t++)
                selfBins.push_back((unsigned char)pData->convertToBin((symCoords[(size_t)s*N + latticeIndex[a]] - symCoords[(size_t)t*N + latticeIndex[a]]).length()));

    // all S x S mate pairs of every pair of lattice positions, if they fit in memory
    const size_t SS = (size_t)S*S;
    const bool useTable = (double)M*M*SS <= 1.0e9;
    std::vector<unsigned char> pairBins;
    if (useTable){
        pairBins.assign((size_t)M*M*SS, 0);
        #pragma omp parallel for schedule(dynamic, 8)
        for (int a=0; a<(int)M; a++){
            for (unsigned int b=a+1; b<M; b++){
                unsigned char * ab = &pairBins[((size_t)a*M + b)*SS];
                unsigned char * ba = &pairBins[((size_t)b*M + a)*SS];
                for (unsigned int s=0; s<S; s++)
                    for (unsigned int t=0; t<S; t++){
                        const unsigned char bin = (unsigned char)pData->convertToBin((symCoords[(size_t)s*N + latticeIndex[a]] - symCoords[(size_t)t*N + latticeIndex[b]]).length());
                        ab[s*S + t] = bin; ba[t*S + s] = bin;
                    }
            }
        }
    }

    std::cout << "*******************             MAP REFINEMENT             *******************" << std::endl;
    logger("SYMMETRY", pModel->getSymmetry());
    logger("MAX BIN", std::to_string(maxbin));
    logger("LATTICE POINTS", std::to_string(M));
    logger("BOUNDS LOWER", std::to_string(minN));
    logger("BOUNDS UPPER", std::to_string(maxN));
    logger("TOTAL ROUNDS", std::to_string(highTempRounds));
    logger("TRIALS PER ROUND", std::to_string(ccmultiple));
    logger("TOPN", std::to_string(topN));
    logger("CONTACT CUTOFF (Angstroms)", formatNumber(contactCutOff, 2));
    logger("VIOLATION LIMIT", formatNumber(violation_limit, 2));
    logger("PAIR-BIN TABLE", useTable ? "YES" : "NO (memory)");
#ifdef _OPENMP
    logger("THREADS", std::to_string(omp_get_max_threads()));
#endif
    std::cout << "*******************       CROSS ENTROPY OPTIMIZATION       *******************" << std::endl;

    const unsigned int trials = ccmultiple;
    const unsigned int keep = std::min(topN, trials);
    std::vector<double> trialEnergy(trials);
    std::vector<std::vector<unsigned int>> trialPositions(trials);
    std::vector<unsigned int> order(trials), useCount(M);
    std::vector<double> prob(M);
    std::vector<double> running_average(3, 0.0);
    const unsigned int seed0 = seedForRun(0);
    const double updateAlpha = 0.61, invTotal = 1.0/(double)topN;

    for (unsigned int round=0; round < highTempRounds; round++){

        for (unsigned int i=0; i<M; i++) prob[i] = lattice[i].prob;
        const double t0 =
#ifdef _OPENMP
            omp_get_wtime();
#else
            (double)std::clock()/CLOCKS_PER_SEC;
#endif

#ifndef CE_SERIAL
        #pragma omp parallel
#endif
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            std::mt19937 gen(seed0 + 7919u*round + 104729u*(unsigned int)tid);
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            std::uniform_int_distribution<unsigned int> randomIndex(minN, maxN-1);
            std::vector<unsigned int> pool(M), chosen, beads; chosen.reserve(maxN); beads.reserve(maxN);
            std::vector<unsigned int> binCount(maxbin);
            std::vector<uint8_t> mark(N, 0);

#ifndef CE_SERIAL
            #pragma omp for schedule(dynamic, 16)
#endif
            for (int trial=0; trial<(int)trials; trial++){
                for (unsigned int i=0; i<M; i++) pool[i] = i;
                std::shuffle(pool.begin(), pool.end(), gen);
                const unsigned int currentN = randomIndex(gen);
                chosen.clear();
                for (unsigned int i=0; i<M; i++){             // one pass, accept with probability, stop at currentN
                    if (distribution(gen) < prob[pool[i]]){
                        chosen.push_back(pool[i]);
                        if (chosen.size() == currentN) break;
                    }
                }
                std::sort(chosen.begin(), chosen.end());
                const unsigned int n = (unsigned int)chosen.size();

                std::fill(binCount.begin(), binCount.end(), 0);
                for (unsigned int x=0; x<n; x++){
                    const unsigned int a = chosen[x];
                    for (size_t k=0; k<(size_t)S*(S-1)/2; k++) binCount[selfBins[(size_t)a*S*(S-1)/2 + k]]++;
                    if (useTable){
                        const unsigned char * row = &pairBins[(size_t)a*M*SS];
                        for (unsigned int y=x+1; y<n; y++){
                            const unsigned char * ab = row + (size_t)chosen[y]*SS;
                            for (size_t k=0; k<SS; k++) binCount[ab[k]]++;
                        }
                    } else {
                        for (unsigned int y=x+1; y<n; y++){
                            const unsigned int b = chosen[y];
                            for (unsigned int s=0; s<S; s++)
                                for (unsigned int t=0; t<S; t++)
                                    binCount[pData->convertToBin((symCoords[(size_t)s*N + latticeIndex[a]] - symCoords[(size_t)t*N + latticeIndex[b]]).length())]++;
                        }
                    }
                }

                beads.resize(n);
                for (unsigned int x=0; x<n; x++) beads[x] = latticeIndex[chosen[x]];
                const unsigned int violations = clashes.countWith(mark, beads.data(), n);

                trialEnergy[trial] = pData->getScore(binCount) + beta*violations;
                trialPositions[trial] = chosen;
            }
        }

        for (unsigned int t=0; t<trials; t++) order[t] = t;
        std::partial_sort(order.begin(), order.begin() + keep, order.end(),
                          [&](unsigned int x, unsigned int y){ return trialEnergy[x] < trialEnergy[y]; });
        const double firstValue = trialEnergy[order[0]], lastValue = trialEnergy[order[keep-1]];

        const double t1 =
#ifdef _OPENMP
            omp_get_wtime();
#else
            (double)std::clock()/CLOCKS_PER_SEC;
#endif
        const double diff = lastValue - firstValue;
        const double average = (running_average[0] + running_average[1] + running_average[2])/3.0;
        const double diffave = std::abs(average - lastValue)/average;

        printf(" RND: %2i 1st: %.1E LAST: %.1E GAP: %.1E EPSLON: %.2f TIME(SEC): %.1f \n", round, firstValue, lastValue, diff, diffave, t1 - t0);

        running_average[0] = running_average[1];
        running_average[1] = running_average[2];
        running_average[2] = lastValue;

        // count how often each lattice point appears among the elite
        std::fill(useCount.begin(), useCount.end(), 0);
        for (unsigned int e=0; e<keep; e++) for (unsigned int pos : trialPositions[order[e]]) useCount[pos]++;

        // smoothed probability update
        double maxprob = 0, minprob = 1;
        for (unsigned int i=0; i<M; i++){
            auto & point = lattice[i];
            point.prob = (updateAlpha*(double)useCount[i])*invTotal + (1.0 - updateAlpha)*point.prob;
            if (point.prob > maxprob) maxprob = point.prob;
            if (point.prob < minprob) minprob = point.prob;
        }

        if (diffave < 0.05 && round > 5){
            std::string line = "EPSILON " + formatNumber(diffave,3) + " < 0.05";
            logger("TERMINATING",line);
            break;
        }
        logger("MAX PROBABILITY", formatNumber(maxprob,3));
        logger("MIN PROBABILITY", formatNumber(minprob,3));
    }
}
