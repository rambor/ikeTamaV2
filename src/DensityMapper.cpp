//
// Created by xos81802 on 10/06/2021.
//
// Copyright 2021 Robert P. Rambo
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



#include "DensityMapper.h"


DensityMapper::DensityMapper(std::string filename, float qmax, float sampling_frequency) :
        sampling_frequency(sampling_frequency),
        delta_r( M_PI/(sampling_frequency*qmax)),
        cutoff((float)(M_PI/qmax)){

    PDBModel tempmodel (filename, true, false);

    kmax = std::ceil(tempmodel.getSMax()*1.2f);

    unsigned int totalCoords = tempmodel.getTotalCoordinates();

    // create lattice point for each point in supporting model (centered)coordinates)
    for(unsigned int i = 0; i<totalCoords; i++) {
        float pX = tempmodel.getCenteredXVec()[i];
        float pY = tempmodel.getCenteredYVec()[i];
        float pZ = tempmodel.getCenteredZVec()[i];
        centered_coordinates.emplace_back(vector3(pX, pY, pZ));
        lattice_points.emplace_back(LatticePoint(i, 3));
    }

    total_centered_coordinates = centered_coordinates.size();
    total_shells = (int)std::ceil(kmax/delta_r) + 1;

    for(int i=0; i<total_shells; i++){
        rvalues.push_back((delta_r*i));
    }

    lmax = (int)std::ceil(rvalues[total_shells-1]*qmax/M_PI) + 4;
    std::cout << "\n total shells " << total_shells << "  " << rvalues[total_shells-1] << std::endl;

    /*
     * calculate smallest distances in set
     * for each lattice point, determine closest point
     * report min and max and average
     */

    DminType dmin_values = SASTOOLS_UTILS_H::getDminValues(centered_coordinates);

    dmin_supremum = dmin_values.dmin_supremum;
    dmin_infimum = dmin_values.dmin_infimum;
    average_dmin = dmin_values.average_dmin;
    stdev_dmin = dmin_values.stdev_dmin;

    fwhm_sigma = 1.5f*cutoff/2.355f; //1.2f*cutoff/2.355f; // using the formula fwhm = 2.355*sigma which should account for 76% area in distribution

    logger("dmin SUPREMUM", formatNumber(dmin_supremum,2));
    logger("dmin INFINUM", formatNumber(dmin_infimum, 2));
    logger("dmin AVERAGE", formatNumber(average_dmin, 2));
    logger("dmin STDEV", formatNumber(stdev_dmin, 3));
    logger("sampling freq", formatNumber(sampling_frequency, 2));
    logger("delta r", formatNumber(delta_r, 3));
    logger("cutoff", formatNumber(cutoff, 2));
    logger("FWHM", formatNumber(fwhm_sigma, 2));

    amplitudes.resize(total_centered_coordinates);
    for(auto & amp : amplitudes){
        amp = 1.0;
    }
    // for each r_value, I can calculate a sphere of specific area :: 4*PI*r^2


//    std::vector<cl::Platform> all_platforms;
//    cl::Platform::get(&all_platforms);
//    if(all_platforms.size()==0){
//        std::cout<<" No platforms found. Check OpenCL installation!\n";
//        exit(1);
//    }
//    cl::Platform default_platform=all_platforms[0];
//    std::cout << "Using platform: "<<default_platform.getInfo<CL_PLATFORM_NAME>()<<"\n";


    /*
     * dispatch library manages the multiple runtimes that might be available on the computer (Intel CPU with NVIDIA GPU)
     *
     */
//    cl_uint platformIdCount = 0;
//    clGetPlatformIDs (0, nullptr, &platformIdCount);
//
//    std::vector<cl_platform_id> platformIds (platformIdCount);
//    clGetPlatformIDs (platformIdCount, platformIds.data (), nullptr);
//
//    cl_uint deviceIdCount = 0;
//    clGetDeviceIDs (platformIds [0], CL_DEVICE_TYPE_ALL, 0, nullptr,
//                    &deviceIdCount);
//
//    std::vector<cl_device_id> deviceIds (deviceIdCount);
//    clGetDeviceIDs (platformIds [0], CL_DEVICE_TYPE_ALL, deviceIdCount,
//                    deviceIds.data (), nullptr);
//
//
//    std::cout <<" " << std::endl;
//    logger("Total OpenCL Platforms", formatNumber(platformIdCount));
    /*
     * device[0] is CPU
     * device[1] is likely GPU
     */


//    clGetDeviceIDs (platformIds[0], CL_DEVICE_TYPE_ALL, 0, nullptr,
//                    &deviceIdCount);
//
//    if (deviceIdCount <2){
//        exit(1);
//    } else {
//
//    }


//    const cl_context_properties contextProperties [] = {
//                    CL_CONTEXT_PLATFORM,
//                    reinterpret_cast<cl_context_properties> (platformIds [0]),
//                    0, 0
//            };
//
//    cl_context context = clCreateContext (
//            contextProperties, deviceIdCount,
//            deviceIds.data (), nullptr,
//            nullptr, &error);
//

//    char *profile = NULL;
//    size_t size;

//    for(unsigned int i=0; i<platformIdCount; i++){
//        clGetDeviceIDs (platformIds[i], CL_DEVICE_TYPE_ALL, 0, nullptr,
//                        &deviceIdCount);
//        logger("Plaform ID", formatNumber(i));
//        clGetPlatformInfo(platformIds[i], CL_PLATFORM_PROFILE, NULL, profile, &size); // get size of profile char array
//        profile = (char*)malloc(size);
//        clGetPlatformInfo(platformIds[i], CL_PLATFORM_PROFILE,size, profile, NULL); // get profile char array
//        std::cout << profile << std::endl;
//        std::string platformVendor;
//        clGetPlatformInfo(platformIds[i], CL_PLATFORM_NAME,NULL, profile, &size); // get profile char array
//        profile = (char*)malloc(size);
//        clGetPlatformInfo(platformIds[i], CL_PLATFORM_NAME,size, profile, NULL); // get profile char array
//
//        std::cout << profile << std::endl;
//        logger("Total Devices", formatNumber(deviceIdCount));
//    }

    //set GPU

    createSphericalFibonacciLattice();
    setModel();
}


/*
 * want area per point to be nearly the same across all shells
 *
 */
void DensityMapper::createSphericalFibonacciLattice() {

    auto initial_area = (float)(delta_r*delta_r*4.0*M_PI);
    auto base_points = (float)std::ceil(delta_r*delta_r*4*M_PI/(sqrtf(3.0f)/4.0f*delta_r*delta_r)); // area of sphere divided by area equilateral triangle with sides delta_r

    float area_per_point = initial_area/base_points;

    for(int i=0; i<total_shells; i++){
        shells.emplace_back( Shell( rvalues[i], area_per_point) );
    }
}

/*
 * Need to set which shell lattice points overlap with model
 * for all others, set to false, we don't need to calculate as P(r at theta and phi) is zero
 */
void DensityMapper::setModel(){
    // which shell coordinates are not within cut-off centered_coordinates
    for(int i=0; i<total_shells; i++){
        /*
         * distance between consecutive shells is delta_r
         *
         */
        shells[i].setModel(centered_coordinates, 2*fwhm_sigma, fwhm_sigma);
        shells[i].populateSphericalHarmonics(lmax);
    }
}


void DensityMapper::calculateDensityCoefficientAtLMR(){
    // for each point in shell, calculate integral of
//    for(auto & shell : shells){
//        shell.updateDensityCoefficients(centered_coordinates, amplitudes, fwhm_sigma);
//    }
}


float DensityMapper::calculateIntensityAtQ(int q_index){

    float intensity = 0.0;
    float value;

    for(int l=0; l<=lmax; l++){
        for(int m=-l; m<=l; m++){
            value =calculateAmplitudeAtQ(q_index, l, m);
            //std::cout << " " << l << " " << m << " :: " << value << std::endl;
            intensity += value;
        }
    }

    return intensity;
}



/*
 * calculate electron density value at each point in the keptHCPlattice points (sampling lattice)
 *
 */
float DensityMapper::populateDensities(std::vector<float> & amplitudesT, std::vector<float> & densities_at_HCP, float * squared_amplitudes){
    unsigned int totalHCPInUse = keptHCPLatticePoints.size();

    const float * const pAmp = amplitudesT.data(); // base amplitude for input MODEL
    float * const pDen = densities_at_HCP.data();

    Neighbors * const  pNeighborhood = neighborhoods.data();

    unsigned int totalNeighbors;
    float amplitude_r_i;
    float sum = 0.0;

    for(unsigned int i=0; i<totalHCPInUse; i++){
        // calculate electron density amplitude based on neighborhood
        auto & neighborhood_i = pNeighborhood[i];
        const float * model_kernel_distances = neighborhood_i.getModelDistances();
        const unsigned int * model_neighbors = neighborhood_i.getModelIndices();

        totalNeighbors = neighborhood_i.getTotalNeighbors();

        amplitude_r_i = 0.0;
        for(int n=0; n<totalNeighbors; n++){ // kernel function
            amplitude_r_i += pAmp[model_neighbors[n]] * model_kernel_distances[n];
        } // pAmp is same length as input PDB model
        sum += amplitude_r_i;
        pDen[i] = amplitude_r_i;
    }


    //auto total_squared = totalHCPInUse*(totalHCPInUse-1)/2 + totalHCPInUse;
    /*
     * pre calculate squared amplitudes
     */
    float amp1;
    unsigned int indexer =0;
//    float * const pSq = squared_amplitudes.data();
    for(unsigned int i=0; i<totalHCPInUse; i++){
        amp1 = pDen[i];
        //ni = i*totalHCPInUse - (i*(i-1)/2);
        //std::cout << indexer << " " << ni << std::endl;
        for(unsigned int j=i; j<totalHCPInUse; j++){
//            squared_amplitudes[ni + (j-i)] = amp1*pDen[j];
//            pSq[indexer] = amp1*pDen[j];
            squared_amplitudes[indexer] = amp1*pDen[j];
            indexer++;
        }
    }

    return sum/(float)totalHCPInUse;
}

float DensityMapper::calculateAmplitudeAtQ(int q_index, int l_index, int m_index){
    // for each point in shell, calculate integral of
    float r_k, bessel_r_k2;
    unsigned int base = qvalues_size*rvalues.size();
    float sum_real = 0.0;
    float sum_imag = 0.0;

    int r=0;

    /*
     * Need to perform Hankel transform
     */
    unsigned int bessel_start = base*l_index + q_index*rvalues.size();

    for(auto & shell : shells){
        r_k = shell.getRadius();
        bessel_r_k2 = r_k*r_k*bessels[r + bessel_start];
        sum_real += shell.getPLMReal(l_index, m_index)*bessel_r_k2;
        sum_imag += shell.getPLMImag(l_index, m_index)*bessel_r_k2;
        r++;
    }

    return (float)(((sum_real*sum_real) + (sum_imag*sum_imag))*2.0/M_PI);
}





void DensityMapper::setBessels(const std::vector<float> &qvalues){

    qvalues_size = (int)(qvalues.size());

    logger("CALCULATING", "SPHERICAL BESSELS");

    unsigned int long total_qr_values = rvalues.size()*qvalues_size;
    bessel_size = (int) (total_qr_values*(lmax+1));

    std::vector<float> qr(total_qr_values);
    float * const pQR = qr.data();

    int index = 0;

    for(int q = 0; q < qvalues_size; q++) {
        const float *pq = &qvalues[q];
        for (float rvalue : rvalues) {
            pQR[ index ] = (*pq * rvalue);
            index++;
        }
    }

    bessels.resize(bessel_size);

    index = 0;
    for (unsigned int l=0; l <= lmax; l++){
        for (int i=0 ; i < total_qr_values; i++) {
            bessels[index] = boost::math::sph_bessel(l,pQR[i]);
            index++;
        }
    }

    logger("Total BESSELS", formatNumber((unsigned int)bessels.size()));
    logger("Total BESSELS (actual)", formatNumber((unsigned int)index));
    logger("Total BESSELS (calc)", formatNumber((unsigned int)bessel_size));
    logger("BESSELS (bytes)", formatNumber((unsigned int)(sizeof(std::vector<float>) + (sizeof(float) * bessels.size()))));
}


void DensityMapper::setAmplitudes(std::vector<float> & amps) {

    float * const ptrNewAmps = amps.data();

    for(unsigned int i = 0; i < amplitudes.size(); i++) {
        amplitudes[i] = ptrNewAmps[i];
    }
}



void DensityMapper::refineModel(int max_rounds, float topPercent, int models_per_round,
                                std::vector<Datum> & workingSet,
                                std::vector<Datum> & workingSetSmoothed ){
    /*
     *
     */
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(0.0,1.0);

    std::clock_t startTime;
    double runtime;

    auto topN = (unsigned int)(std::ceil((float)models_per_round*topPercent));

    auto total_amplitudes_from_lattice_points = lattice_points.size();
    std::vector<float> prior_model_amplitudes(total_amplitudes_from_lattice_points);
    amplitudes.resize(total_amplitudes_from_lattice_points);
    float * const pAmp = amplitudes.data();

    std::vector<Trial> topTrials;
    topTrials.reserve(topN);
    topTrials.resize(topN);

    const unsigned int last = topN-1;

    std::vector<float> qvalues;
    std::vector<float> sigmas_squared;

    for(auto & datum : workingSet){
            qvalues.push_back(datum.getQ());
            sigmas_squared.push_back(1.0f/(datum.getSigma()*datum.getSigma()));
    }

    const unsigned int total_data = workingSet.size();
    Datum * const pWorkingSet = workingSet.data();
    Datum * const pWorkingSetSmoothed = workingSetSmoothed.data();

    std::vector<float> i_calc(total_data);
    std::vector<float> trial_residuals(total_data);

    int total_sets_in_residuals=3;
    int total_residuals = total_data*total_sets_in_residuals; // only look at top 5 models
    std::vector<float> residuals(total_residuals);
    float * const pRes = residuals.data();

    float * const pICalc = i_calc.data();
    float * const pSigma = sigmas_squared.data();
    const LatticePoint * pLattice = lattice_points.data();

    float chi2, scale, ave, dw;

    //this->setBessels(qvalues);
    this->createHCPGrid(); //
    unsigned int long totalKept = keptHCPLatticePoints.size();
    std::vector<float> hcp_electron_densities(totalKept);

    unsigned int long size_of_dbf = (totalKept*(totalKept-1)/2 + totalKept)*workingSet.size();
//    debye_factors = new alignas(16) float[size_of_dbf];
    debye_factors = (float *)_aligned_malloc(sizeof(float)*size_of_dbf + 16, 16); // 16 byte aligned

    // use aligned memory for squared amplitudes and debye_factors
    //std::vector<float> squared_amplitudes(keptHCPLatticePoints.size()*(keptHCPLatticePoints.size()-1)/2 + keptHCPLatticePoints.size());
    auto * squared_amplitudes = (float *)_aligned_malloc(sizeof(float)*((totalKept*(totalKept-1)/2 + totalKept)), 16);
    this->setDebyeFactors(workingSet);
    //this->setYlms();

    logger("STARTING", "REFINEMENT");
    logger("Total Trials per Round", formatNumber((unsigned)models_per_round));
    logger("TopN", formatNumber((unsigned)topN));
    logger("10000 models at 0.01 per trial", formatNumber(10000*0.01/60, 2) + " minutes");
    logger("50000 models at 0.01 per trial", formatNumber(50000*0.01/60, 2) + " minutes");
    logger("UPDATE", "EACH ROUND ~ 2 to 10 minutes");

    int base_models_per_round = models_per_round;
    models_per_round *= 3;

    std::vector<float>durbinWatsons(max_rounds);
    std::vector<float>variations(max_rounds);
    int stretch = 4;
    std::vector<double>best(max_rounds); // average variation?
    double variation = 0;
    int round=0;

    std::ofstream logfile;
    logfile.open("logs_cemap.txt", std::ios_base::app);
    logfile << "# round best DW variations\n";
    logfile.close();
    char buffer[80];

    for(; round<max_rounds; round++){
        unsigned int topAdded=0;

        startTime = std::clock();
        // should be parallelized
        for(int model=1; model<models_per_round; model++){
            // populate amplitudes
            float * const pPriors = prior_model_amplitudes.data();
            for(int i=0; i < total_amplitudes_from_lattice_points; i++){
                pPriors[i] = pLattice[i].guessAmplitude(distribution(gen));
            }

//            for(auto & shell : shells){
//                shell.updateDensityCoefficients(centered_coordinates, prior_density_amplitudes, fwhm_sigma);
//            }

            //takes less than 100 usec
//            startTime = std::clock();
            ave = populateDensities(prior_model_amplitudes, hcp_electron_densities, squared_amplitudes);
//            runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//            logger("populateDensities TIME", formatNumber((float)runtime,8));

            /*
             * score the quality of the model
             * any densities isolated?
             * should be connected at average?
             *
             * look at all points above average, if any isolated, penalize by furthest distance
             *
             */

            // calculate intensities
            /*
             * by far the slowest part -
             */
//            startTime = std::clock();
            populateICalc(total_data, i_calc, squared_amplitudes);
//            runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//            logger("ICALC TIME", formatNumber((float)runtime,8));

//            startTime = std::clock();
            //takes less than 100 usec
            scale = calculateScaleFactor(total_data, pICalc, pSigma, pWorkingSetSmoothed);
            chi2 = getChiSquare(total_data, scale, pICalc, pSigma, pWorkingSet, pWorkingSetSmoothed, trial_residuals.data());
            //chi2 += alpha*(ave - 0.4f)*(ave - 0.4f);
            // update topN
            if (topAdded < topN){
                Trial * pTrial = &topTrials[topAdded];
                pTrial->value = chi2;
                pTrial->counter = ave;
                pTrial->scale = scale;
                pTrial->model_amplitudes.swap(prior_model_amplitudes);
                pTrial->residuals.swap(trial_residuals);

                prior_model_amplitudes.resize(total_amplitudes_from_lattice_points);
                trial_residuals.resize(total_data);

                topAdded++;
                if (topAdded == topN){
                    std::sort(topTrials.begin(), topTrials.begin()+topAdded);
                }
            } else {
                if (chi2 < topTrials[last].value){
                    /*
                     * replace last entry and sort
                     */
                    Trial * pTrial = &topTrials[last];
                    pTrial->value = chi2;
                    pTrial->counter = ave;
                    pTrial->scale = scale;
                    pTrial->model_amplitudes.swap(prior_model_amplitudes);
                    pTrial->residuals.swap(trial_residuals);

                    prior_model_amplitudes.resize(total_amplitudes_from_lattice_points);
                    trial_residuals.resize(total_data);
                    std::sort(topTrials.begin(), topTrials.end());
                }
            }
//            runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//            logger("Trial TIME", formatNumber((float)runtime,8));
        }
        //update probability model
        /*
         * for each bead in topN, count occurrences
         */
        for(auto & point : lattice_points){
            point.resetCounter();
        }

        for(auto & trial : topTrials){
            int lattice_index = 0;
            for(auto & assignedDensity : trial.model_amplitudes){
                lattice_points[lattice_index].addToCounter(assignedDensity);
                lattice_index++;
            }
        }

        // update probabilities
        for(auto & point : lattice_points){
            point.updateProbabilities();
        }

        // calculate DurbinWatson for the top5
        int indexer=0;
        for(int top5=0; top5<total_sets_in_residuals; top5++){
            // copy residuals into total_residuals
            auto & trial = topTrials[top5];
            for(auto & res : trial.residuals){
                pRes[indexer] = res;
                indexer++;
            }
        }

        dw = calculateDurbinWatson(total_residuals, residuals.data());
        durbinWatsons[round] = dw;
        best[round] = topTrials[0].value;

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
        std::cout << "==================== ROUND :: " << round << std::endl;
        logger("ROUND TIME", formatNumber((float)runtime/60.0f,2) + " minutes ");
        logger("CHI2 TOP BEST", formatNumber((float)topTrials[0].value, 4));
        logger("CHI2 2nd BEST", formatNumber((float)topTrials[1].value, 4));
        logger("CHI2 TOP LAST", formatNumber((float)topTrials[last].value, 4));

        float timeremaining = ((float)runtime*(float)(max_rounds - round)/60.0f);

        // update amplitudes with best
        float * pDens = topTrials[0].model_amplitudes.data();
        // update amplitudes with best
        for(int i=0; i<lattice_points.size(); i++){
            pAmp[i] = pDens[i];
        }

        this->createXPLORMap("best_"+std::to_string(round));
        writeLatticePoints("best_"+std::to_string(round));

        ave = populateDensities(amplitudes, hcp_electron_densities, squared_amplitudes);
        populateICalc(total_data, i_calc, squared_amplitudes);
        scale = topTrials[0].scale;
        writeICalc(total_data, scale, pICalc, pWorkingSet, "best_" + std::to_string(round));

        logger("DurbinWatson", formatNumber(dw, 5));
        logger("EST TIME REMAINING", formatNumber(timeremaining, 2) + " minutes");

        /*
         * determine maximum variation in last 5 best
         */
        if (round > stretch){ // check stopping criteria
            int start_at = round - (stretch-1);
            for(int k=start_at; k<= round; k++){
                float val1 = best[k];
                int l=k+1;
                for(; l<=round; l++){
                    float tempval = std::abs(val1 - best[l]);
                    if (tempval > variation){
                        variation = tempval;
                    }
                }
            }
            logger("Last 4 MAX DIFF", formatNumber(variation, 3));
            variations[round] = variation;
            if ( (variation < 0.01 && dw < 0.1) || round == (max_rounds-1)){
                for(int i=0; i<lattice_points.size(); i++){
                    pAmp[i] = 0;
                }

               // make the average of top 3
                for(int t=0; t<3; t++){
                    std::cout << "trial : " << t << " " << topTrials[t].value << std::endl;
                   // update amplitudes with best
                    for(int i=0; i<lattice_points.size(); i++){
                        pAmp[i] += pDens[i];
                    }
                }

                for(int i=0; i<lattice_points.size(); i++){
                    pAmp[i] /= 3.0f;
                }

                this->createXPLORMap("top3_averaged");
                writeLatticePoints("top3_averaged");
                break;
            }
            variation = 0;
        }

        std::snprintf(buffer, 80, "%8i %.3f %.3f %.3f\n", round, best[round], durbinWatsons[round], variations[round]);
        logfile.open("logs_cemap.txt", std::ios_base::app);
        logfile << buffer;
        logfile.close();

        std::cout << " " << std::endl;

        if (round > 3 && models_per_round >= base_models_per_round){
            models_per_round = base_models_per_round;
        }
    }

    // for plotting
//    std::string tempHeader = "# round best durbinwatson\n";
//    char buffer[80];
//    for(int r=0; r<round; r++){
//        std::snprintf(buffer, 80, "%8i %.3f %.3f %.3f\n", r, best[r], durbinWatsons[r], variations[r]);
//        tempHeader.append(buffer);
//    }
//
//    std::string logname = "logs_cemap.txt";
//    const char * outputFileName = logname.c_str();
//    FILE * pFile = fopen(outputFileName, "w");
//    fprintf(pFile, tempHeader.c_str());
//    fclose(pFile);

    _aligned_free(squared_amplitudes);
    _aligned_free(debye_factors);
}


void DensityMapper::printLatticePointsInfo(){

    for(auto & point : lattice_points){
        //point.printProbabilities();
        point.PDBLineForOutput(centered_coordinates[point.getIndex()]);
    }
}

/*
 * Use to make final density map based on amplitudes of model
 * delta_r => PI/(qmax * sampling_frequency)
 * cutoff => PI/qmax
 */
void DensityMapper::createDensityMapGrid() {

    // build enclosing box 2*kmax, spacing is delta_r
    int na;// = (int)std::ceil(2*kmax/delta_r) + 1;
    int nb;// = na;
    int nc;// = na;

    std::string tempHeader = headerParametersForXPLOR(na, nb, nc);

    // 0,0,0 index is (-kmax, -kmax, -kmax)
    // create grid
    for(int a=0; a<na; a++){
        float a_deltar = a*delta_r - kmax;
        for(int b=0; b<nb; b++){
            float b_deltar = b*delta_r - kmax;
            for(int c=0; c<nc; c++){
                grid_for_map.emplace_back(vector3(a_deltar, b_deltar, c*delta_r - kmax));
            }
        }
    }

    // which grid points overlap with model
    int total_in_grid = grid_for_map.size();
    unsigned int long total_in_coords = centered_coordinates.size();
    useIt.resize(total_in_grid);

    for(int i=0; i<total_in_grid; i++){
        useIt[i] = false;
    }

    for(unsigned int long i=0; i<total_in_grid; i++){

        vector3 * grid_vec = &grid_for_map[i];

        for(unsigned int long j_model=0; j_model<total_in_coords; j_model++){
            vector3 * coor_vec = &centered_coordinates[j_model];
            if ((*(grid_vec) - *(coor_vec)).length() < cutoff){
                auto it = neighbors.find(i);
                if (it == neighbors.end()){
                    auto temp = neighbors.emplace(i, std::vector<unsigned int long>() );
                    (*temp.first).second.push_back(j_model);
                    useIt[i] = true;
                } else {
                    it->second.push_back(j_model);
                }
            }
        }
    }
}


std::string DensityMapper::headerParametersForXPLOR(int & na, int & nb, int & nc){

    float grid_spacing = delta_r;
//    std::cout << " BOUNDING BOX "<< std::endl;
//    std::cout << "   AXIS   MIN      MAX     LENGTH"  << std::endl;
    float cminx = -kmax;
    float cmaxx = kmax;
    float cminy = -kmax;
    float cmaxy = kmax;
    float cminz = -kmax;
    float cmaxz = kmax;
    float a_side = cmaxx - cminx;
    float b_side = cmaxy - cminy;
    float c_side = cmaxz - cminz;

//    printf("  => X %8.4f %8.4f %8.4f \n", cminx, cmaxx, a_side);
//    printf("  => Y %8.4f %8.4f %8.4f \n", cminy, cmaxy, b_side);
//    printf("  => Z %8.4f %8.4f %8.4f \n", cminz, cmaxz, c_side);

    auto startingNA = (int)std::round(cminx/grid_spacing);
    auto startingNB = (int)std::round(cminy/grid_spacing);
    auto startingNC = (int)std::round(cminz/grid_spacing);

    auto stoppingNA = (int)std::round(cmaxx/grid_spacing);
    auto stoppingNB = (int)std::round(cmaxy/grid_spacing);
    auto stoppingNC = (int)std::round(cmaxz/grid_spacing);

    na = std::abs(startingNA) + std::abs(stoppingNA) + 1;
    nb = std::abs(startingNB) + std::abs(stoppingNB) + 1;
    nc = std::abs(startingNC) + std::abs(stoppingNC) + 1;

    std::string tempHeader = "\n";
    tempHeader += "        4 !NTITLE\n";
    tempHeader += "REMARK 265 EXPERIMENTAL DETAILS\n";
    tempHeader += "REMARK 265\n";
    tempHeader += "REMARK 265 EXPERIMENT TYPE : X-RAY SOLUTION SCATTERING\n";
    tempHeader += "REMARK 265 DATA ACQUISITION\n";

    char buffer[80];
    std::snprintf(buffer, 80, "%8i%8i%8i%8i%8i%8i%8i%8i%8i\n", na, startingNA, stoppingNA, nb, startingNB, stoppingNB, nc, startingNC, stoppingNC);
    tempHeader.append(buffer);

    std::snprintf(buffer, 80, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E \n",a_side, b_side, c_side, 90.0, 90.0, 90.0);

    tempHeader.append(buffer);
    // std::cout << tempHeader << std::endl;
    // write the matrix
    tempHeader += "ZYX\n";

    return tempHeader;
}

std::string DensityMapper::headerParametersForXPLORFlipped(int & na, int & nb, int & nc){

    float grid_spacing = delta_r;

    float cminx = -kmax;
    float cmaxx = kmax;
    float cminy = -kmax;
    float cmaxy = kmax;
    float cminz = -kmax;
    float cmaxz = kmax;
    float a_side = cmaxx - cminx;
    float b_side = cmaxy - cminy;
    float c_side = cmaxz - cminz;

    auto startingNA = (int)std::round(cminx/grid_spacing);
    auto startingNB = (int)std::round(cminy/grid_spacing);
    auto startingNC = (int)std::round(cminz/grid_spacing);

    auto stoppingNA = (int)std::round(cmaxx/grid_spacing);
    auto stoppingNB = (int)std::round(cmaxy/grid_spacing);
    auto stoppingNC = (int)std::round(cmaxz/grid_spacing);

    na = std::abs(startingNA) + std::abs(stoppingNA) + 1;
    nb = std::abs(startingNB) + std::abs(stoppingNB) + 1;
    nc = std::abs(startingNC) + std::abs(stoppingNC) + 1;

    std::string tempHeader = "\n";
    tempHeader += "        4 !NTITLE\n";
    tempHeader += "REMARK 265 EXPERIMENTAL DETAILS\n";
    tempHeader += "REMARK 265\n";
    tempHeader += "REMARK 265 EXPERIMENT TYPE : X-RAY SOLUTION SCATTERING\n";
    tempHeader += "REMARK 265 DATA ACQUISITION\n";

    char buffer[80];
    std::snprintf(buffer, 80, "%8i%8i%8i%8i%8i%8i%8i%8i%8i\n", na, -stoppingNA, -startingNA, nb, startingNB, stoppingNB, nc, startingNC, stoppingNC);
    tempHeader.append(buffer);

    std::snprintf(buffer, 80, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E \n",a_side, b_side, c_side, 90.0, 90.0, 90.0);

    tempHeader.append(buffer);
    // write the matrix
    tempHeader += "ZYX\n";

    return tempHeader;
}


void DensityMapper::createXPLORMap(std::string name){

    int na, nb, nc;

    float grid_spacing = delta_r;

    //float inv_bandwidth = 1.0f/(cutoff/2.355f); // should be same as fwhm_sigma?
    float inv_bandwidth = 1.0f/fwhm_sigma;
    float inv_bandwidth_squared = inv_bandwidth*inv_bandwidth;
    float inv_bandwidth_2PI = 1.0f/sqrtf(2.0f*M_PI)*inv_bandwidth;
    float limit = 2.0f*cutoff;
    float mapSum = 0.0f, mapSumSquared = 0.0f;
    float mapCount = 0.0f;

    const vector3 * pVec = centered_coordinates.data();
    const float * pAmp = amplitudes.data();

    std::string tempHeader = headerParametersForXPLOR(na, nb, nc);

    char buffer[80];

    std::vector<float> averages;
    std::vector<vector3> coords_averages;
    std::map<std::string, float > kernel_mapping;

    for(int z=0; z<nc; z++){

        float zsection = -kmax + z*grid_spacing;

        std::snprintf(buffer, 80, "%8i\n", z);
        tempHeader.append(buffer);
        unsigned int stringIndex = 0;

        for(int y=0; y<nb; y++){

            float ysection = -kmax + y*grid_spacing;

            for(int x=0; x<na; x++){

                std::string key = std::to_string(x)+"-"+std::to_string(y)+"-"+std::to_string(z);

                vector3 locale(-kmax+x*grid_spacing, ysection, zsection);

                float kernelSum = 0.0f;

                for(unsigned int c=0; c<total_centered_coordinates; c++){ // calculate distances to vectors in compiled model
                    /*
                     * centered_coordinates -> model
                     * lattice_points -> amplitudes assigned to each point in model
                     */
                    vector3 vec1 = *(pVec+c) - locale;
                    float length = vec1.length();

                    if (length < limit){
                        kernelSum += *(pAmp + c)*exp(-length*length*0.5*inv_bandwidth_squared);
                    }
                    // calculate amplitude based on
                }

                kernelSum *= inv_bandwidth_2PI;
                averages.push_back(kernelSum);
                coords_averages.emplace_back(vector3(locale));

                std::snprintf(buffer, 80, "%12.5E", kernelSum);
                tempHeader.append(buffer);
                kernel_mapping.emplace(key, kernelSum);

                stringIndex += 1;

                if (stringIndex%6 == 0 ){
                    tempHeader += "\n";
                }

                mapSum += kernelSum;
                mapSumSquared += kernelSum*kernelSum;
                mapCount += 1.0;
            }
        }

        if (stringIndex % 6 != 0){
            tempHeader += "\n";
        }
    }

    std::snprintf(buffer, 80, "%8i\n",-9999);
    tempHeader.append(buffer);

    float ave = mapSum/mapCount;
    std::snprintf(buffer, 80, "%12.4E %12.4E\n", ave, std::sqrt(mapSumSquared/mapCount - ave*ave));
    tempHeader.append(buffer);

    // write to file
    std::string map = name + "_map.xplor";
    const char * outputFileName = map.c_str();

    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);

    // write out PDB cooordinates
    tempHeader = "";

    float average = (mapSumSquared/mapCount);
    float stdev = std::sqrt(mapSumSquared/mapCount - ave*ave);
    float aboveit = average + 2.5f*stdev;

    int index = 1;
    for(int i=0; i<averages.size(); i++){
        float value = averages[i];
        if(value > aboveit){ // print
            vector3 & vec = coords_averages[i];
            std::snprintf(buffer, 80, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00%5.2f\n", "ATOM", 1," CA ", "ALA", "A", index, vec.x, vec.y, vec.z, value );
            tempHeader.append(buffer);
            index++;
        }
    }

    tempHeader.append("END\n");
    std::string out = name+"_above_average.pdb";
    outputFileName = out.c_str();

    pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);

    /*
     * calculate flipped map
     * flip the map x -> -x
     * reverse the order for each section
     */
    tempHeader = headerParametersForXPLORFlipped(na, nb, nc);
    out = name + "_flipped";

    for(int z=0; z<nc; z++){

        std::snprintf(buffer, 80, "%8i\n", z);
        tempHeader.append(buffer);
        unsigned int stringIndex = 0;

        for(int y=0; y<nb; y++){
            for(int x=(na-1); x>=0; x--){ // reverse order

                std::string key = std::to_string(x)+"-"+std::to_string(y)+"-"+std::to_string(z);

                auto allIt = kernel_mapping.find(key);
                /*
                 * write to file
                 */
                std::snprintf(buffer, 80, "%12.5E", (*allIt).second);
                tempHeader.append(buffer);
                stringIndex += 1;

                if (stringIndex%6 == 0 ){
                    tempHeader += "\n";
                }
            }
        }

        if (stringIndex % 6 != 0){
            tempHeader += "\n";
        }
    }

    map = out + "_map.xplor";
    outputFileName = map.c_str();

    pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);

    // write out flipped PDB coordinates
    tempHeader = "";

    index = 1;
    for(int i=0; i<averages.size(); i++){
        float value = averages[i];
        if(value > aboveit){ // print
            vector3 & vec = coords_averages[i];
            std::snprintf(buffer, 80, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00%5.2f\n", "ATOM", 1," CA ", "ALA", "A", index, -vec.x, vec.y, vec.z, value );
            tempHeader.append(buffer);
            index++;
        }
    }
    tempHeader.append("END\n");
    out = name+"flipped_above_average.pdb";
    outputFileName = out.c_str();

    pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);
}


void DensityMapper::writeLatticePoints(std::string name){

    float sum = 0.0f;

    for(auto & lat : lattice_points){
        sum += lat.getMostProbableAmplitude();
    }
    float avg = sum/(float)lattice_points.size();

    char buffer[80];
    std::string tempHeader = "";
    std::snprintf(buffer, 80, "REMARK AVG %.3f\n", avg);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK Total Points in Base %5i\n", lattice_points.size());
    tempHeader.append(buffer);
    int index =1;
    /*
     * centered_coordinates holds coordiniates lattice_point
     */
    for(auto & lat : lattice_points){
        if (lat.getMostProbableAmplitude() > avg){
            vector3 & vec = centered_coordinates[lat.getIndex()];
            std::snprintf(buffer, 80, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00%5.2f\n", "ATOM", 1," CA ", "ALA", "A", index, vec.x, vec.y, vec.z, lat.getMostProbableAmplitude() );
            tempHeader.append(buffer);
            index++;
        }
    }

    std::string out = name + "_base.pdb";
    const char * outputFileName = out.c_str();

    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);

}

void DensityMapper::createHCPGrid(){
    // read in pdb file

    modelDensityHCP = PointSetModel((2*kmax), delta_r); // grid is more closely spaced than input model

    const float beadradius_limit = (delta_r*std::sqrt(3.0f/2.0f));

    float inv_sigma = 1.0f/(2.0f*fwhm_sigma*fwhm_sigma);
    float inv_sqrt = 1.0f/(fwhm_sigma*sqrtf(2.0f*M_PI));

    auto totalBeads = modelDensityHCP.getTotalNumberOfBeadsInUniverse();

    vector3 * tempvec;

    for(unsigned int i=0; i < totalBeads; i++){
        auto bead = modelDensityHCP.getBead(i);
        auto bead_vec = bead->getVec();

        for (unsigned int a=0; a < total_centered_coordinates; a++){

            tempvec = &centered_coordinates[a];

            if ( (bead_vec - *tempvec).length() < beadradius_limit){ // if bead overlaps with atom, keep it
                keptHCPLatticePoints.push_back(i);
                break;
            }
        }
    }

    // kepHCPLatticePoints should contain only unique indices
    //rvalues.clear();
    //neighborhoods.clear();

    double phi_calc;
    float ncoount;
    std::string tempHeader;
    char buffer[80];

    // for each bead, get all atoms that overlap
    int indexer=1, ind=0;
    const char *colour[6] = { "A", "B", "C", "D", "E", "F" };

    for(auto & index : keptHCPLatticePoints){

        // make new neighborhood for each kept point
        neighborhoods.emplace_back(Neighbors(index));
        // get neighborhood
        Neighbors & lastOne = neighborhoods[neighborhoods.size() - 1];

        auto bead = modelDensityHCP.getBead(index);
        auto bead_vec = bead->getVec();

        rvalues.push_back(bead_vec.length());

        // for SHE (spherical coordinates)
        phi_calc = (bead_vec.x == 0 || bead_vec.x == 0.0) ? 0.0 : atan2(bead_vec.y,bead_vec.x);
        phis.push_back(phi_calc);
        cos_thetas.push_back( cos(acos(bead_vec.z/bead_vec.length())) );

        /*
         * calculate how far away each selected lattice point is to neighboring point in MODEL
         * if acceptable, calculate Gaussian kernel and put into neighborhood
         */
        ncoount = 0.0;
        for (unsigned int a=0; a < total_centered_coordinates; a++){
            tempvec = &centered_coordinates[a];
            float dis = (bead_vec - *tempvec).length();
            if (dis < 1.2*cutoff){ // if bead overlaps with atom, keep it -> cutoff = PI/qmax
                /*
                 * inv_sigma = 1.0f/(2.0f*fwhm_sigma*fwhm_sigma);
                 * fwhm_sigma = 1.5f*cutoff/2.355f;
                 */
                lastOne.add_neighbor(a, inv_sqrt*expf(-dis*dis*inv_sigma));
                ncoount+=1.0;
            }
        }
        //std::cout << lastOne.getHCPLatticeIndex() << " " << lastOne.getTotalNeighbors() << std::endl;

        std::snprintf(buffer, 80, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00%5.2f\n", "ATOM", 1," CA ", "ALA", colour[ind], indexer, bead_vec.x, bead_vec.y, bead_vec.z, ncoount );
        tempHeader.append(buffer);
        indexer++;

        if (indexer%9999==0){
            indexer=1;
            ind += 1;
        }
    }

    logger("Total KEPT HCP Lattice", formatNumber((unsigned int)keptHCPLatticePoints.size()));

    std::string name = "hcp_grid_model.pdb";
    const char * outputFileName = name.c_str();
    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);
}



void DensityMapper::tester(std::string pdbfile, std::vector<Datum> & workingSet){

    // read in pdb file
    PDBModel pdbModel(pdbfile, true, true); // coordinates are centered
    float searchSpace = pdbModel.getDmax() + pdbModel.getDmax()*0.2f;
    const float beadradius = cutoff/8.0f;

    float volume = 4.0f/3.0f*M_PI*beadradius*beadradius*beadradius;

    PointSetModel model(searchSpace, beadradius);
    unsigned int totalAtoms = pdbModel.getTotalCoordinates();

    const float beadradius_limit = (beadradius*std::sqrt(3.0f/2.0f));

    auto totalBeads = model.getTotalNumberOfBeadsInUniverse();

    std::vector<unsigned int> indicesToCheck(totalAtoms);

    unsigned int * const ptr = (totalAtoms != 0) ? indicesToCheck.data() : nullptr;
    for(unsigned int i = 0; i < totalAtoms; i++) {
        ptr[i] = i;
    }

    vector3 tempvec;
    std::vector<unsigned int> keptbeads;
    for(unsigned int i=0; i<totalBeads; i++){
        auto bead = model.getBead(i);
        auto bead_vec = bead->getVec();

        for (unsigned int a=0; a < totalAtoms; a++){
            tempvec.x = *(pdbModel.getCenteredXVec() + ptr[a]);
            tempvec.y = *(pdbModel.getCenteredYVec() + ptr[a]);
            tempvec.z = *(pdbModel.getCenteredZVec() + ptr[a]);

            if ((bead_vec - tempvec).length() < beadradius_limit){ // if bead overlaps with atom, keep it
                keptbeads.push_back(i);
                break;
            }
        }
    }

    std::map<unsigned int, std::vector< int >> atom_types_per_bead;
    std::set<unsigned int> has_CA;
    // for each bead, get all atoms that overlap
    for(auto & index : keptbeads){
        auto bead = model.getBead(index);
        auto bead_vec = bead->getVec();
        for (unsigned int a=0; a < totalAtoms; a++){
            tempvec.x = *(pdbModel.getCenteredXVec() + ptr[a]);
            tempvec.y = *(pdbModel.getCenteredYVec() + ptr[a]);
            tempvec.z = *(pdbModel.getCenteredZVec() + ptr[a]);

            if ((bead_vec - tempvec).length() < beadradius_limit){ // if bead overlaps with atom, keep it
                auto it = atom_types_per_bead.find(index);
                if (it == atom_types_per_bead.end()){
                    auto temp = atom_types_per_bead.emplace(index, std::vector< int >() );
                    (*temp.first).second.push_back(pdbModel.getAtomicNumberByIndex(a));
                } else {
                    it->second.push_back(pdbModel.getAtomicNumberByIndex(a));
                }

                std::string tempAtom = pdbModel.getAtomTypeByIndex(a);
                trimWhiteSpace(tempAtom);
                if (tempAtom == "CA"){
                    has_CA.insert(index);
                }
            }
        }
    }


    std::vector<float> total_electrons_per_bead;
    for(auto & pair : atom_types_per_bead){
        auto numbers = pair.second;
        float sum=0.0f;
        for(auto & num : numbers){
            sum += num;
        }
        total_electrons_per_bead.push_back(sum/volume);
        //std::cout << index << " => " << " average " << (sum/numbers.size()) << " " << sum << " d: " << (sum/volume) << std::endl;
    }

    // set amplitudes for all
    std::string tempHeader = "";
    char buffer[80];

    int index = 1;
    for(auto & bd : keptbeads){
        auto bead = model.getBead(bd);
        auto vec = bead->getVec();
        float occ = total_electrons_per_bead[index-1];
        amplitudes[index-1] = occ;

        std::snprintf(buffer, 80, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00%5.2f\n", "ATOM", 1," CA ", "ALA", "A", index, vec.x, vec.y, vec.z, occ );
        tempHeader.append(buffer);
        index++;
    }

    std::cout << " AMPS " << amplitudes.size() << " kept " << keptbeads.size() << " " << centered_coordinates.size() << std::endl;

    std::string out = "centered_SHE_test_model.pdb";
    const char * outputFileName  = out.c_str();

    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);

    for(auto & shell : shells){
        shell.updateDensityCoefficients(centered_coordinates, amplitudes, fwhm_sigma);
    }

    int total_data = workingSet.size();
    std::vector<float> i_calc(total_data);
    float * const pICalc = i_calc.data();

    for(int i=0; i<total_data; i++){
        pICalc[i] = this->calculateIntensityAtQ(i);
    }


    int totalKeptBeads = keptbeads.size();
    float qr;
    for(int i=0; i<total_data; i++){

        float qvalue = workingSet[i].getQ();
        float intensity = 0.0;

        for(int j=0; j<totalKeptBeads; j++){
            auto bead = model.getBead(keptbeads[j]);
            auto vec = bead->getVec();
            for(int k=(j+1); k<totalKeptBeads; k++){
                qr = qvalue*((vec - model.getBead(keptbeads[k])->getVec()).length());
                intensity += 2*sinf(qr)/qr;
            }
        }

        intensity += totalKeptBeads;

        std::cout << workingSet[i].getQ() << " " << pICalc[i] << " " << workingSet[i].getI() << " " << intensity << std::endl;
    }

    // make a mask
    // for any bead overlapping with CA, give amp 1, not 0.3
    // use mask as input - see if we recover bead-model
}


void DensityMapper::trimWhiteSpace(std::string &text) {
    text.erase(text.begin(), std::find_if(text.begin(), text.end(), [](int ch) {
        return !std::isspace(ch);
    }));

    text.erase(std::find_if(text.rbegin(), text.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), text.end());
}

/*
 * Calculate Debye factors for each point kept in the HCP lattice
 *
 * keptHCPLatticePoints must be populated before use
 *
 */
void DensityMapper::setDebyeFactors(std::vector<Datum> &workingSet) {

    Bead * bead2;
    float qr;

    //int totalKept = neighborhoods.size();
    int totalKept = keptHCPLatticePoints.size();
    // total keptHCPLatticePoints is same size as neighborhoods

    //debye_factors;
    std::vector<float> distances;

    // precompute all distances
    for (int i=0; i<totalKept; i++){
        auto bead = modelDensityHCP.getBead(keptHCPLatticePoints[i]);
        auto bead_vec = bead->getVec();
        distances.push_back(0.0f);
        int next = i+1;
        for (int j=next; j<totalKept; j++){
            bead2 = modelDensityHCP.getBead(keptHCPLatticePoints[j]);
            distances.push_back((bead_vec - bead2->getVec()).length());
        }
    }
    //std::cout << " Total Distances " << (totalKept*(totalKept-1)/2 + totalKept) << " " << distances.size() << std::endl;
    // populate debye factors


//    unsigned int long size_of_dbf = (totalKept*(totalKept-1)/2 + totalKept)*workingSet.size();
//    debye_factors = (float *)_aligned_malloc(sizeof(float)*size_of_dbf, 16); // 16 byte aligned
//    debye_factors = new float[size_of_dbf]; // 16 byte aligned

    unsigned int long counter = 0;
    for(auto & data : workingSet){

        float q_value = data.getQ();
        unsigned int indexer = 0;
        for (int i=0; i<totalKept; i++){
            debye_factors[counter] = 1.0f;
            counter++;
//            debye_factors.push_back(1.0f);
//            int next = i*(i-1)/2;
            indexer++;
            for (int j=(i+1); j<totalKept; j++){
//                next++;
//                qr = q_value*distances[next];
                qr = q_value*distances[indexer];
                indexer++;
                //debye_factors.push_back(sinf(qr)/qr);
                debye_factors[counter] = sinf(qr)/qr;
                counter++;
            }
        }
    }
}

void DensityMapper::populateICalc(unsigned int total_q, std::vector<float> & iCalc, float * squared_amplitudes) {

    unsigned int totalHCPInUse = keptHCPLatticePoints.size(); // shared memory
//    const float * const pAmp = &squared_amplitudes.front();
//    const float * const pDebye = &debye_factors[0];//.front(); // shared memory
//    const simde_float32 * const pDebye = &debye_factors.front(); // shared memory

    float intensity_r_at_zero=0.0f;
    for (int i = 0; i < totalHCPInUse; i++) { // calculate diagonal term of Debye matrix
        intensity_r_at_zero += squared_amplitudes[ i*totalHCPInUse - (i*(i-1)/2)];// pDebye is 1 for self-amplitude squared
    }

    const unsigned int total_width = totalHCPInUse*(totalHCPInUse-1)/2+totalHCPInUse;

    float intensity2;
    unsigned int locale;

    float * pICalc = iCalc.data();
    unsigned int window = 8, next_i; // window is 8 for _m128 and 16 for _m256
    auto end_N = (unsigned int)std::floor(total_width/window);
    unsigned int start_of_tail = end_N*window;

    simde__m128 mmSum1, mmSum2, v0, v1, v2, v3, s01, s23;
//    simde__m256 mmSum1, mmSum2; // vector of 8 floating point numbers

    float icalc_sse[4]; // size is 4 for m128 and 8 for m256

    for(unsigned int q_index = 0; q_index < total_q; q_index++) {

        next_i  = q_index * total_width;
        mmSum1 = simde_mm_setzero_ps();
        mmSum2 = simde_mm_setzero_ps();

        for (unsigned int i=0; i < start_of_tail; i+=window ){ // SSE intrinsic here would multiply and then sum
            next_i += window;
            //_m128 is four floating point numbers (4x32 = 128 bits) => 16 bytes total
            // one byte is equivalent to eight bits
            v0 = simde_mm_load_ps(&squared_amplitudes[i]); // loading 4 floats at a time or 16 bytes total
            v1 = simde_mm_load_ps(&debye_factors[next_i]);
//            s01 = simde_mm_mul_ps(v0,v1);
            mmSum1 = simde_mm_add_ps(mmSum1, simde_mm_mul_ps(v0,v1));

            v2 = simde_mm_load_ps(&squared_amplitudes[i + 4]); // loading next 4 floats for total 8 per cycle
            v3 = simde_mm_load_ps(&debye_factors[next_i + 4]);
//            s23 = simde_mm_mul_ps(v2,v3);
            mmSum2 = simde_mm_add_ps(mmSum2, simde_mm_mul_ps(v2,v3));
            // unroll 2-times in windows of 8 long for _m256
//            simde__m256 v0 = simde_mm256_loadu_ps(&squared_amplitudes[i]);
//            simde__m256 v1 = simde_mm256_loadu_ps(&debye_factors[next_i]);
//            simde__m256 s01 = simde_mm256_mul_ps(v0, v1);
//            mmSum1 = simde_mm256_add_ps(mmSum1, s01);
//
//            simde__m256 v2 = simde_mm256_loadu_ps(&squared_amplitudes[i + 8]);
//            simde__m256 v3 = simde_mm256_loadu_ps(&debye_factors[next_i + 8]);
//            simde__m256 s23 = simde_mm256_mul_ps(v2, v3);
//            mmSum1 = simde_mm256_add_ps(mmSum1, s23);
        }
        mmSum1 = simde_mm_hadd_ps(mmSum1, mmSum2);
        simde_mm_storeu_ps(icalc_sse, mmSum1);
//        simde_mm_store_ps(icalc_sse, mmSum1);
//        simde_mm256_store_ps(icalc_sse, mmSum1);
        intensity2=icalc_sse[0] + icalc_sse[1] + icalc_sse[2] + icalc_sse[3];// + icalc_sse[4] + icalc_sse[5] + icalc_sse[6] + icalc_sse[7];

        // add remaining bits since we can't do in full width of 8
        for (unsigned int i=start_of_tail; i < total_width; i++ ){ // SSE intrinsic here would multiple and then sum
            intensity2 += squared_amplitudes[i] * debye_factors[locale + i];
        }

        *(pICalc) = intensity_r_at_zero + 2*(intensity2 - intensity_r_at_zero);
        ++pICalc;
    }
}

float DensityMapper::calculateScaleFactor(unsigned int total, float * const pICalc, float * const pSigma, Datum * const pWorkingSet) {

    float scale_top = 0.0;
    float scale_bottom = 0.0;
    float sigma_value, i_calc_value;

    for(int i=0; i<total; i++){
        i_calc_value = pICalc[i];
        sigma_value = pSigma[i]*i_calc_value;
        scale_top += sigma_value*pWorkingSet[i].getI();
        scale_bottom += sigma_value*i_calc_value;
    }

    scale_top /= scale_bottom;

    return scale_top;
}

/*
 * calculate chi2 score using smoothed dataset and Durbin watson using the raw working set.
 */
float DensityMapper::getChiSquare(unsigned int total, float scale, float * const pICalc, float * const pSigma,
                                  Datum * const pWorkingSet,
                                  Datum * const pWorkingSetSmoothed,
                                  float * const res){
    float diff_value;
    float chi2 = 0.0;
    float icalc;

    for(unsigned int i=0; i <total; i++){
        icalc = pICalc[i];
        diff_value = pWorkingSetSmoothed[i].getI() - scale*icalc;
//        diff_value = pWorkingSet[i].getI() - scale*icalc;
//        res[i] = diff_value;
        diff_value *= diff_value;
        diff_value *= pSigma[i];
        chi2 += diff_value;

        // use the workset residuals to calculate durbin watson
        res[i] = pWorkingSet[i].getI() - scale*icalc;
    }
    return chi2/(float)total;
}


void DensityMapper::writeICalc(unsigned int total, float scale, float * const pICalc, Datum * const pWorkingSet, std::string name){

    std::string tempHeader = "";
    char buffer[80];

    for(int i=0; i <total; i++){
        std::snprintf(buffer, 80, "%5i %.5E %.5E %.5E \n", i, pWorkingSet[i].getQ(), pWorkingSet[i].getI(), (scale*pICalc[i]));
        tempHeader.append(buffer);
    }

    std::string mname = name + "_fit.txt";
    const char * outputFileName = mname.c_str();
    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);
}


float DensityMapper::calculateDurbinWatson(unsigned int total, float * const residuals){

    float top=0.0f;
    float e_t, diff;

    e_t = residuals[0];
    float bottom = e_t*e_t;

    for(unsigned int i=1; i <total; i++){
        e_t = residuals[i];
        diff = (e_t - residuals[i-1]);
        top += diff*diff;
        bottom += e_t*e_t;
    }

    diff = 2.0f - top/bottom;
    return diff*diff;
}

int DensityMapper::getTotal_centered_coordinates() const {
    return total_centered_coordinates;
}



//void DensityMapper::openMP() {
//
//    const double N1 = 1.0 / 4294967291;
//    double res = 0.0;
//    double wall0 = omp_get_wtime();
//
//    omp_set_num_threads(1);
//#pragma omp parallel reduction(+ : res)
//    {
//        double tmp = 0.0;
//#pragma omp for
//        for (long i = 0; i < 4294967291; i++)
//        {
//            const double x = (i + 0.5) * N1;
//            tmp += 4.0 / (x * x + 1);
//        }
//        res += tmp;
//    }
//
//    double wall1 = omp_get_wtime();
//    std::cout << "THREADS " << omp_get_num_threads() << "elapsed time= " << wall1 - wall0 << " s :: " << res*N1 << std::endl;
//
//    for(int i = 0; i<141; i++){
//        std::cout << i << " " << i%21 << std::endl;
//    }
//}



