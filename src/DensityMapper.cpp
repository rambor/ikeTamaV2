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
    float max_r = 0.0f;
    for(unsigned int i = 0; i<totalCoords; i++) {
        float pX = tempmodel.getCenteredXVec()[i];
        float pY = tempmodel.getCenteredYVec()[i];
        float pZ = tempmodel.getCenteredZVec()[i];
        centered_coordinates.emplace_back(vector3(pX, pY, pZ));
        lattice_points.emplace_back(LatticePoint(i, 4));
        float value = pX*pX + pY*pY + pZ*pZ;
        if (value > max_r){
            max_r = value;
        }
    }

    this->setLatticePointNeighbors();

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

//    fwhm_sigma = cutoff*sqrtf(3.0f/8.0f); // => 0.61 seems to converge quickly
    float frac = 1.1f;//1.2*sqrtf(3.0f/8.0f);
    //fwhm_sigma = (1.5f-0.71f)*cutoff/2.355f;
    //fwhm_sigma = cutoff/2.355f; // 0.5x causes steep cut-off, density looks very isolated
    fwhm_sigma = sqrtf(-(0.5*0.5*cutoff*cutoff)*0.5/logf(0.01));
    /*
     * increasing fwhm_sigma_factor increases correlations between lattice points and makes it smoother
     * frac    fwhm_sigma_factor
     * 0.37     0.5        too bumpy
     * 0.37     0.612      better connectivity but bumpy, lots of features like connected tubes of density
     * 0.37     0.67       not converging but much smoother
     * 0.5      0.67       works quite well on the average
     */
//    within_limit =frac*cutoff; // 0.8 worked well
//    upper_limit = 1.2*cutoff;
//cut-off is PI/q_max (bin-width)
    within_limit = frac*cutoff; // 0.8 worked well
    upper_limit = frac*cutoff;

    // sharper the cutoff which is difference between 1.0f and frac, the more bumpy it is
    /*
     * cutoff*sqrtf(3.0f/8.0f) + 0.25*delta_r seems to converge quickly, if too large convergence is longer and
     * model is smoother.
     *
     * if fwhm is too small, like 0.2*cutoff, each point in Nyquist grid acts independent
     */

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

    for (unsigned int i=0; i< total_centered_coordinates; i++){
        neighboringHCPPointsOfModelLattice.insert(std::make_pair(i, std::set<unsigned int>()));
    }

    // use PointSetModel to make neighborhood for each point in centered coordinates
    inputBaseModel = PointSetModel(centered_coordinates, 2.1f*sqrtf(max_r), dmin_supremum*0.5);
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

//    createSphericalFibonacciLattice();
//    setModel();
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
void DensityMapper::updateDensities(std::vector<unsigned int> & indicesInUse, std::vector<float> & amplitudesT, std::vector<float> & densities_at_HCP, float * hcp_squared_amplitudes){
    unsigned int totalHCPInUse = keptHCPSamplingPoints.size();

    const float * const pAmp = amplitudesT.data(); // base amplitude for input MODEL
    float * const pDen = densities_at_HCP.data();

    Neighbors * const  pNeighborhood = neighborhoods.data();

    unsigned int totalNeighbors;
    float amplitude_r_i;
    float sum = 0.0;

    // make list of indices in keptHCPSamplingPoints that need to be update
    std::set<unsigned int> samplingPointsToUpdate;
    for(auto & ind : indicesInUse){
        auto iter = neighboringHCPPointsOfModelLattice.find(ind);

        if (iter != neighboringHCPPointsOfModelLattice.end()){
            for(auto useMe : iter->second){
                samplingPointsToUpdate.insert(useMe);
            }
        } else {
            std::cout << " not found ! " << ind << std::endl;
        }
    }

    unsigned int totalSamplingPoints = samplingPointsToUpdate.size();
    std::vector<unsigned int> indicesInpDen(totalSamplingPoints);

    int ii=0;
    for(auto & ind : samplingPointsToUpdate){
        // ind is the index with respect to total sampling lattice array
        // find neighborhood that corresponds to the index
        auto it = find_if(neighborhoods.begin(), neighborhoods.end(), [&ind](Neighbors & obj) {return obj.getHCPLatticeIndex() == ind;});
        auto index = std::distance(neighborhoods.begin(), it);
        indicesInpDen[ii] = index;
        ii++;

        const float * model_kernel_distances = it->getModelDistances();
        const unsigned int * model_neighbors = it->getModelIndices();

        totalNeighbors = it->getTotalNeighbors();

        amplitude_r_i = 0.0;
        for(int n=0; n<totalNeighbors; n++){ // kernel function
            amplitude_r_i += pAmp[model_neighbors[n]] * model_kernel_distances[n];
        } // pAmp is same length as input PDB model
        pDen[index] = amplitude_r_i/(float)totalNeighbors;
    }

    /*
     * update pre calculate squared amplitudes
     */
    float amp1;
    for(auto & ind : indicesInpDen){
        amp1 = pDen[ind];
        // down the column
        for(unsigned int j=0; j<ind; j++){
            unsigned int indexer = j*totalHCPInUse - (j*(j-1)/2) + (ind - j);
            hcp_squared_amplitudes[indexer] = amp1*pDen[j];
        }
        // across the row
        unsigned int indexer = ind*totalHCPInUse - (ind*(ind-1)/2);
        for(unsigned int j=ind; j<totalHCPInUse; j++){
            hcp_squared_amplitudes[indexer] = amp1*pDen[j];
            indexer++;
        }
    }

}


void DensityMapper::populateSquaredAmplitudes(unsigned int totalHCPInUse, std::vector<float> & densities_at_HCP, float * squared_amplitudes){
    float * const pDen = densities_at_HCP.data();
    //auto total_squared = totalHCPInUse*(totalHCPInUse-1)/2 + totalHCPInUse;
    /*
     * pre calculate squared amplitudes
     */
    float amp1;
    unsigned int indexer =0;
    // total => (totalHCPInUse*(totalHCPInUse-1)/2 + totalHCPInUse)
    for(unsigned int i=0; i<totalHCPInUse; i++){
        amp1 = pDen[i];
        for(unsigned int j=i; j<totalHCPInUse; j++){
            squared_amplitudes[indexer] = amp1*pDen[j];
            indexer++;
        }
    }
}

/*
 * calculate electron density value at each point in the keptHCPlattice points (sampling lattice)
 *
 */
void DensityMapper::populateDensitiesOMP(unsigned int totalHCPInUse, std::vector<Neighbors> neighbors, std::vector<float> & amplitudesT, std::vector<float> & densities_at_HCP){

    const float * const pAmp = amplitudesT.data(); // base amplitude for input MODEL
    float * const pDen = densities_at_HCP.data();

    Neighbors * const  pNeighborhood = neighbors.data();

    unsigned int totalNeighbors;
    float amplitude_r_i;
    float val;

    for(unsigned int i=0; i<totalHCPInUse; i++){
        // calculate electron density amplitude based on neighborhood
        auto & neighborhood_i = pNeighborhood[i];
        const float * model_kernel_distances = neighborhood_i.getModelDistances();
        const unsigned int * model_neighbors = neighborhood_i.getModelIndices();

        totalNeighbors = neighborhood_i.getTotalNeighbors();

        amplitude_r_i = 0.0;
        float inv_sum = 0.0f;
        for(int n=0; n<totalNeighbors; n++){ // kernel function
            val = model_kernel_distances[n];
            inv_sum += val;
            amplitude_r_i += pAmp[model_neighbors[n]] * val;
        } // pAmp is same length as input PDB model
        /*
         * can do as a sum or average of values of base lattice
         * doing inverse distance weighting interpolation
         */
        pDen[i] = amplitude_r_i/inv_sum;//max_neighbors/(float)totalNeighbors;
    }
}

/*
 * calculate electron density value at each point in the keptHCPlattice points (sampling lattice)
 *
 */
float DensityMapper::populateDensities(std::vector<float> & amplitudesT, std::vector<float> & densities_at_HCP){
    unsigned int totalHCPInUse = keptHCPSamplingPoints.size();

    const float * const pAmp = amplitudesT.data(); // base amplitude for input MODEL
    float * const pDen = densities_at_HCP.data();

    Neighbors * const  pNeighborhood = neighborhoods.data();

    unsigned int totalNeighbors;
    float amplitude_r_i;
    float sum = 0.0;
    float val;

    for(unsigned int i=0; i<totalHCPInUse; i++){
        // calculate electron density amplitude based on neighborhood
        auto & neighborhood_i = pNeighborhood[i];
        const float * model_kernel_distances = neighborhood_i.getModelDistances();
        const unsigned int * model_neighbors = neighborhood_i.getModelIndices();

        totalNeighbors = neighborhood_i.getTotalNeighbors();

        amplitude_r_i = 0.0;
        float inv_sum = 0.0;
        for(int n=0; n<totalNeighbors; n++){ // kernel function
            val = model_kernel_distances[n];
            inv_sum += val;
            amplitude_r_i += pAmp[model_neighbors[n]] * val;
        } // pAmp is same length as input PDB model
        sum += amplitude_r_i;
        /*
         * can do as a sum or average of neighbors
         */
        pDen[i] = amplitude_r_i/inv_sum;
//        if (val != 0 && pDen[i] != val){
//            std::cout << "HCP INDEX " << i << " " << val << " " << pDen[i] << " " << neighborhood_i.getHCPLatticeIndex() << " " << neighborhood_i.getKeptIndex() << std::endl;
//        }
    }

    //auto total_squared = totalHCPInUse*(totalHCPInUse-1)/2 + totalHCPInUse;
    /*
     * pre calculate squared amplitudes
     */
//    float amp1;
//    unsigned int indexer =0;
//    // (totalKept*(totalKept-1)/2 + totalKept)
//    for(unsigned int i=0; i<totalHCPInUse; i++){
//        amp1 = pDen[i];
//        for(unsigned int j=i; j<totalHCPInUse; j++){
//            squared_amplitudes[indexer] = amp1*pDen[j];
//            indexer++;
//        }
//    }

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

/*
 * Adaptive Simulated Annealing
 */
void DensityMapper::refineModelASA(unsigned int total_annealing_steps,
                                std::vector<Datum> & workingSet,
                                std::vector<Datum> & workingSetSmoothed ){
    /*
     *
     */
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(0.0,1.0);

    const unsigned int total_data = workingSet.size();
    Datum * const pWorkingSet = workingSet.data();
    Datum * const pWorkingSetSmoothed = workingSetSmoothed.data();
    std::vector<float> i_calc(total_data);
    std::vector<float> residuals(total_data);
    std::vector<float> qvalues;
    std::vector<float> sigmas_squared;

    for(auto & datum : workingSet){
        qvalues.push_back(datum.getQ());
        sigmas_squared.push_back(1.0f/(datum.getSigma()*datum.getSigma()));
    }

    float * const pICalc = i_calc.data();
    float * const pSigmaSquared = sigmas_squared.data();
    const LatticePoint * pLattice = lattice_points.data();
    auto total_lattice_points = lattice_points.size();
    amplitudes.resize(total_lattice_points); // holds the amplitudes in use

    unsigned int total_amplitudes_per_lattice_point = lattice_points[0].getTotalAmplitudes();

    std::vector<unsigned int> lattice_indices(total_lattice_points);
    std::uniform_int_distribution<unsigned int> randomAmpIndex(0,total_amplitudes_per_lattice_point-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomLatticePt(0,total_lattice_points-1); // guaranteed unbiased

    std::vector<unsigned int> amplitude_indices(total_amplitudes_per_lattice_point);
    for(unsigned int i=0; i<total_amplitudes_per_lattice_point; i++){
        amplitude_indices[i] = i;
    }

    unsigned int amp_index, amp_index2, amp_index3, selected_lattice_pt, selected_lattice_pt2, selected_lattice_pt3;
    std::vector<float> test_model_amplitudes(total_lattice_points);

    // set up initial random density model
    float * const pTestModelAmplitudes = test_model_amplitudes.data();
    float * const pAmp = amplitudes.data();

    for(int i=0; i < total_lattice_points; i++){
        pTestModelAmplitudes[i] = pLattice[i].getAmplitudeByIndex(randomAmpIndex(gen));
        pAmp[i] = pTestModelAmplitudes[i];
        lattice_indices[i] = i;
    }

    this->createHCPGrid(); //
    this->setMaxNeighbors();
    unsigned int long totalKept = keptHCPSamplingPoints.size();
    std::vector<float> hcp_electron_densities(totalKept);
    std::vector<float> prior_hcp_electron_densities(totalKept);

    unsigned int long size_of_dbf = (totalKept*(totalKept-1)/2 + totalKept)*workingSet.size();
    debye_factors = (float *)_aligned_malloc(sizeof(float)*size_of_dbf + 16, 16); // 16 byte aligned

    // use aligned memory for squared amplitudes and debye_factors
    auto * hcp_squared_amplitudes = (float *)_aligned_malloc(sizeof(float)*((totalKept*(totalKept-1)/2 + totalKept)), 16);

    this->setDebyeFactors(workingSet);

    //takes less than 100 usec
    float ave = populateDensities(test_model_amplitudes, hcp_electron_densities);
    populateSquaredAmplitudes(totalKept, hcp_electron_densities, hcp_squared_amplitudes);
    /*
    * by far the slowest part -
    */
    populateICalc(total_data, i_calc, hcp_squared_amplitudes);
    // need prior values from hcp_electron_densities in order to updateIcalc
    /*
     * takes less than 100 usec for scale and chi2
     */
    float scale = calculateScaleFactor(total_data, pICalc, pSigmaSquared, pWorkingSetSmoothed);
    float chi2, bestchi2 = getChiSquare(total_data, scale, pICalc, pSigmaSquared, pWorkingSet, pWorkingSetSmoothed, residuals.data());
    float dw, bestdw = calculateDurbinWatson(total_data, residuals.data());

    float dw_weight = 5.0f;
    float test_energy, current_energy = bestchi2+dw_weight*bestdw;

    auto temperature =  (double)0.1;
    double inv_kb_temp = 1.0/temperature;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;
    int counter = 1;
    std::vector<float> acceptanceRateDuringRun((unsigned long int)total_annealing_steps);
    std::vector<double> tempDuringRun((unsigned long int)total_annealing_steps);
    std::vector<float> scores((unsigned long int)total_annealing_steps);
    std::vector<float> chis((unsigned long int)total_annealing_steps);
    std::vector<int> flips((unsigned long int)total_annealing_steps);
    /*
     * in the first 25% of the annealing, randomly change 33%
     */
    unsigned int stop_limit, rampLimit = std::round(0.85*total_annealing_steps);
    unsigned int constantLimit = std::round(0.65*total_annealing_steps);

    auto fraction = (unsigned int)std::round(0.47*total_lattice_points);
    float slope = (1.0f - (float)fraction)/((float)(rampLimit - constantLimit));
    float intercept = 1.0f - slope*rampLimit;

    std::string status = "fixed";
//    std::vector<unsigned int> indicesInUse(1);
//    std::vector<unsigned int> selectedHCPSamplingPoints;//keptHCPSamplingPoints
//    std::vector<Density> priorHCPdensities;
    // step - temp - acc - score

    logger("TOTAL ANNEALING STEPS", formatNumber(total_annealing_steps));
    LatticePoint * pSelectedLattice;

    // parallel programming of this would divide total_annealing_steps into different processes and using MPI to update best solution
    int threads = 6;//omp_get_max_threads();
    // find whole number multiple of constantLimit divided by number threads
    unsigned int chunk = constantLimit / threads;
    double result = 0;

    unsigned int total_step = 0;
    double runtime;
    std::clock_t startTime = std::clock();
    double wall0 = omp_get_wtime();

    #pragma omp parallel for num_threads(threads)
    for(int t = 0; t < threads; t++){
        std::random_device omprd;
        std::mt19937 ompgen(omprd());
        std::uniform_real_distribution<> omp_distribution(0.0,1.0);

        unsigned int step = 0, omp_total_annealing_steps = total_annealing_steps;
        double sum = 0;

        float omp_scale, omp_chi2, omp_dw, omp_test_energy;
        float omp_dw_weight = dw_weight;
        float omp_current_energy = current_energy;
        auto inv_kb_temp_omp = inv_kb_temp;
        float omp_bestchi2, omp_bestdw;

        auto omp_temperature =  temperature;
        float omp_acceptRate = acceptRate, omp_inv500 = inv500, omp_inv500slash499 = inv500slash499;

        int sa_stop_limit = fraction;
        unsigned int totalHCPInUse = totalKept;
        const unsigned int omp_total_data = total_data;

        std::vector<Datum> omp_workingSet;
        std::vector<Datum> omp_workingSetSmoothed;

        std::vector<float> omp_i_calc(omp_total_data);
        std::vector<float> omp_residuals(omp_total_data);
        std::vector<float> omp_sigmas_squared(omp_total_data);
        for(auto & data : workingSet){
            omp_workingSet.emplace_back(Datum(data));
        }

        for(auto & data : workingSetSmoothed){
            omp_workingSetSmoothed.emplace_back(Datum(data));
        }

//        std::copy(workingSet.begin(), workingSet.end(), omp_workingSet.begin());
//        std::copy(workingSetSmoothed.begin(), workingSetSmoothed.end(), omp_workingSetSmoothed.begin());
//        std::copy(i_calc.begin(), i_calc.end(), omp_i_calc.begin());
//        std::copy(residuals.begin(), residuals.end(), omp_residuals.begin());
        std::copy(sigmas_squared.begin(), sigmas_squared.end(), omp_sigmas_squared.begin());
        // make copies of amplitudes, lattice indices, etc.
        std::vector<float> omp_amplitudes(total_lattice_points);
        std::vector<LatticePoint> omp_lattice_points;
        for(auto & lp : lattice_points){
            omp_lattice_points.emplace_back(LatticePoint(lp));
        }

        std::vector<unsigned int> omp_lattice_indices(total_lattice_points);
        std::vector<float> omp_hcp_electron_densities(totalKept);
        std::vector<float> omp_test_model_amplitudes(total_lattice_points);
        std::vector<Neighbors> omp_neighborhoods(neighborhoods);

        std::copy(amplitudes.begin(), amplitudes.end(), omp_amplitudes.begin());
        std::copy(amplitudes.begin(), amplitudes.end(), omp_test_model_amplitudes.begin());
//        std::copy(lattice_points.begin(), lattice_points.end(), omp_lattice_points.end());
        std::copy(lattice_indices.begin(), lattice_indices.end(), omp_lattice_indices.begin());
//        std::copy(neighborhoods.begin(), neighborhoods.end(), omp_neighborhoods.end());

        auto * omp_hcp_squared_amplitudes = (float *)_aligned_malloc(sizeof(float)*((totalHCPInUse*(totalHCPInUse-1)/2 + totalHCPInUse)), 16);

        Datum * const pOMPWorkingSet = omp_workingSet.data();
        Datum * const pOMPWorkingSetSmoothed = omp_workingSetSmoothed.data();

        float * omp_debye_factors = (float *)_aligned_malloc(sizeof(float)*size_of_dbf + 16, 16); // 16 byte aligned
        for(unsigned int m=0; m<size_of_dbf; m++){
            omp_debye_factors[m] = debye_factors[m];
        }


        while(step < chunk){ // perform constant temperature search broken into separate threads
            std::shuffle(omp_lattice_indices.begin(), omp_lattice_indices.end(), ompgen);
            unsigned int lat_index;
            for(int i=0; i < sa_stop_limit; i++){ // generate new state
                lat_index = omp_lattice_indices[i];
                omp_test_model_amplitudes[lat_index] = omp_lattice_points[lat_index].guessAmplitude(omp_distribution(ompgen));
            }

            populateDensitiesOMP(totalHCPInUse, omp_neighborhoods, omp_test_model_amplitudes, omp_hcp_electron_densities);
            populateSquaredAmplitudes(totalHCPInUse, omp_hcp_electron_densities, omp_hcp_squared_amplitudes);
            populateICalcOpenMP(omp_total_data, totalHCPInUse, omp_i_calc, omp_hcp_squared_amplitudes, omp_debye_factors);

            omp_scale = calculateScaleFactor(omp_total_data, omp_i_calc.data(), omp_sigmas_squared.data(), pOMPWorkingSetSmoothed);
            omp_chi2 = getChiSquare(omp_total_data, omp_scale, omp_i_calc.data(), omp_sigmas_squared.data(), pOMPWorkingSet, pOMPWorkingSetSmoothed, omp_residuals.data());
            omp_dw = calculateDurbinWatson(omp_total_data, omp_residuals.data());
            omp_test_energy = omp_chi2 + omp_dw_weight*omp_dw;

            if (omp_test_energy < omp_current_energy || (std::exp((omp_current_energy - omp_test_energy)*inv_kb_temp_omp) > omp_distribution(ompgen))) {
                omp_acceptRate = omp_inv500slash499*omp_acceptRate+omp_inv500;
                omp_current_energy = omp_test_energy;
                omp_bestchi2 = omp_chi2;
                omp_bestdw = omp_dw_weight*omp_dw;

                std::copy(omp_test_model_amplitudes.begin(), omp_test_model_amplitudes.end(), omp_amplitudes.begin());

            } else { // undo changes
                omp_acceptRate *= omp_inv500slash499;
                std::copy(omp_amplitudes.begin(), omp_amplitudes.end(), omp_test_model_amplitudes.begin());

            }

            updateASATemp(step, omp_total_annealing_steps, omp_acceptRate, omp_temperature, inv_kb_temp_omp);

            if ((step % 3999) == 0){
                // update amplitudes with best
                #pragma omp critical
                {
                    std::cout << "______________________________________________________________________________" << std::endl;
                    printf("T%-2i STEP %-6i %%CMPLTE %5.3f %-4i ENRGY %.3E TEMP %-.2E ACCPTRATE %5.3f  \n",
                           t,
                           step,
                           (float)total_step/(float)constantLimit,
                           sa_stop_limit,
                           omp_current_energy,
                           omp_temperature,
                           omp_acceptRate);
                    if (step > 0)
                        total_step += 3999;
                }
            }
            step++;
        }

        _aligned_free(omp_hcp_squared_amplitudes);
        _aligned_free(omp_debye_factors);
        // update best answer
        #pragma omp critical
        {
            if (omp_test_energy < current_energy){
                current_energy = omp_test_energy;
                bestchi2 = omp_bestchi2;
                bestdw = omp_bestdw;
                inv_kb_temp = inv_kb_temp_omp;
                acceptRate = omp_acceptRate;
                temperature = omp_temperature;
                std::copy(omp_amplitudes.begin(), omp_amplitudes.end(), amplitudes.begin());
            }
        }
    }

    // start cooling
    std::cout << " finished omp " << threads << std::endl;
    runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
    std::cout << "THREADS " << omp_get_num_threads() << "elapsed time= " << omp_get_wtime() - wall0 << " s :: "  << runtime << std::endl;
    std::copy(amplitudes.begin(), amplitudes.end(), test_model_amplitudes.begin());

    for (unsigned int step=constantLimit; step < total_annealing_steps; step++) {
        // grab a random lattice point
            std::shuffle(lattice_indices.begin(), lattice_indices.end(), gen);
//            stop_limit = (step < constantLimit) ? fraction : (unsigned int)std::round(slope*step + intercept);
            stop_limit = (step < rampLimit) ? (step < constantLimit) ? fraction : (unsigned int)std::round(slope*step + intercept) : 1;
            status = (step < constantLimit) ? "CONSTANT" : "RAMPING ";

            if (stop_limit > 7){
                unsigned int lat_index;
                for(int i=0; i < stop_limit; i++){ // generate new state
                    lat_index = lattice_indices[i];
                    pTestModelAmplitudes[lat_index] = pLattice[lat_index].guessAmplitude(distribution(gen));
                }
            } else { // when stop_limit gets low, changes of guessAmplitude remaining same per point is 1/4
                unsigned int lat_index;
                float test_amp, current_amp;
                for(int i=0; i < stop_limit; i++){ // generate new state
                    lat_index = lattice_indices[i];
                    current_amp = amplitudes[lat_index];
                    pSelectedLattice = &lattice_points[lat_index];
                    test_amp = pSelectedLattice->guessAmplitude(distribution(gen));
                    while(fabs(test_amp - current_amp) < 0.0001f ){
                        test_amp = pSelectedLattice->guessAmplitude(distribution(gen));
                    }
                    pTestModelAmplitudes[lat_index] = test_amp;
                }
            }

        // modify density about 100x faster than populateICalc
//        startTime = std::clock();
// about 100x faster than populateIcalc
        ave = populateDensities(test_model_amplitudes, hcp_electron_densities);
        populateSquaredAmplitudes(totalKept, hcp_electron_densities, hcp_squared_amplitudes);
//        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//        logger("populateDensities TIME", formatNumber((float)runtime,8));

//        startTime = std::clock();
//        updateICalc(priorHCPdensities, total_data, i_calc, prior_hcp_electron_densities.data(), hcp_electron_densities.data());
//        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//        logger("updateICalc TIME", formatNumber((float)runtime,8));

//        startTime = std::clock();
        populateICalc(total_data, i_calc, hcp_squared_amplitudes);
//        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//        logger("populateICalc TIME", formatNumber((float)runtime,8));

        scale = calculateScaleFactor(total_data, pICalc, pSigmaSquared, pWorkingSetSmoothed);
        chi2 = getChiSquare(total_data, scale, pICalc, pSigmaSquared, pWorkingSet, pWorkingSetSmoothed, residuals.data());
        dw = calculateDurbinWatson(total_data, residuals.data());
        test_energy = chi2 + dw_weight*dw;

        if (test_energy < current_energy || (std::exp((current_energy - test_energy)*inv_kb_temp) > distribution(gen))) {
            acceptRate = inv500slash499*acceptRate+inv500;
            current_energy = test_energy;
            bestchi2 = chi2;
            bestdw = dw_weight*dw;

            std::copy(test_model_amplitudes.begin(), test_model_amplitudes.end(), amplitudes.begin());
//            for(int i=0; i < stop_limit; i++){ // generate new state
//                lat_index = lattice_indices[i];
//                amplitudes[lat_index] = test_model_amplitudes[lat_index];
//            }
        } else { // undo changes
            acceptRate *= inv500slash499;
            std::copy(amplitudes.begin(), amplitudes.end(), test_model_amplitudes.begin());
//            for(int i=0; i < stop_limit; i++){ // generate new state
//                lat_index = lattice_indices[i];
//                test_model_amplitudes[lat_index] = amplitudes[lat_index];
//            }
        }


        // if accepted update index count in lattice point
        if ((step % 7999) == 0){
            // update amplitudes with best
            std::cout << "______________________________________________________________________________" << std::endl;
            printf("   STEP : %-7i %% COMPLETE : %5.3f  %8s  %-5i   TEMP : %-.3E \n", step, (float)step/(float)total_annealing_steps, status.c_str(), stop_limit, temperature);
            printf(" ENERGY : %.5f CHI2 : %-.3E DW : %-.3E ACCEPT_RATE : %5.3f  \n", current_energy, bestchi2, bestdw, acceptRate);

            this->createXPLORMap("best_"+std::to_string(counter));
            writeLatticePoints("best_"+std::to_string(counter));
//            ave = populateDensities(amplitudes, hcp_electron_densities, hcp_squared_amplitudes);
//            populateICalc(total_data, i_calc, hcp_squared_amplitudes);
//            writeICalc(total_data, scale, pICalc, pWorkingSet, pWorkingSetSmoothed, "best_" + std::to_string(counter));
            counter++;
        }

//        acceptanceRateDuringRun[step-constantLimit] = acceptRate;
//        scores[step-constantLimit] = current_energy;
//        chis[step-constantLimit] = bestchi2;
//        tempDuringRun[step-constantLimit] = bestdw;
//        flips[step-constantLimit] = stop_limit;
        updateASATemp(step, total_annealing_steps, acceptRate, temperature, inv_kb_temp);
    }

    printParameters(acceptanceRateDuringRun, tempDuringRun, scores, chis, flips);

    for (auto & pLP : lattice_points){
        pLP.calculateProbabilitiesFromOccurrences();
    }

    this->createXPLORMap("final");
    writeLatticePoints("final");

    ave = populateDensities(amplitudes, hcp_electron_densities);
    populateSquaredAmplitudes(totalKept, hcp_electron_densities, hcp_squared_amplitudes);
    populateICalc(total_data, i_calc, hcp_squared_amplitudes);

    scale = calculateScaleFactor(total_data, pICalc, pSigmaSquared, pWorkingSetSmoothed);

    writeICalc(total_data, scale, pICalc, pWorkingSet, pWorkingSetSmoothed, "final");

}


void DensityMapper::refineModelOPENMP(int max_threads, int max_rounds, float topPercent, int models_per_round,
                                std::vector<Datum> & workingSet,
                                std::vector<Datum> & workingSetSmoothed ){

    // parallel programming of this would divide total_annealing_steps into different processes and using MPI to update best solution
    models_per_round *= 7;
    int threads = max_threads;
    const unsigned int total_data = workingSet.size();
    Datum * const pWorkingSet = workingSet.data();
    Datum * const pWorkingSetSmoothed = workingSetSmoothed.data();

    float dw, updateAlpha = 0.67;

    auto topN = (unsigned int)(std::ceil((float)models_per_round*topPercent));
    std::vector<Trial> topTrials(topN);
    for(auto & trial : topTrials){
        trial.value = DBL_MAX;
    }

    // find whole number multiple of constantLimit divided by number threads
    unsigned int round = 0, chunk = models_per_round / threads;
    this->createHCPGrid(); //
    this->setMaxNeighbors();
    unsigned int long totalKept = keptHCPSamplingPoints.size();
    std::vector<float> hcp_electron_densities(totalKept);

    unsigned int long size_of_dbf = (totalKept*(totalKept-1)/2 + totalKept)*workingSet.size();
    debye_factors = (float *)_aligned_malloc(sizeof(float)*size_of_dbf + 16, 16); // 16 byte aligned

    // use aligned memory for squared amplitudes and debye_factors
    auto * squared_amplitudes = (float *)_aligned_malloc(sizeof(float)*((totalKept*(totalKept-1)/2 + totalKept)), 16);
    this->setDebyeFactors(workingSet);

    auto total_amplitudes_from_lattice_points = lattice_points.size();
    std::vector<float> prior_best_amplitudes(total_amplitudes_from_lattice_points, 0.25f);

    std::vector<float> i_calc(total_data);
    std::vector<float> sigmas_squared;

    for(auto & datum : workingSet){
        sigmas_squared.push_back(1.0f/(datum.getSigma()*datum.getSigma()));
    }
    float * const pSigmaSquared = sigmas_squared.data();

    int total_sets_in_residuals=1;
    int total_residuals = total_data*total_sets_in_residuals; // only look at top 5 models
    std::vector<float> residuals(total_residuals);
    float * const pRes = residuals.data();

    std::vector<float>durbinWatsons(max_rounds);
    std::vector<double>best(max_rounds); // average variation?
    std::clock_t startTime;
    double runtime;

    std::cout << "PER THREAD " << chunk << std::endl;
    std::cout << " PER ROUND " << models_per_round << std::endl;

    for(; round<max_rounds; round++){

        bool added_already = false;

        #pragma omp parallel for num_threads(threads)
        for(int th=0; th<threads; th++){
            // make copies of everything
            auto omp_topN = topN;
            std::vector<Trial> omp_topTrials(topN);
            const unsigned int last = topN-1;

            float scale, chi2, omp_dw;
            std::random_device omprd;
            std::mt19937 ompgen(omprd());
            std::uniform_real_distribution<> omp_distribution(0.0,1.0);
            unsigned int totalHCPInUse = totalKept;
            const unsigned int omp_total_data = total_data;

            std::vector<float> omp_i_calc(omp_total_data);
            std::vector<float> omp_residuals(omp_total_data);
            std::vector<float> omp_sigmas_squared(sigmas_squared);

            std::vector<Datum> omp_workingSet;
            std::vector<Datum> omp_workingSetSmoothed;
            for(auto & data : workingSet){
                omp_workingSet.emplace_back(Datum(data));
            }

            for(auto & data : workingSetSmoothed){
                omp_workingSetSmoothed.emplace_back(Datum(data));
            }

            std::vector<Neighbors> omp_neighborhoods(neighborhoods);

            auto omp_total_amplitudes_from_lattice_points = total_amplitudes_from_lattice_points;
            std::vector<float> omp_test_model_amplitudes(omp_total_amplitudes_from_lattice_points);
            std::vector<float> omp_hcp_electron_densities(totalHCPInUse);
            auto * omp_hcp_squared_amplitudes = (float *)_aligned_malloc(sizeof(float)*((totalHCPInUse*(totalHCPInUse-1)/2 + totalHCPInUse)), 16);

            std::vector<LatticePoint> omp_lattice_points(lattice_points);
//            for(auto & lp : lattice_points){
//                omp_lattice_points.emplace_back(LatticePoint(lp));
//            }

            const LatticePoint * pLattice = omp_lattice_points.data();
            Datum * const pOMPWorkingSet = omp_workingSet.data();
            Datum * const pOMPWorkingSetSmoothed = omp_workingSetSmoothed.data();

            float * omp_debye_factors = (float *)_aligned_malloc(sizeof(float)*size_of_dbf + 16, 16); // 16 byte aligned
            for(unsigned int m=0; m<size_of_dbf; m++){
                omp_debye_factors[m] = debye_factors[m];
            }

            unsigned int omp_topAdded=0;

            for(int model=1; model<chunk; model++){
                float * pPriors = omp_test_model_amplitudes.data();
                for(int i=0; i < omp_total_amplitudes_from_lattice_points; i++){
                    pPriors[i] = pLattice[i].guessAmplitude(omp_distribution(ompgen));
                }

                // based on guess, anything with Euler tour > 1 should be ingored
                std::set<unsigned int> selectedIndices;

                for(int i=0; i < omp_total_amplitudes_from_lattice_points; i++){
                    if (pPriors[i] > 0){
                        selectedIndices.insert(i);
                    }
                }

                populateDensitiesOMP(totalHCPInUse, omp_neighborhoods, omp_test_model_amplitudes, omp_hcp_electron_densities);
                populateSquaredAmplitudes(totalHCPInUse, omp_hcp_electron_densities, omp_hcp_squared_amplitudes);
                populateICalcOpenMP(omp_total_data, totalHCPInUse, omp_i_calc, omp_hcp_squared_amplitudes, omp_debye_factors);
//
                scale = calculateScaleFactor(omp_total_data, omp_i_calc.data(), omp_sigmas_squared.data(), pOMPWorkingSetSmoothed);
                chi2 = getChiSquare(omp_total_data, scale, omp_i_calc.data(), omp_sigmas_squared.data(), pOMPWorkingSet, pOMPWorkingSetSmoothed, omp_residuals.data());
                omp_dw = calculateDurbinWatson(omp_total_data, omp_residuals.data());
                chi2 += 2.0f*omp_dw;

                if (!selectedIndices.empty() && !checkConnectivity(selectedIndices, omp_lattice_points)) {
                    chi2 *= 1.5f;
                }

                // update topN
                if (omp_topAdded < omp_topN){
                    Trial * pTrial = &omp_topTrials[omp_topAdded];
                    pTrial->value = chi2;
                    pTrial->dw = omp_dw;
                    pTrial->counter = 1;
                    pTrial->scale = scale;
                    pTrial->model_amplitudes.swap(omp_test_model_amplitudes);
                    pTrial->residuals.swap(omp_residuals);

                    omp_test_model_amplitudes.resize(omp_total_amplitudes_from_lattice_points);
                    omp_residuals.resize(total_data);

                    omp_topAdded++;
                    if (omp_topAdded == omp_topN){
                        std::sort(omp_topTrials.begin(), omp_topTrials.begin()+omp_topAdded);
                    }
                } else {
                    if (chi2 < omp_topTrials[last].value){
                        /*
                         * replace last entry and sort
                         */
                        Trial * pTrial = &omp_topTrials[last];
                        pTrial->value = chi2;
                        pTrial->dw = omp_dw;
                        pTrial->counter = 1;
                        pTrial->scale = scale;
                        pTrial->model_amplitudes.swap(omp_test_model_amplitudes);
                        pTrial->residuals.swap(omp_residuals);

                        omp_test_model_amplitudes.resize(omp_total_amplitudes_from_lattice_points);
                        omp_residuals.resize(total_data);
                        std::sort(omp_topTrials.begin(), omp_topTrials.end());
                    }
                }
            }

            _aligned_free(omp_hcp_squared_amplitudes);
            _aligned_free(omp_debye_factors);
            // at the end, update the master topN list
            #pragma omp critical
            {
                if (added_already){
                    for(auto & testTrial : omp_topTrials){

                        Trial * pTrialLast = &topTrials[last];

                        if (testTrial.value < pTrialLast->value){
                            /*
                             * replace last entry and sort
                             */
                            pTrialLast->value = testTrial.value;
                            pTrialLast->dw = testTrial.dw;
                            pTrialLast->counter = 1;
                            pTrialLast->scale = testTrial.scale;
                            pTrialLast->model_amplitudes.swap(testTrial.model_amplitudes);
                            pTrialLast->residuals.swap(testTrial.residuals);

                            std::sort(topTrials.begin(), topTrials.end());
                        }
                    }
                } else {
                    std::cout << " INITIALIZING TopN THREAD :: " << th << std::endl;
                    std::copy(omp_topTrials.begin(), omp_topTrials.end(), topTrials.begin()); // already sorted
                    added_already = true;
                }
            }
        } // end of threads
        added_already = false;
        /*
         * update probability model
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
            point.updateProbabilities(updateAlpha);
        }

        // calculate RMSD of amplitudes
//        float rmsd = calculate_rmsd(prior_best_amplitudes);
//        // reassign prior to current
//        for(unsigned int i=0; i < total_amplitudes_from_lattice_points; i++){
//            prior_best_amplitudes[i] = lattice_points[i].getWeightedAmplitude();
//        }

        auto counter = (float)probability_count(0.67)/(float)lattice_points.size();

        // calculate DurbinWatson for the top5
        auto & trial = topTrials[0];
//        unsigned int indexer = 0;
//        for(auto & res : trial.residuals){
//            pRes[indexer] = res;
//            indexer++;
//        }

        dw = trial.dw;//calculateDurbinWatson(total_residuals, residuals.data());
        durbinWatsons[round] = dw;
        best[round] = topTrials[0].value;

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
        startTime=std::clock(); // reset clock

        logger("ROUND", formatNumber(round));
        logger("CHI2 TOP BEST", formatNumber((float)topTrials[0].value, 4));
        logger("CHI2 TOP LAST", formatNumber((float)topTrials[topN-1].value, 4));
        logger("% > than 50", formatNumber(counter, 4));

        // update amplitudes with best
        float * pDens = topTrials[0].model_amplitudes.data();
        for(int i=0; i<lattice_points.size(); i++){
            amplitudes[i] = pDens[i];
        }

        this->createXPLORMap("best_"+std::to_string(round));
        writeLatticePoints("best_"+std::to_string(round));

        float ave = populateDensities(amplitudes, hcp_electron_densities);
        populateSquaredAmplitudes(totalKept, hcp_electron_densities, squared_amplitudes);
        populateICalc(total_data, i_calc, squared_amplitudes);

        float scale = topTrials[0].scale;
        writeICalc(total_data, scale, i_calc.data(), pWorkingSet, pWorkingSetSmoothed, "best_" + std::to_string(round));

        logger("DurbinWatson", formatNumber(dw, 5));
        logger("Time Per Round", formatNumber((float)runtime/60.0f,1));

        // write top 3
        for (int m=1; m<4; m++){
            pDens = topTrials[m].model_amplitudes.data();
            for(int i=0; i<lattice_points.size(); i++){
                amplitudes[i] = pDens[i];
            }
            this->createXPLORMap("top_"+std::to_string(m));
        }
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
    std::vector<float> prior_best_amplitudes(total_amplitudes_from_lattice_points, 0.25f);
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

    int total_sets_in_residuals=1;
    int total_residuals = total_data*total_sets_in_residuals; // only look at top 5 models
    std::vector<float> residuals(total_residuals);
    float * const pRes = residuals.data();

    float * const pICalc = i_calc.data();
    float * const pSigma = sigmas_squared.data();
    const LatticePoint * pLattice = lattice_points.data();
    /*
     * choose random lattice_point set
     */
//    float cv_fraction = 0.2;
//    int total_random_lattice_points = cv_fraction*lattice_points.size();
//    std::vector<unsigned int> lattice_indices = std::vector<unsigned int>(total_random_lattice_points);

    float chi2, scale, ave, dw;

    //this->setBessels(qvalues);
    this->createHCPGrid(); //
    this->setMaxNeighbors();
    unsigned int long totalKept = keptHCPSamplingPoints.size();
    std::vector<float> hcp_electron_densities(totalKept);

    unsigned int long size_of_dbf = (totalKept*(totalKept-1)/2 + totalKept)*workingSet.size();
//    debye_factors = new alignas(16) float[size_of_dbf];
    debye_factors = (float *)_aligned_malloc(sizeof(float)*size_of_dbf + 16, 16); // 16 byte aligned

    // use aligned memory for squared amplitudes and debye_factors
    //std::vector<float> squared_amplitudes(keptHCPSamplingPoints.size()*(keptHCPSamplingPoints.size()-1)/2 + keptHCPSamplingPoints.size());
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
    models_per_round *= 7;
    logger("Initial Search Trials per Round", formatNumber((unsigned)models_per_round));
    std::vector<float>durbinWatsons(max_rounds);
    std::vector<float>variations(max_rounds);
    std::vector<double>best(max_rounds); // average variation?

    int stretch = 4;
    std::vector<Trial> topTrialsForAveraging(stretch);
    double variation = 0;
    unsigned int round=0;

    std::ofstream logfile;
    logfile.open("logs_cemap.txt", std::ios_base::app);
    logfile << "# round best DW variations\n";
    logfile.close();
    char buffer[80];

//    float alpha = 100;
    float cutoff_limit, temp_amp, updateAlpha = 0.67;//, inv_total_lattice = 1.0f/(float)total_amplitudes_from_lattice_points;
    cutoff_limit = lattice_points[0].getLastAmplitude();
    std::set<unsigned int> selected_indices;
    unsigned int neighborLimit = inputBaseModel.getNeighborLimit();

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

            //takes less than 100 usec
//            startTime = std::clock();
            ave = populateDensities(prior_model_amplitudes, hcp_electron_densities);
            populateSquaredAmplitudes(totalKept, hcp_electron_densities, squared_amplitudes);

            // calculate intensities
//            startTime = std::clock();
            populateICalc(total_data, i_calc, squared_amplitudes);
//            runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//            logger("ICALC TIME", formatNumber((float)runtime,8));

            //takes less than 100 usec
            scale = calculateScaleFactor(total_data, pICalc, pSigma, pWorkingSetSmoothed);
            chi2 = getChiSquare(total_data, scale, pICalc, pSigma, pWorkingSet, pWorkingSetSmoothed, trial_residuals.data());

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
            point.updateProbabilities(updateAlpha);
        }

        // calculate RMSD of amplitudes
        float rmsd = calculate_rmsd(prior_best_amplitudes);
        // reassign prior to current
        for(unsigned int i=0; i < total_amplitudes_from_lattice_points; i++){
//            prior_best_amplitudes[i] = pLattice[i].getAssignedAmplitude();
            prior_best_amplitudes[i] = pLattice[i].getWeightedAmplitude();
        }

        auto counter = probability_count(0.7);

        // calculate DurbinWatson for the top5
        int indexer=0;
        for(int top=0; top<total_sets_in_residuals; top++){
            // copy residuals into total_residuals
            auto & trial = topTrials[top];
            for(auto & res : trial.residuals){
                pRes[indexer] = res;
                indexer++;
            }
        }

        dw = calculateDurbinWatson(total_residuals, residuals.data());
        durbinWatsons[round] = dw;
        best[round] = topTrials[0].value;

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

//        printf("             ROUND : %-7i %% COMPLETE : %5.3f MODELS_PER_ROUND %i \n", round, (float)round/(float)max_rounds, models_per_round);
//        printf("              TIME : %.1f (min) threshold count %i \n", (float)runtime/60.0f, counter);
//        printf("       CHI2 (best) : %-.4f (last) : %-.4f  \n", (float)topTrials[0].value, (float)topTrials[last].value);

        logger("ROUND", formatNumber(round));
        logger("CHI2 TOP BEST", formatNumber((float)topTrials[0].value, 4));
        logger("RMSD", formatNumber(rmsd, 6));
        logger("CHI2 TOP LAST", formatNumber((float)topTrials[last].value, 4));
        logger("COUNTER", formatNumber(counter));
        //lattice_points[lattice_points.size()/2].printProbabilities();

        float timeremaining = ((float)runtime*(float)(max_rounds - round)/60.0f);

        // update amplitudes with best
        float * pDens = topTrials[0].model_amplitudes.data();
        // update amplitudes with best
        for(int i=0; i<lattice_points.size(); i++){
            pAmp[i] = pDens[i];
        }

        this->createXPLORMap("best_"+std::to_string(round));
        writeLatticePoints("best_"+std::to_string(round));

        ave = populateDensities(amplitudes, hcp_electron_densities);
        populateSquaredAmplitudes(totalKept, hcp_electron_densities, squared_amplitudes);
        populateICalc(total_data, i_calc, squared_amplitudes);

        scale = topTrials[0].scale;
        writeICalc(total_data, scale, pICalc, pWorkingSet, pWorkingSetSmoothed, "best_" + std::to_string(round));

        logger("DurbinWatson", formatNumber(dw, 5));
        logger("Time Per Round", formatNumber((float)runtime/60.0f,1));
        logger("EST TIME REMAINING", formatNumber(timeremaining, 2) + " minutes");

        // track the last 4 best for averaging
        topTrialsForAveraging[round%stretch].model_amplitudes.swap(topTrials[0].model_amplitudes);
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
                    if (tempval > variation){ // find max variation
                        variation = tempval;
                    }
                }
            }

            logger("Last 4 MAX DIFF", formatNumber(variation, 3));
            variations[round] = variation;
            if ( (rmsd < 0.1 && dw < 0.1 && topTrials[0].value < 0.001) || round == (max_rounds-1)){
                pDens = topTrialsForAveraging[0].model_amplitudes.data();
                for(int i=0; i<lattice_points.size(); i++){
                    pAmp[i] = pDens[i];
                }

               // make the average of best from last 4 rounds
                for(int t=1; t<stretch; t++){
                   // update amplitudes with best
                    pDens = topTrialsForAveraging[t].model_amplitudes.data();
                    for(int i=0; i<lattice_points.size(); i++){
                        pAmp[i] += pDens[i];
                    }
                }

                for(int i=0; i<lattice_points.size(); i++){
                    pAmp[i] /= (float)stretch;
                }

                this->createXPLORMap("best_averaged");
                writeLatticePoints("best_averaged");
                break;
            }
            variation = 0;
        }

        std::snprintf(buffer, 80, "%8i %.3f %.4f %.5f %.5f\n", round, best[round], durbinWatsons[round], rmsd, variation);
        logfile.open("logs_cemap.txt", std::ios_base::app);
        logfile << buffer;
        logfile.close();

        std::cout << " " << std::endl;

        /*
         * reset models per round to base value after 3 rounds
         */
        if (round == 4 && models_per_round >= base_models_per_round){
            models_per_round = base_models_per_round;
            updateAlpha = 0.67;
        }

    }

    _aligned_free(squared_amplitudes);
    _aligned_free(debye_factors);
}


void DensityMapper::printLatticePointsInfo(){

    for(auto & point : lattice_points){
        //point.printProbabilities();
        point.PDBLineForOutput(centered_coordinates[point.getIndex()]);
    }
}



std::string DensityMapper::headerParametersForXPLOR(int & na, int & nb, int & nc, float grid_spacing){

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

std::string DensityMapper::headerParametersForXPLORFlipped(int & na, int & nb, int & nc, float grid_spacing){

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

    float grid_spacing = delta_r; // should match HCP grid

    //float inv_bandwidth = 1.0f/(cutoff/2.355f); // should be same as fwhm_sigma?
//    float inv_bandwidth = 1.0f/fwhm_sigma;
//    float inv_bandwidth_squared = 0.5f*inv_bandwidth*inv_bandwidth;
//    float inv_bandwidth_2PI = 1.0f/sqrtf(2.0f*M_PI)*inv_bandwidth;

    float mapSum = 0.0f, mapSumSquared = 0.0f;
    float mapCount = 0.0f;

    const vector3 * pVec = centered_coordinates.data();
    const float * pAmp = amplitudes.data();

    std::string tempHeader = headerParametersForXPLOR(na, nb, nc, grid_spacing);

    char buffer[80];

    std::vector<float> averages;
    std::vector<vector3> coords_averages;
    std::map<std::string, float > kernel_mapping;
    float val;

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
                float divisor = 0.0f;

                for(unsigned int c=0; c<total_centered_coordinates; c++){ // calculate distances to vectors in compiled model
                    /*
                     * centered_coordinates -> model
                     * lattice_points -> amplitudes assigned to each point in model
                     */
                    vector3 vec1 = *(pVec+c) - locale;
                    float length = vec1.length();
                    /*
                     * Consider point that coincides with base lattice at density of 0.3
                     * if we average over all lattice points, then the point will be much less than 0.3
                     *
                     */

                    if (length < upper_limit){
                        val = convolutionFunction(length);
                        kernelSum += *(pAmp + c)*val;
                        divisor += val;
                    }
                }

                divisor = divisor < 1 ? 1.0f : divisor;
                kernelSum /= divisor;
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
    tempHeader = "REMARK 256  4 sigma filter\n";

    float average = (mapSumSquared/mapCount);
    float stdev = std::sqrt(mapSumSquared/mapCount - ave*ave);
    float aboveit = average + 2.3f*stdev;

    int index = 1;
    char segids[] = { 'A', 'B', 'C', 'D'};
    int seg = 0;

    for(int i=0; i<averages.size(); i++){
        float value = averages[i];
        if(value >= aboveit){ // print
            vector3 & vec = coords_averages[i];
            std::snprintf(buffer, 80, "%-6s%5i %4s %3s %1c%4i    %8.3f%8.3f%8.3f  1.00%5.2f\n", "ATOM", index," CA ", "ALA", segids[seg], index, vec.x, vec.y, vec.z, value );
            tempHeader.append(buffer);
            index++;
            if (index > 9999){
                index = 1;
                seg += 1;
            }
        }
    }

    tempHeader.append("END\n");
    std::string out = name+"_filtered.pdb";
    outputFileName = out.c_str();

    pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);

    /*
     * calculate flipped map
     * flip the map x -> -x
     * reverse the order for each section
     */
    tempHeader = headerParametersForXPLORFlipped(na, nb, nc, grid_spacing);
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
    tempHeader = "REMARK 256  4 sigma filter\n";
    tempHeader.append("REMARK 256  Flipped (inverted) by changing X to -X for each point\n");

    index = 1;
    seg = 1;
    for(int i=0; i<averages.size(); i++){
        float value = averages[i];
        if(value >= aboveit){ // print
            vector3 & vec = coords_averages[i];
            std::snprintf(buffer, 80, "%-6s%5i %4s %3s %1c%4i    %8.3f%8.3f%8.3f  1.00%5.2f\n", "ATOM", 1," CA ", "ALA", segids[seg], index, -vec.x, vec.y, vec.z, value );
            tempHeader.append(buffer);
            index++;
            if (index > 9999){
                index = 1;
                seg += 1;
            }
        }
    }
    tempHeader.append("END\n");
    out = name+"_filtered_flipped.pdb";
    outputFileName = out.c_str();

    pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);
}


void DensityMapper::writeLatticePoints(std::string name){

    float sum = 0.0f;
    float sum2 = 0.0f;
    float tamp;
    for(auto & lat : lattice_points){
//        lat.setWeightedAmplitude();
//        tamp = lat.getWeightedAmplitude();
//        tamp = lat.getMostProbableAmplitude();
        tamp = lat.getMaxProbability();
        sum += tamp;
        sum2 += tamp*tamp;
    }

    float avg = sum/(float)lattice_points.size();
    float stdev = sqrt(sum2/(float)lattice_points.size() - avg*avg);

    char buffer[80];
    std::string tempHeader = "";
    std::snprintf(buffer, 80, "REMARK AVG %.3f\n", avg);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK Total Points in Base %5i\n", lattice_points.size());
    tempHeader.append(buffer);
    int index =1;
    /*
     * centered_coordinates holds coordinates lattice_point
     */
    for(auto & lat : lattice_points){
        //if (lat.getMaxProbability() < 0.5){ // outputs beads of highest certainty - regardless of amplitude
            vector3 & vec = centered_coordinates[lat.getIndex()];
            std::snprintf(buffer, 80, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00%5.2f\n", "ATOM", 1," CA ", "ALA", "A", index, vec.x, vec.y, vec.z, lat.getMaxProbability() );
            tempHeader.append(buffer);
            index++;
        //}
    }

    std::string out = name + "_prob.pdb";
    const char * outputFileName = out.c_str();

    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);
}

void DensityMapper::createHCPGrid(){
    // read in pdb file
    modelHCPSamplingDensityPoints = PointSetModel((2 * kmax), delta_r); // grid is more closely spaced than input model
    modelHCPSamplingDensityPoints.clearNeighborsAndBins();

    //const float beadradius_limit = (delta_r*std::sqrt(3.0f/2.0f));
    const float beadradius_limit = (delta_r+0.5f*cutoff);

    auto totalBeads = modelHCPSamplingDensityPoints.getTotalNumberOfBeadsInUniverse();

    // clear bins and neighbors in PointSetModel

    vector3 * tempvec;
    // keep lattice points that overlap with input bead model
    for(unsigned int i=0; i < totalBeads; i++){
        auto bead = modelHCPSamplingDensityPoints.getBead(i);
        auto bead_vec = bead->getVec();

        for (unsigned int a=0; a < total_centered_coordinates; a++){

            tempvec = &centered_coordinates[a];

            if ( (bead_vec - *tempvec).length() < beadradius_limit){ // if bead overlaps with atom, keep it
                keptHCPSamplingPoints.push_back(i);
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


    for(auto & index : keptHCPSamplingPoints){

        // make new neighborhood for each kept point
        neighborhoods.emplace_back(Neighbors(index));
        // get neighborhood
        Neighbors & lastOne = neighborhoods[neighborhoods.size() - 1];
        lastOne.setKeptIndex(neighborhoods.size() - 1);

        auto bead = modelHCPSamplingDensityPoints.getBead(index);
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
            auto iter = neighboringHCPPointsOfModelLattice.find((a));

            if (dis < upper_limit){
                float val = convolutionFunction(dis);
                lastOne.add_neighbor(a, val);
                iter->second.insert(index);
                ncoount+=val;
            }
        }
        std::snprintf(buffer, 80, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00%5.2f\n", "ATOM", 1," CA ", "ALA", colour[ind], indexer, bead_vec.x, bead_vec.y, bead_vec.z, ncoount );
        tempHeader.append(buffer);
        indexer++;

        if (indexer%9999==0){
            indexer=1;
            ind += 1;
        }
    }

    logger("Total KEPT HCP Lattice", formatNumber((unsigned int)keptHCPSamplingPoints.size()));

    std::string name = "hcp_grid_model.pdb";
    const char * outputFileName = name.c_str();
    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);
}

void DensityMapper::setMaxNeighbors(){
    max_neighbors = 0.0f;
    for(auto & ng : neighborhoods){
        if (ng.getTotalNeighbors() > max_neighbors){
            max_neighbors = ng.getTotalNeighbors();
        }
    }
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
 * keptHCPSamplingPoints must be populated before use
 *
 */
void DensityMapper::setDebyeFactors(std::vector<Datum> &workingSet) {

    Bead * bead2;
    float qr;

    int totalKept = keptHCPSamplingPoints.size();
    // total keptHCPSamplingPoints is same size as neighborhoods

    //debye_factors;
    std::vector<float> distances;

    // precompute all distances
    for (int i=0; i<totalKept; i++){
        auto bead = modelHCPSamplingDensityPoints.getBead(keptHCPSamplingPoints[i]);
        auto bead_vec = bead->getVec();
        distances.push_back(0.0f); // distances[ i*totalKept - (i*(i-1)/2)] => 0
        int next = i+1;
        for (int j=next; j<totalKept; j++){
            bead2 = modelHCPSamplingDensityPoints.getBead(keptHCPSamplingPoints[j]);
            distances.push_back((bead_vec - bead2->getVec()).length());
        }
    }

    // populate debye factors
    unsigned int long counter = 0;
    for(auto & data : workingSet){

        float q_value = data.getQ();
        unsigned int indexer = 0;
        for (int i=0; i<totalKept; i++){
            debye_factors[counter] = 1.0f;
            counter++;

            indexer++;
            for (int j=(i+1); j<totalKept; j++){
                qr = q_value*distances[indexer];
                indexer++;

                debye_factors[counter] = sinf(qr)/qr;
                counter++;
            }
        }
    } // very last entry of debye_factors should be 1
}

void DensityMapper::populateICalcOpenMP(unsigned int total_q, unsigned int totalHCPInUse, std::vector<float> & iCalc, float * squared_amplitudes, float * db) {

    float intensity_r_at_zero=0.0f;

    // diagonal entries will be hard to optimize as they are spaced out by large blocks in the aligned memory
//    for (int i = 0; i < totalHCPInUse; i++) { // calculate diagonal term of Debye matrix
//        intensity_r_at_zero += squared_amplitudes[ i*totalHCPInUse - (i*(i-1)/2)];// pDebye is 1 for self-amplitude squared
//    }
//    const unsigned int total_width = totalHCPInUse*(totalHCPInUse-1)/2+totalHCPInUse;
//    float * pICalc = iCalc.data();
//    for(unsigned int q_index = 0; q_index < total_q; q_index++) {
//        iCalc[q_index] = -intensity_r_at_zero;
//    }
//    #pragma omp parallel for num_threads(4)
//    for(unsigned int qw_index = 0; qw_index < total_q*total_width; qw_index++) {
//        iCalc[qw_index/total_width] += 2*squared_amplitudes[qw_index%total_width]*debye_factors[qw_index];
//    }

    //totalKept*(totalKept-1)/2 + totalKept
    for (int i = 0; i < totalHCPInUse; i++) { // calculate diagonal term of Debye matrix
        intensity_r_at_zero += squared_amplitudes[ i*totalHCPInUse - (i*(i-1)/2)];// pDebye is 1 for self-amplitude squared
    }

    const unsigned int total_width = totalHCPInUse*(totalHCPInUse-1)/2 + totalHCPInUse;

    float intensity2;
    unsigned int locale;

    float * pICalc = iCalc.data();
    unsigned int window = 8, next_i; // window is 8 for _m128 and 16 for _m256
    auto end_N = (unsigned int)std::floor(total_width/window);
    const unsigned int start_of_tail = end_N*window;

    simde__m128 mmSum1, mmSum2, v0, v1, v2, v3;//, s01, s23;

    auto * icalc_sse = (float *)_aligned_malloc(sizeof(float)*4, 16);
    float * pDF;

    for(unsigned int q_index = 0; q_index < total_q; q_index++) {

        next_i  = q_index * total_width;
        locale = next_i + start_of_tail;
        mmSum1 = simde_mm_setzero_ps();
        mmSum2 = simde_mm_setzero_ps();

        pDF = &db[next_i];

        for (unsigned int i=0; i < start_of_tail; i+=window ){ // SSE intrinsic here would multiply and then sum
            //_m128 is four floating point numbers (4x32 = 128 bits) => 16 bytes total
            // one byte is equivalent to eight bits
            v0 = simde_mm_load_ps(&squared_amplitudes[i]); // loading 4 floats at a time or 16 bytes total
            v1 = simde_mm_load_ps(pDF);
            mmSum1 = simde_mm_add_ps(mmSum1, simde_mm_mul_ps(v0,v1));

            v2 = simde_mm_load_ps(&squared_amplitudes[i + 4]); // loading next 4 floats for total 8 per cycle

            pDF += 4;
            v3 = simde_mm_load_ps(pDF);
            mmSum2 = simde_mm_add_ps(mmSum2, simde_mm_mul_ps(v2,v3));
            pDF += 4;
        }

        mmSum1 = simde_mm_add_ps(mmSum1, mmSum2);;
        simde_mm_store_ps(icalc_sse, mmSum1);

        intensity2 = icalc_sse[0] + icalc_sse[1] + icalc_sse[2] + icalc_sse[3];// + icalc_sse[4] + icalc_sse[5] + icalc_sse[6] + icalc_sse[7];

        // add remaining bits since we can't do in full width of 8
        pDF = &db[locale];
        for (unsigned int i=start_of_tail; i < total_width; i++ ){ // SSE intrinsic here would multiple and then sum
            //intensity2 += squared_amplitudes[i] * debye_factors[locale + (i - start_of_tail)];
            intensity2 += squared_amplitudes[i] * *(pDF);
            ++pDF;
        }

        // diagonal terms are included in intensity2, must remove their contribution
        *(pICalc) = -intensity_r_at_zero + 2.0f*intensity2;
        ++pICalc;
    }
}




void DensityMapper::populateICalc(unsigned int total_q, std::vector<float> & iCalc, float * squared_amplitudes) {

    unsigned int totalHCPInUse = keptHCPSamplingPoints.size(); // shared memory

    float intensity_r_at_zero=0.0f;
    //totalKept*(totalKept-1)/2 + totalKept
    for (int i = 0; i < totalHCPInUse; i++) { // calculate diagonal term of Debye matrix
        intensity_r_at_zero += squared_amplitudes[ i*totalHCPInUse - (i*(i-1)/2)];// pDebye is 1 for self-amplitude squared
    }

    const unsigned int total_width = totalHCPInUse*(totalHCPInUse-1)/2 + totalHCPInUse;

    float intensity2;
    unsigned int locale;

    float * pICalc = iCalc.data();
    unsigned int window = 8, next_i; // window is 8 for _m128 and 16 for _m256
    auto end_N = (unsigned int)std::floor(total_width/window);
    const unsigned int start_of_tail = end_N*window;

    simde__m128 mmSum1, mmSum2, v0, v1, v2, v3;//, s01, s23;
//    simde__m256 mmSum1, mmSum2; // vector of 8 floating point numbers

    auto * icalc_sse = (float *)_aligned_malloc(sizeof(float)*4, 16);
    float * pDF;

    for(unsigned int q_index = 0; q_index < total_q; q_index++) {

        next_i  = q_index * total_width;
        locale = next_i + start_of_tail;
        mmSum1 = simde_mm_setzero_ps();
        mmSum2 = simde_mm_setzero_ps();

        pDF = &debye_factors[next_i];

        for (unsigned int i=0; i < start_of_tail; i+=window ){ // SSE intrinsic here would multiply and then sum
            //_m128 is four floating point numbers (4x32 = 128 bits) => 16 bytes total
            // one byte is equivalent to eight bits
            v0 = simde_mm_load_ps(&squared_amplitudes[i]); // loading 4 floats at a time or 16 bytes total
//            v1 = simde_mm_load_ps(&debye_factors[next_i]);
            v1 = simde_mm_load_ps(pDF);
            mmSum1 = simde_mm_add_ps(mmSum1, simde_mm_mul_ps(v0,v1));

            v2 = simde_mm_load_ps(&squared_amplitudes[i + 4]); // loading next 4 floats for total 8 per cycle
            //v3 = simde_mm_load_ps(&debye_factors[next_i + 4]);
            pDF += 4;
            v3 = simde_mm_load_ps(pDF);
            mmSum2 = simde_mm_add_ps(mmSum2, simde_mm_mul_ps(v2,v3));
//            next_i += window;
            pDF += 4;
        }

        mmSum1 = simde_mm_add_ps(mmSum1, mmSum2);;
        simde_mm_store_ps(icalc_sse, mmSum1);
//        simde_mm256_store_ps(icalc_sse, mmSum1);
        intensity2 = icalc_sse[0] + icalc_sse[1] + icalc_sse[2] + icalc_sse[3];// + icalc_sse[4] + icalc_sse[5] + icalc_sse[6] + icalc_sse[7];

        // add remaining bits since we can't do in full width of 8
        pDF = &debye_factors[locale];
        for (unsigned int i=start_of_tail; i < total_width; i++ ){ // SSE intrinsic here would multiple and then sum
            //intensity2 += squared_amplitudes[i] * debye_factors[locale + (i - start_of_tail)];
            intensity2 += squared_amplitudes[i] * *(pDF);
            ++pDF;
        }

        // diagonal terms are included in intensity2, must remove their contribution
        *(pICalc) = -intensity_r_at_zero + 2.0f*intensity2;
        ++pICalc;
    }

//    for(unsigned int q_index = 0; q_index < total_q; q_index++) {
//        next_i  = q_index * total_width;
//        intensity2 = 0.0f;
//        for (unsigned int i=0; i < total_width; i++ ){ // SSE intrinsic here would multiply and then sum
//            intensity2 += squared_amplitudes[i]*debye_factors[next_i + i];
//        }
//        // diagonal terms are included in intensity2, must remove their contribution
//        float itemp = intensity_r_at_zero + 2.0*(intensity2 - intensity_r_at_zero);
//        std::cout << q_index << " :: " << iCalc[q_index] << " " << itemp << std::endl;
//        //iCalc[q_index] = itemp;
//    }
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
float DensityMapper::getChiSquare(unsigned int total, float scale, float * const pICalc, float * const pSigmaSquared,
                                  Datum * const pWorkingSet,
                                  Datum * const pWorkingSetSmoothed,
                                  float * const res){
    float diff_value;
    float chi2 = 0.0;
    float icalc, iobs;

    for(unsigned int i=0; i <total; i++){
        icalc = pICalc[i];
        iobs = pWorkingSetSmoothed[i].getI();
        diff_value = iobs - scale*icalc;
        diff_value *= diff_value;
        diff_value *= pSigmaSquared[i];
        chi2 += diff_value;

        // use the workset residuals to calculate durbin watson
        res[i] = pWorkingSet[i].getI() - scale*icalc;
        //res[i] = iobs - scale*icalc;
    }
    return chi2/(float)total;
}


void DensityMapper::writeICalc(unsigned int total, float scale,
                               float * const pICalc,
                               Datum * const pWorkingSet,
                               Datum * const pWorkingSetSmoothed,
                               std::string name){

    std::string tempHeader = "";
    char buffer[80];

    for(int i=0; i <total; i++){
        std::snprintf(buffer, 80, "%5i %.5E %.5E %.5E %.5E \n", i, pWorkingSet[i].getQ(), pWorkingSet[i].getI(), pWorkingSetSmoothed[i].getI(), (scale*pICalc[i]));
        tempHeader.append(buffer);
    }

    std::string mname = name + "_fit.txt";
    const char * outputFileName = mname.c_str();
    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);
}


float DensityMapper::calculateDurbinWatson(unsigned int total, const float * const residuals){

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

unsigned int DensityMapper::probability_count(float limit){
    float smallest = lattice_points[0].getMaxProbability();
    unsigned int total_amplitudes_from_lattice_points = lattice_points.size();
    LatticePoint * pLattice = lattice_points.data();
    float prior;
    unsigned int count=0;

    for(unsigned int i=0; i < total_amplitudes_from_lattice_points; i++){
        prior = pLattice[i].getMaxProbability();
        if (prior > limit){
            count++;
        }
    }

    return count;
}

float DensityMapper::calculate_rmsd(std::vector<float> & priors){

    unsigned int total_amplitudes_from_lattice_points = lattice_points.size();
    LatticePoint * pLattice = lattice_points.data();
    LatticePoint * pL;
    float diff, prior;
    float sum = 0.0f;

    float smallest = lattice_points[0].getMaxProbability();
    for(unsigned int i=1; i < total_amplitudes_from_lattice_points; i++){
        pL = &pLattice[i];
        prior = pL->getMaxProbability();
        if (smallest > prior){
            smallest = prior;
        }
    }

    float sum2 = 0.0f;
    for(unsigned int i=0; i < total_amplitudes_from_lattice_points; i++){
        pL = &pLattice[i];
        pL->setWeightedAmplitude();
//        pL->getMostProbableAmplitude();
        prior = priors[i];
//        diff = std::abs(prior - pL->getMaxProbability());
//        diff = std::abs(prior - pL->getAssignedAmplitude());
        diff = std::abs(prior - pL->getWeightedAmplitude());
        sum += diff;//*diff;
        sum2 += diff*diff;
    }

//    return sum/(float)total_amplitudes_from_lattice_points;
    return smallest;
    //return sqrtf(sum/(float)total_amplitudes_from_lattice_points);

}

double DensityMapper::f1(double x, double y){
    double result = 0.002;
    for(int i = -2; i <= 2; i++)
    {
        for(int j = -2; j<=2; j++)
        {
            double temp_res = 1 / (5*(i+2) + j + 3 + std::pow(x - 16*j, 6) + std::pow(y - 16*i, 6));
            result += std::abs(temp_res);
        }
    }
    return 1.0 / result;
}

void DensityMapper::openMP() {
    double x_from = 0.0;
    double x_to = 10000.0;
    double y_from = 0.0;
    double y_to = 10.0;
    double precision = 0.01;

    int threads = omp_get_max_threads();
    std::cout << " Threads " << threads << std::endl;
    double xInterval = std::abs((x_from - x_to)) / (double)threads;
    double result = 0;

    double wall0 = omp_get_wtime();

    #pragma omp parallel for
    for(int i = 0; i < threads; i++){
        double x_from_val = i * xInterval;
        double x_to_val = (i + 1) * xInterval;
        double y_from_val = y_from; // we create these variables, to avoid race condtions between different threads and moreover braking the data. Now, this is a thread-local variable.
        double y_to_val = y_to;
        double sum = 0;
        while(x_from_val <= x_to_val){
            double y0 = y_from_val;
            while(y0 <= y_to_val) {
                sum += f1((2 * x_from_val + precision) / 2, (2 * y0 + precision) / 2) * precision * precision;
                y0 += precision;
            }
            x_from_val += precision;
        }
        #pragma omp critical
        result += sum;
    }

    double wall1 = omp_get_wtime();
    std::cout << "THREADS " << omp_get_num_threads() << "elapsed time= " << wall1 - wall0 << " s :: "  << std::endl;
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
}

/**
 * from Vincent A. Cicirello
 * On the Design of an Adaptive Simulated Annealing Algorithm
 *
 * @param index
 * @param evalMax
 * @param acceptRate
 * @param temp
 * @param inv_temp
 */
void DensityMapper::updateASATemp(unsigned int index, unsigned int evalMax, float acceptRate, double &temp, double &inv_temp){

    double stepEval = (double)index/(double)evalMax;
    double lamRate=0.44;
    double complementASAAcceptanceRate = 1.0 - lamRate;
    double intASAAcceptanceRate = std::floor(1000*lamRate);
    double intComplementASAAcceptanceRate = std::floor(1000*complementASAAcceptanceRate);

    if (stepEval < 0.15) {
        //lamRate = 0.44+0.56*pow(560, -stepEval*6.666667);
        lamRate += complementASAAcceptanceRate*pow(intComplementASAAcceptanceRate, -stepEval*6.666667);
    } else if (stepEval >= 0.65){
        lamRate *= pow(intASAAcceptanceRate, -(stepEval - 0.65)*2.857142857);
    }

    temp = (acceptRate > lamRate) ? (0.999*temp) : (temp*1.001001001001);

    inv_temp = 1.0/temp;
}

void DensityMapper::getIndicesOfHCPSamplingGrid(std::vector<unsigned int> & selected_model_lattice_pts,
                                                std::vector<unsigned int> & keptHCPSamplePointsToUpdate) {

    std::set<unsigned int> hcpSamplingPointsToUpdate;
    for(auto & ind : selected_model_lattice_pts){ // for each model lattice point, get its neighboring points in sampling lattice
        auto iter = neighboringHCPPointsOfModelLattice.find(ind);
        if (iter != neighboringHCPPointsOfModelLattice.end()){
            for(auto useMe : iter->second){
                hcpSamplingPointsToUpdate.insert(useMe);
            }
        }
    }

    keptHCPSamplePointsToUpdate.resize(hcpSamplingPointsToUpdate.size());

    int ii=0;
    for(auto & ind : hcpSamplingPointsToUpdate){
        // ind is the index of the total lattice array, not in order of the pDen
        // find neighborhood that corresponds to the index
        auto it = find_if(neighborhoods.begin(), neighborhoods.end(), [&ind](Neighbors & obj) {return obj.getHCPLatticeIndex() == ind;});
        keptHCPSamplePointsToUpdate[ii] = it->getKeptIndex();
        ii++;
    }
}


void DensityMapper::updateICalc(std::vector<Density> & priorHCPdensities,
                                         const unsigned int total_q,
                                         std::vector<float> & i_calc,
                                         float * prior_hcp_electron_densities,
                                         float * updated_hcp_electron_densities){

    unsigned int totalPriorHCPDensities = priorHCPdensities.size();

    // correct each q value in icalc
    const unsigned int totalHCPInUse = keptHCPSamplingPoints.size(); // shared memory
    const unsigned int total_width = totalHCPInUse*(totalHCPInUse-1)/2 + totalHCPInUse;
    const unsigned int total_debye_factors = total_width*total_q;

    float debye_factor, intensity2;

//    const unsigned int window = 8; // window is 8 for _m128 and 16 for _m256
//    auto end_N = (unsigned int)std::floor(totalHCPInUse/window);
//    const unsigned int start_of_tail = end_N*window;
//    simde__m128 mmSum1, mmSum2, v0, v1, v2, v3;

    float * pICalc = i_calc.data();
    unsigned int next_i;
    unsigned int modified_index, locale, column_i;
//    simde__m128 mmSum1, mmSum2, v0, v1, v2, v3;//

    std::clock_t startTime;
    double runtime;
    float prior_en, new_en, cross_en, diff, temp_sum;
    float intensity_r_at_zero=0.0f;
    for(auto & pModified : priorHCPdensities) {
        prior_en = pModified.density;
        new_en = updated_hcp_electron_densities[pModified.index];
        intensity_r_at_zero += new_en*new_en - prior_en*prior_en;
    }

    startTime = std::clock();
    for(unsigned int q_index = 0; q_index < total_q; q_index++) {
        // for each q_value recalculate squared amplitudes for each density point
        float sum = 0.0f;
        next_i  = q_index * total_width;

        for(auto & pModified : priorHCPdensities) {
            modified_index = pModified.index;
            prior_en = pModified.density;
            new_en = updated_hcp_electron_densities[modified_index];
            temp_sum = 0.0f;
//            mmSum1 = simde_mm_setzero_ps();
//            mmSum2 = simde_mm_setzero_ps();
//
//            for (unsigned int i=0; i < start_of_tail; i+=window ) { // SSE intrinsic here would multiply and then sum
//                //_m128 is four floating point numbers (4x32 = 128 bits) => 16 bytes total
//                // one byte is equivalent to eight bits
//                v0 = simde_mm_load_ps(&squared_amplitudes[i]); // loading 4 floats at a time or 16 bytes total
//                v1 = simde_mm_load_ps(&debye_factors[next_i]);
//                mmSum1 = simde_mm_add_ps(mmSum1, simde_mm_mul_ps(v0, v1));
//
//                v2 = simde_mm_load_ps(&squared_amplitudes[i + 4]); // loading next 4 floats for total 8 per cycle
//                v3 = simde_mm_load_ps(&debye_factors[next_i + 4]);
//                mmSum2 = simde_mm_add_ps(mmSum2, simde_mm_mul_ps(v2, v3));
//            }

            for(unsigned int row=0; row < modified_index; row++){
                temp_sum += prior_hcp_electron_densities[row]*debye_factors[next_i + row*totalHCPInUse - (row*(row-1)/2) + modified_index - row];
            }

            column_i = next_i + modified_index*totalHCPInUse - (modified_index*(modified_index-1)/2) - modified_index;
            for(unsigned int col=modified_index+1; col < totalHCPInUse; col++){
                temp_sum += prior_hcp_electron_densities[col]*debye_factors[column_i + col];
            }
            sum += (new_en - prior_en)*temp_sum;
        }
        // correct icalc
        //*(pICalc) = intensity_r_at_zero + 2*(intensity2 - intensity_r_at_zero);
        i_calc[q_index] += intensity_r_at_zero + 2.0*sum;
    }

//    runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//    logger("per q TIME", formatNumber((float)runtime,8));
}


float DensityMapper::convolutionFunction(float length){
    //return (cutoff - length)/cutoff; // linear model

    if (length <= upper_limit){ // if bead overlaps with atom, keep it -> cutoff = PI/qmax
        /*
         * inv_sigma = 1.0f/(2.0f*fwhm_sigma*fwhm_sigma);
         * fwhm_sigma = 1.5f*cutoff/2.355f;
         */
        return 1.0f/std::pow(length, 1.2);
        //return (cutoff - length)/cutoff;
    } else { // sets neighbors to include in calculation
        return 0;
    }

   // return expf(-0.5f*length*length/(fwhm_sigma*fwhm_sigma));
}



void DensityMapper::printParameters(std::vector<float> & accept, std::vector<double> & temp, std::vector<float> & score, std::vector<float> & chis, std::vector<int> & flips){

    int total = accept.size();

    FILE * pFile;
    const char *outputFileName;
    std::string nameOf = "run_parameters.txt";
    outputFileName = nameOf.c_str() ;
    pFile = std::fopen(outputFileName, "w");

    std::string index;

    for (int i=0; i < total; i++){
        fprintf(pFile, "%i %5i %.4f %.4f %.4f %.4f\n", (i+1), flips[i], accept[i], temp[i], score[i], chis[i] );
    }

    fclose(pFile);
}

bool DensityMapper::checkConnectivity(std::set<unsigned int> & selectedIndices, std::vector<LatticePoint> & omp_lattice_points){

    /*
     * randomly grab from set
     * start adding if neighbor
     */
    unsigned int total = selectedIndices.size();
    std::set<unsigned int> tour;
    std::vector<unsigned int> toCheck(total-1);

    auto it = selectedIndices.begin();
    tour.insert(*it);
    it++;
    std::copy(it, selectedIndices.end(), toCheck.begin());
    total-=1;

    // starting with first selectedIndices, check if any of the remaining are a neighbor
    // if true, add to list, and check again
    bool added = true;

    while(added){

        added = false;
        for(it = tour.begin(); it != tour.end(); it++){
            LatticePoint * pLat = &omp_lattice_points[*it];

            for(unsigned int j=0; j<total; j++){ // check if neighbor is present in remaining
                unsigned int neighbor = toCheck[j];

                if (pLat->isNeighbor(neighbor)) { // if true, add to Euler Tour
                    //swap to
                    tour.insert(neighbor);
                    // remove from to check
                    std::iter_swap(toCheck.begin() + j, toCheck.begin() + total - 1);
                    toCheck.pop_back();
                    total--;
                    added = true;
                    break;
                } // else check the next one
            }

            if (added){
                break;
            }

            // if I make it all the way through with no neighbors, remove selected so we don't keep checking the remaining
        }
    }

    // if toCheck has any points remaining, suggests more than one tour
    return toCheck.empty(); // if empty means one tour
}


void DensityMapper::setLatticePointNeighbors(){

    unsigned int total = centered_coordinates.size();
    float dis;
    float min = FLT_MAX;
    // get closest packing distance from input
    for(unsigned int i = 0; i<total; i++){
        vector3 * pVec = &centered_coordinates[i];
        for(unsigned int j = i+1; j<total; j++){
            vector3 * pVec2 = &centered_coordinates[j];
            dis = (*pVec - *pVec2).sqlength();
            if (dis < min){
                min = dis;
            }
        }
    }

    min *= 1.05;

    for(unsigned int i = 0; i<total; i++){
        vector3 * pVec = &centered_coordinates[i];
        LatticePoint & pLat = lattice_points[i];

        // find all neighbors within minimum distance
        for(unsigned int j = 0; j<i; j++){
            vector3 * pVec2 = &centered_coordinates[j];
            dis = (*pVec - *pVec2).sqlength();
            if (dis < min){
                pLat.addNeighbor(j);
            }
        }

        for(unsigned int j = (i+1); j<total; j++){
            vector3 * pVec2 = &centered_coordinates[j];
            dis = (*pVec - *pVec2).sqlength();
            if (dis < min){
                pLat.addNeighbor(j);
            }
        }
    }

//    for(auto & pLat : lattice_points){
//        std::cout << pLat.getIndex() << " neighbors => " << pLat.getTotalNeighbors() << std::endl;
//    }
}