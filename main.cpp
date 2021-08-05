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
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <sys/stat.h>
#include <regex>
#include "src/utils/Example.h"
#include "src/PointSetModel.h"
#include "src/KDE.h"
#include "src/Objective.h"
#include "src/DensityMapper.h"

using namespace std::regex_constants;
namespace po = boost::program_options;

namespace {
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;
}

bool fileExists(const std::string& file) {
    struct stat buf;
    return (stat(file.c_str(), &buf) == 0);
}

int main(int argc, char** argv) {
    std::string datFile;
    std::string multiphaseFile;
    std::string chemicalType, example;
    std::vector<std::string> pdbFiles;
    std::vector<std::string> datFiles;

    float bead_radius;

    unsigned int totalSteps = 10000;

    //float volume;
    float highT, alpha=0.43, beta=0.1;
    float percentAddRemove = 0.931, lambda = 0.1, stov=0.001;//0.000051;

    float acceptanceRate = 0.44, eta = 0.001, mu = 0.0001f; // 10^-4 to 10^-5 seems to be acceptable
    unsigned int highTempRounds;
    float xaxis, yaxis, zaxis, cheight=0, cradius=0, toppercent=0.00875;
    unsigned int multiple = 37, trials=21371;

    unsigned int totalModels=1, totalPhasesForSeeded=1;

    std::string sym, prefix, maskFile, seedFile, phaseFile, hFile, anchorFile, consFile, fileList, histoFile;
    bool reduce=false, cemap=false, kdemap=false, makeSeed=false, mask=false, useCylinder=false, histo=false;
    bool useHelical=false, axesSet=false, isSeeded = false, fast=false, refine=false, useUnSmoothed=false;

    std::string descText =
            "\n              USAGE : sugatama myRefinedScatterData_pr.dat\n";
    descText += "\n       FOR SYMMETRY : -o D2 --beta 0.01";
    descText += "\n     FOR REFINEMENT : use either damstart.pdb or damfilt.pdb some other averaged aligned set";
    descText += "\n     FOR REFINEMENT : --refine --seed infile.pdb --alpha 0.003 -r 0.01 -c 3 -p refined";
    descText += "\n ADAPTIVE SIMULATED ANNEALING (ASA)";
    descText += "\n ====>  length of ASA is set using ccmultiple";
    descText += "\n ====>  acceptanceRate controls the starting temperature of the ASA run\n";
    descText += "\n ====>  default is 0.44\n";
    descText += "\n CONSTRAINTS";
    descText += "\n ====>  alpha : (only for contacts) weight for Contacts Potential";
    descText += "\n ====>   beta : Weight volume constraint for defined component or violations";
    descText += "\n ====>        : with symmetry, beta should be around 0.1 ";
    descText += "\n ====>    eta : ";
    descText += "\n ====>     mu : controls strength of volume constraint for the CVX hull in initial search";
    descText += "\n ====>   \n\n";
    descText += "\n SEQUENTIAL RUNS IN BASH SHELL";
    descText += "\n ====>  for i in {0..11} ; do iketama file_pr.dat -p run_${i} ; done\n\n";
    descText += "\n ====>  ";
    descText += "\n MODELS WITH SYMMETRY";
    descText += "\n ====>  increase mu to 0.01 using --mu 0.01 and set -u to a higher value";
    descText += "\n ====>  ";
    descText += "\n HIGH RESOLUTION MODELING AND REFINEMENT";
    descText += "\n ====>  start with a low resolution model where q_max of data is restricted to 0.2";
    descText += "\n ====>  perform 7 to 13 runs, average and use this average as a seed for high resolution";
    descText += "\n ====>  ";
    descText += "\n ====>  hires_data_pr.dat --seed aligned_set.pdb --mask";
    descText += "\n ====>  ";
    descText += "\n MAP GENERATION AND REFINEMENT";
    descText += "\n ====>  ";
    descText += "\n ====>  kde.inp is a list of pdb files that are aligned to a reference";
    descText += "\n ====>  ** REFERENCE pdb model must be the first file in the list";
    descText += "\n ====>  ";
    descText += "\n ====>  file_pr.dat --cemap --list kde.inp  --toppercent 0.001";
    descText += "\n ====>  ";
    descText += "\n ====>  with SYMMETRY";
    descText += "\n ====>  file_pr.dat -o D7 --cemap --list kde.inp --toppercent 0.001";
    descText += "\n ====>  ";
    descText += "\n ====>  If a refined model is made from the aligned set,";
    descText += "\n ====>  the seed can be used to seed the map";
    descText += "\n ====>  The seed will bias the map";
    descText += "\n ====>  ** SEED model must be aligned to the reference model in kde.inp";
    descText += "\n ====>  ";
    descText += "\n ====>  file_pr.dat --cemap --list kde.inp --seed refined_model.pdb";
    descText += "\n ====>  ";
    descText += "\n REFINEMENT STRATEGIES WITH SYMMETRY";
    descText += "\n ====>  After a single modeling run, the final subunit model can be used as input using :";
    descText += "\n ====>  ikeTama pr_file.dat -o D7 --refine --seed subunit_annealed.pdb --alpha 0.003 -c 3 -p run\n";
    descText += "\n ====>  Alternatively, a low temp run can be performed with -r as :";
    descText += "\n ====>  iketama pr_file.dat -o D7 --refine --seed damstart_subunit.pdb --alpha 0.003 -r 0.01 -c 3 -p refined\n";
    descText += "\n MASKS";
    descText += "\n ====>  A mask can be made that restricts the search space.  roku9 will create an aligned_set.pdb which is";
    descText += "\n ====>  the entire set of aligned models, this can be read in to subselect the lattice and restrict the search";
    descText += "\n";
    po::options_description desc(descText);

    desc.add_options()
            ("help,h", "Print help messages")
            ("dat", po::value<std::vector<std::string> >(&datFiles), "ScAtter *.dat files from P(r) refinement (_sx and _pr)")
            ("sym,o", po::value<std::string>(&sym)->default_value("C1"), "Sets symmetry operator for single model run")
            ("ccmultiple,c", po::value<unsigned int>(&multiple)->default_value(37), "Multiple of the Coupon Collector Sampling, increase u to increase number of steps in ASA")
            //("totalModels,d", po::value<unsigned int>(&totalModels)->default_value(1))
            ("highTempForSearch", po::value<float>(&highT)->default_value(0.0000027)) //0.0001f
            ("totalCoolingSteps", po::value<unsigned int>(&totalSteps), "Default is 10000 steps")
            ("highTempRounds,g", po::value<unsigned int>(&highTempRounds)->default_value(41269))
            ("percentAddRemove", po::value<float>(&percentAddRemove), "Sets probability of Add/Remove versus positional refinement")
            ("alpha", po::value<float>(&alpha)->default_value(alpha), "weight Contacts Distribution")
            ("beta", po::value<float>(&beta)->default_value(beta), "weight to minimize symmetry violations ")
            ("eta,e", po::value<float>(&eta)->default_value(eta), "compactness weight")
            ("lambda,l", po::value<float>(&lambda), "connectivity weight, default is 10^-2")
            ("mu,m", po::value<float>(&mu)->default_value(mu), "convex hull weight, percentage of KL divergence")
            ("components,a", po::value<std::string>(&multiphaseFile), "File describing mapping of components with dat files")
            ("nametag,p", po::value<std::string>(&prefix)->default_value("run"), "Name to tag output files")
            ("seed", po::value<std::string>(&seedFile), "Use specified input PDB as seed")
            ("mask", po::value<std::string>(&maskFile), "Use to mask the search universe")
            ("totalPhasesForSeeded", po::value<unsigned int>(&totalPhasesForSeeded), "Total number of unique, interconnected phases to be modeled")
            ("phaseFile", po::value<std::string>(&phaseFile), "File that maps data to components")
            ("refine,f", po::bool_switch(&refine), "Refine input PDB model, file is specified using seed flag")
            ("anchor", po::value<std::string>(&anchorFile), "txt file containing description of anchors")
            ("acceptanceRate,r", po::value<float>(&acceptanceRate), "ASA acceptance rate, use a rate < 0.1 for refinement")
            ("consFile", po::value<std::string>(&consFile), "Constraint File for use with --sym Xn")
            ("updateCD", po::value<std::string>(&histoFile), "Histogram of the Contacts Distribution")
            ("xaxis,x", po::value<float>(&xaxis), "length of x-axis, Angstroms (integer)")
            ("yaxis,y", po::value<float>(&yaxis), "length of y-axis, Angstroms (integer)")
            ("zaxis,z", po::value<float>(&zaxis), "length of z-axis, Angstroms (integer)")
            ("cylinder", po::value<bool>(&useCylinder), "create cylindrical search space")
            ("height", po::value<float>(&cheight), "height of cylinder")
            ("radius", po::value<float>(&cradius), "radius of cylinder")
            ("hFile", po::value<std::string>(&hFile), "file with helical parameters")
            ("direct,d", po::bool_switch(&fast), "direct")
            ("reduce", po::bool_switch(&reduce), "Reduce bead radius defining the lattice")

            ("list", po::value<std::string>(&fileList), "List of files to make estimated map")
            ("KDEmap", po::bool_switch(&kdemap), "List of files to make estimated map")
            ("cemap", po::bool_switch(&cemap), "CROSS-ENTROPY MAP from list of aligned models")
            ("makeSeed", po::bool_switch(&makeSeed), "create seed model from list of aligned models")
            ("useUnSmoothed", po::bool_switch(&useUnSmoothed), "Use Smoothed IofQ data (default)")
            ("toppercent", po::value<float>(&toppercent), "TOP PERCENT of trials to select per round for CE optimization")
            ("example", po::value<std::string>(&example), "anchor, helical")
            ("histogram", po::bool_switch(&histo), "")
            ;

    po::positional_options_description positionalOptions;
    positionalOptions.add("dat", 1);
    po::variables_map vm;

    try {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(positionalOptions).run(), vm);

        // checking options list
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return SUCCESS;
        }


        if (vm["dat"].empty()){
            std::cout << "**********************************" << std::endl;
            std::cout << " => NO DAT FILE LOADED or DETECTED" << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

        // print example input files
        if (vm.count("example")){
            std::cout << " " << vm.count("example") << std::endl;
            SASTOOLS_UTILS_H::logger("WRITING EXAMPLE INPUT", vm["example"].as<std::string >());
            Example examplar(vm["example"].as<std::string >());
            exit(0);
        }

        // can't use if seeded, setups a rectangular box
        if (vm.count("xaxis") && vm.count("yaxis") && vm.count("zaxis")){
            axesSet = true;
        }

        if (eta < 0 || lambda < 0 || mu < 0 || acceptanceRate < 0 || stov < 0 || percentAddRemove < 0){
            throw std::invalid_argument("  ERROR => INVALID SIGNED VALUE :");
        }

        if (vm.count("seed")){
            isSeeded = true;
            SASTOOLS_UTILS_H::logger("=> ", "SEEDED BUILD");
            SASTOOLS_UTILS_H::logger("INPUT PDB MODEL (SEED)", seedFile);
        }

        /*
         * mask is used for refining, reads in model and defines a limited search space
         */
        if (vm.count("mask")){
            if (!fileExists(maskFile)){ // check file exists
                // must include
                // volume of subunit
                // pitch (estimated from no symmetry model)
                // subunits per turn
                // width fiber (cylinder)
                // additonal symmetry : default C1
                logger("*** ERROR CANNOT READ FILE", maskFile);
                return 0;
            }
            mask = true;
        }

        /*
         * search space for elongated objects, should run after doing a spherical run to determine search
         * space parameters
         */
        if (vm.count("cylinder")){
            if (useCylinder){
                axesSet = false;
                if (cheight < 10 || cradius < 10){
                    throw std::invalid_argument("  ERROR => height and radius of cylinder too small");
                }
            }
        }


        /*
         * Kernel Density Approximations to aligned set of models
         */
        if (vm["KDEmap"].as<bool>()){
            std::cout << "  => creating KDE map " << std::endl;
            std::vector<std::string> datfiles = vm["dat"].as<std::vector<std::string> >();
            PofRData data(datfiles[0], false);

            bead_radius = (float)(0.4999999f * (data.getBinWidth())); // half bin width

            if (toppercent < 0.001){
                throw std::invalid_argument("** ERROR => TOPPERCENT : TOO SMALL MUST BE GREATER THAN 0.001 : " + std::to_string(toppercent) );
            }

            KDE kde(fileList, bead_radius, toppercent);
            std::cout << " BOUNDING BOX X " << kde.getCenteredMaxX() << std::endl;
            std::cout << " BOUNDING BOX Y " << kde.getCenteredMaxY() << std::endl;
            std::cout << " BOUNDING BOX Z " << kde.getCenteredMaxZ() << std::endl;

            PointSetModel model(bead_radius, 2*kde.getCenteredMaxX(), 2*kde.getCenteredMaxY(), 2*kde.getCenteredMaxZ());

            kde.createKDE(&model);

            return SUCCESS;
        }


        /*
         * make seed from list of input files (akin to damstart.pdb)
         */
        if (vm["makeSeed"].as<bool>()){
            std::cout << "  => creating SEED FROM LIST OF FILES " << std::endl;
            std::vector<std::string> datfiles = vm["dat"].as<std::vector<std::string> >();
            PofRData data(datfiles[0], false);
            bead_radius = (float)(0.4999999999 * (data.getBinWidth()));
            KDE kde(fileList, bead_radius, 0.1);
            PointSetModel model(kde.getAverageDmin(), 2*kde.getCenteredMaxX(), 2*kde.getCenteredMaxY(), 2*kde.getCenteredMaxZ());
            kde.generateSeedFromInputFile(&model, &data, prefix);
            return SUCCESS;
        }

        /*
         * Cross-Entropy Refinement and Kernel Density Approximations to aligned set of models
         */
        if (vm["cemap"].as<bool>()){

            std::cout << "  => creating KDE CE map " << std::endl;
            std::vector<std::string> datfiles = vm["dat"].as<std::vector<std::string> >();

            IofQData iofqdata = IofQData(datfiles[0], false);
            iofqdata.extractData();
            iofqdata.makeWorkingSet();
            auto workingset = iofqdata.getWorkingSet();
            auto workingsetSmoothed = iofqdata.getWorkingSetSmoothed();
            float qmax = iofqdata.getQmax();

            DensityMapper dm(maskFile, qmax, 2.5);

            // maxBin should be adjusted to 19 for ribo data, 150 is too small for dmax
            //dm.setBessels(iofqdata_bsa.getQvalues());
            //dm.calculateDensityCoefficientAtLMR();
            //int totalroounds = dm.getTotalCenteredCoordinates()*10/0.01;

//            dm.openMP();

    dm.refineModel(50, toppercent, highTempRounds,
                   const_cast<std::vector<Datum> &>(iofqdata.getWorkingSet()),
                   const_cast<std::vector<Datum> &>(iofqdata.getWorkingSetSmoothed()));


//            auto qvalues = iofqdata_bsa.getQvalues();
//            int total_in = qvalues.size();
//
//            for(int i=0; i<total_in; i++){
//                std::cout << i << " " << qvalues[i] << " " << dm.calculateIntensityAtQ(i) << std::endl;
//            }


return 1;

            PofRData data(datfiles[0], false);
            bead_radius = (float)(0.4999999999 * (data.getBinWidth()));

            if (toppercent < 0.001){
                throw std::invalid_argument("** ERROR => TOPPERCENT : TOO SMALL MUST BE GREATER THAN 0.001 : " + std::to_string(toppercent) );
            }

            KDE kde(fileList, bead_radius, toppercent);

            logger("Distribution-based lattice ", std::to_string(bead_radius));

            if (sym.compare("C1") != 0) {
                if (std::regex_match(sym, std::regex("(C|D)[0-9]+", ECMAScript | icase))) {
                    // make an offset bound box based
                    PointSetModel model(kde.getAverageDmin(), kde.getMinx(), kde.getMiny(), kde.getMinz(), kde.getMaxx(), kde.getMaxy(), kde.getMaxz(), sym);
//                    if (isSeeded){ // if refined model is available, use to bias initial map
//                        kde.add_prior(seedFile);
//                    }
                    kde.map_refineSym(&model, &data, prefix);
//                    if (vm.count("updateCD")){
//                        mainAnneal.updateContactsDistributionToFile(vm["updateCD"].as<std::string>());
//                    }
                    return SUCCESS;
                } else if (std::regex_match(sym, std::regex("X[0-9]+", ECMAScript | icase))){

                }

            } else {
                // create my universe using a spacing based on
                //float kde_bead_radius = (bead_radius - kde.getDminStDev());
                float kde_bead_radius = 0.75f*bead_radius;

                PointSetModel model(kde_bead_radius, 2*kde.getCenteredMaxX(), 2*kde.getCenteredMaxY(), 2*kde.getCenteredMaxZ());
                if (isSeeded){ // if refined model is available, use to bias initial map
                    kde.add_prior(seedFile);
                }
                kde.map_refine(&model, &data, prefix);

                return SUCCESS;
            }
        }


        Objective minFunction;

        // load data file into Objective Function
        datFiles = vm["dat"].as<std::vector<std::string> >();
        if (datFiles.size() == 1) { // single phase modeling or refinement
            /*
             * should read in a cif file that is complete with Intensity data and P(r)-distribution
             * if user wants to model reciprocal space data, then program requires dmax estimate, qmax
             *
             * Reciprocal Space Modeling
             * 1. bin_width = qmax/PI
             * 2. dmax search space or box dimensions
             *
             * Real Space Modeling
             *
             *
             */
            minFunction.addDataObject(datFiles[0]);
            /*
             * two paths to take, modeling reciprocal- or real-space data
             */
            PofRData * mainDataset = dynamic_cast<PofRData*>(minFunction.getMainDataset());

            bead_radius = reduce ? (float)(mainDataset->getBinWidth() * std::sqrt(3) / 6.0d) : (float)(mainDataset->getBinWidth()/2.0f);
            float interconnectivityCutOff = bead_radius * 2.001f;
            /*
             * SQUARE LATTICE
             */
//            bead_radius = mainDataset->getBinWidth()/(2*std::sqrt(2));
//            interconnectivityCutOff = bead_radius * 2.001f * std::sqrt(2);
            std::cout << "  =>        BEAD RADIUS :  " << bead_radius << std::endl;
            std::cout << "  =>          BIN WIDTH :  " << mainDataset->getBinWidth() << std::endl;

            /*
             * check that the specified axes make sense
             */
            if (axesSet) {
                if (xaxis / bead_radius < 10) {
                    throw "x-axis too small, increase (UNITS ANGSTROMS)";
                }
                if (yaxis / bead_radius < 10) {
                    throw "y-axis too small, increase (UNITS ANGSTROMS)";
                }
                if (zaxis / bead_radius < 10) {
                    throw "z-axis too small, increase (UNITS ANGSTROMS)";
                }
            }

            /*
             * for tetrahedron, radius of circumscribing sphere is 1.22*d_max
             * As a guess, 1.3 should work in most cases;
             */
            float searchSpace = (isSeeded && !refine) ? (mainDataset->getDmax() * 1.47f) : mainDataset->getDmax()*1.4f;

            /*
             * create Search Space Universe : PointSetModel
             */
            PointSetModel model;
            //model = PointSetModel(maskFile, searchSpace, bead_radius); // use this for high resolution modeling

            if (sym.compare("C1") != 0) {

                if (std::regex_match(sym, std::regex("(C|D)[0-9]+", ECMAScript | icase)) ||
                    std::regex_match(sym, std::regex("X[0-9]+", ECMAScript | icase))) {

                    if (axesSet) {
                        model = PointSetModel(bead_radius, xaxis, yaxis, zaxis, sym);
                    } else if (useCylinder) {
                        model = PointSetModel(bead_radius, cheight, cradius, sym);
                    } else if (vm.count("seed") && (mask || refine)) {
                        /*
                         * use the PDBs in seedFile to mask the universe
                         */
                        model = PointSetModel(maskFile, searchSpace, bead_radius, sym, mask); // use this for high resolution modeling
                    } else {
                        searchSpace = mainDataset->getDmax()*1.3f;
                        model = PointSetModel(searchSpace, bead_radius, sym);
                    }

                } else if (useHelical) { // unrecognized symmetry
                    model = PointSetModel(hFile, bead_radius);
                }

            } else if (sym == "C1" || sym == "c1") { // no symmetry

                if (axesSet) {
                    model = PointSetModel(bead_radius, xaxis, yaxis, zaxis);
                } else if (useCylinder) {
                    model = PointSetModel(bead_radius, cheight, cradius);
                } else if (mask) {
                    model = PointSetModel(maskFile, searchSpace, bead_radius); // use this for high resolution modeling
                } else {
                    model = PointSetModel(searchSpace, bead_radius);
                }
            }


            /*
             * SetUp Annealer
             */
            Anneal mainAnneal(highT,
                              percentAddRemove,
                              highTempRounds,
                              prefix,
                              alpha,
                              beta,
                              eta,
                              lambda,
                              mu,
                              multiple,
                              acceptanceRate,
                              interconnectivityCutOff
            );


            /*
             *
             * de novo lattice modeling has two main approaches
             * Fastest optimized method creates the search space and precalculates the entire distance space
             * The fastest method will exceed the memory capactiy of the CPU
             * If too large, we switch to slow mode (Direct Method)
             */
            std::string nameTo = "initial_" + prefix;

            /*
             * Following is for Direct Method
             * Does not precalculate the distances in the Universe
             * Should not be limited by size
             *
             */
            if (fast || (model.isUseDirectMethod() && !(std::regex_match(sym, std::regex("(C|D)[0-9]+", ECMAScript | icase)) &&
                                                        sym.compare("C1") !=
                                                        0))) { // set to true if too many beads in search space

                /*
                 *
                 * 1. refinement
                 *     - input PDB model is searched at low temp to produce a better model
                 * 2. seeded
                 *     - input PDB model is used as a scaffold for building new model
                 * 3. de novo
                 *     - no starting model, just build from scratch
                 */
                if (isSeeded && refine) { // direct method refine model
                    /*
                    * align all models using a program like damaver -a *.pdb
                    * realign damfilt or damstart to reference model if necessary
                    * Easy check is to open damfilt against the reference, should superimpose
                    * use this damfilt/damstart model as input for refined.
                    */
//                    if (!mainAnneal.initializeModelToRefine(&model, mainDataset, prefix, seedFile)) {
//                        return ERROR_UNHANDLED_EXCEPTION;
//                    }
//                    mainAnneal.refineHomogenousBodyASAHybridDirect(&model, mainDataset, prefix);

                } else if (isSeeded && !refine) { // SEEDED


                } else { // DE NOVO
                    if (!mainAnneal.createInitialModelCVXHullDirect(&model, mainDataset, nameTo)) {
                        return ERROR_UNHANDLED_EXCEPTION;
                    }
                    mainAnneal.refineHomogenousBodyASAHybridDirect(&model, mainDataset, prefix);
                }

            } else if (std::regex_match(sym, std::regex("X[0-9]+", ECMAScript | icase))) { // IDENTICAL SUBUNITS

                if ((refine && isSeeded)) {
                    /*
                     * align all models using a program like damaver -a *.pdb
                     * realign damfilt or damstart to reference model if necessary
                     * Easy check is to open damfilt against the reference, should superimpose
                     * use this damfilt/damstart model as input for refined.
                     *
                     * output file should have rotation matrix, cofm point for each subunit
                     *
                     */

                } else {

//                    if (!mainAnneal.createInitialModelSymmetryX(&model, mainDataset)) {
//                        return ERROR_UNHANDLED_EXCEPTION;
//                    }
//
//                    mainAnneal.estimateInterSubUnitDistancesX(&model, mainDataset);
//                    mainAnneal.refineHomogenousBodyASAHybridSymmetryX(&model, mainDataset, prefix);

                }

            }  else if (useHelical) {

//                if (!mainAnneal.createInitialModelHelicalSymmetry(&model, mainDataset)) {
//                    return ERROR_UNHANDLED_EXCEPTION;
//                }
//
//                mainAnneal.refineSymModelHelical(&model, mainDataset, prefix);


            } else { // This section uses PreCalculated Distances from Universe

                if (std::regex_match(sym, std::regex("(C|D)[0-9]+", ECMAScript | icase )) && sym.compare("C1") != 0 ) { // SYMMETRY CONSTRAINT

                    logger("SYMMETRY SET", sym);

                    /*
                     *
                     * 1. refinement
                     *     - no starting model, just build from scratch
                     * 2. de novo
                     *     - input PDB model is used as a scaffold for building new model
                     *
                     */
                    if ((refine && isSeeded)) { // with symmetry

                        /*
                         * align all models using a program like roku9 or damaver -a *.pdb
                         * realign damfilt or damstart to reference model if necessary
                         * Easy check is to open damfilt against the reference, should superimpose
                         * use this damfilt/damstart model as input for refined.
                         *
                         * MASKED search performs a normal search but within a highly defined search space
                         *
                         * Refinement takes a model and performs a low energy optimization since we are
                         * assuming we a close to the correct answer
                         *
                         */

                        if (!mainAnneal.initializeModelToRefineSym(&model, mainDataset, prefix, seedFile)) {
                            return ERROR_UNHANDLED_EXCEPTION;
                        }

                        // make masked universe from input seed and run refineSymModel
                        mainAnneal.refineSymModelRefine(&model, mainDataset, prefix);

                    } else { // masked Universe or not will be the set up for the following
                      
                        if (!mainAnneal.createInitialModelSymmetry(&model, mainDataset)) {
                            return ERROR_UNHANDLED_EXCEPTION;
                        }

                        mainAnneal.refineSymModel(&model, mainDataset, prefix);
                    }

                } else { // no symmetry
                    /*
                     *
                     * 1. refinement
                     *     - input PDB model is searched at low temp to produce a better model
                     * 2. seeded
                     *     - input PDB model is used as a scaffold for building new model
                     * 3. de novo
                     *     - no starting model, just build from scratch
                     */
                    if (isSeeded && refine) {

                        /*
                        * align all models using a program like damaver -a *.pdb
                        * realign damfilt or damstart to reference model if necessary
                        * Easy check is to open damfilt against the reference, should superimpose
                        * use this damfilt/damstart model as input for refined.
                        */
//                        if (!mainAnneal.initializeModelToRefine(&model, mainDataset, prefix, seedFile)){
//                            return ERROR_UNHANDLED_EXCEPTION;
//                        }
//
//                        mainAnneal.refineHomogenousBodyInputEx(&model, mainDataset, prefix);

                    } else if (isSeeded && !refine && !mask) { // SEEDED
                        /*
                         * use seeded if you have a starting atomistic model and would like to add to it
                         * for instance, structure of A and complex of AB, assuming in the bound state A is the same
                         * or if a domain is missing in A' and SAXS is of A (A' + missing domain)
                         */

                        if (fileExists(anchorFile)){ // check that anchor points exists in seedFile
                            std::cout << "*** READING ANCHOR FILE ***" << std::endl;
                            mainAnneal.setAnchorPoints(anchorFile, seedFile, &model);
                        }

                        std::cout << "*** CREATING INITIAL MODEL FROM SEED ***" << std::endl;
                        mainAnneal.createSeedFromPDB(&model, mainDataset, "reduced_seed", seedFile);
exit(0);
                        mainAnneal.refineHomogenousBodyASACVXSeeded(&model, mainDataset, prefix);

                    } else { // DE NOVO

                        if (totalModels == 1){

                            if (!mainAnneal.createInitialModelCVXHull(&model, mainDataset, nameTo)){
                                return ERROR_UNHANDLED_EXCEPTION;
                            }
                            mainAnneal.refineHomogenousBodyASAHybridEx(&model, mainDataset, prefix);
                        }
                    }
                }


            }


        } else if (vm.count("components") == 1) {

            // read in file

            // validate multiphase file format
            return 1;

        }


    } catch (boost::program_options::required_option& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return ERROR_IN_COMMAND_LINE;
    } catch (boost::program_options::error& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return ERROR_IN_COMMAND_LINE;
    } catch (const std::invalid_argument& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr<<"Type "<<typeid(e).name()<<std::endl;
    }


    return 0;
}