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
#include <sastools/PDBModel.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>

DensityMapper::DensityMapper(std::string filename) : inputFilename(std::move(filename)) {
    PDBModel model(inputFilename, true, false);
    const unsigned int total = model.getTotalCoordinates();
    if (total == 0) throw std::invalid_argument("DensityMapper: no coordinates in " + inputFilename);
    centered_coordinates.reserve(total);
    for (unsigned int i = 0; i < total; i++)
        centered_coordinates.emplace_back(model.getCenteredXVec()[i], model.getCenteredYVec()[i], model.getCenteredZVec()[i]);
    model.writeCenteredCoordinatesToFile("CenteredInputModel");   // the frame of the refined maps
}

/*
 * bead radius from the "EDGE RADIUS : r" remark that Kaiten/IKETAMA write into model files (0 if absent)
 */
float DensityMapper::readEdgeRadius(const std::string & pdbfile){
    std::ifstream in(pdbfile); std::string line;
    while (std::getline(in, line)){
        if (line.compare(0, 6, "ATOM  ") == 0) break;
        auto pos = line.find("EDGE RADIUS");
        if (pos == std::string::npos) continue;
        auto colon = line.find(':', pos);
        if (colon == std::string::npos) continue;
        try { return std::stof(line.substr(colon + 1)); } catch (...) { return 0.0f; }
    }
    return 0.0f;
}

/*
 * See DensityMapper.h and DensityRefiner.h.
 */
densrefine::LambdaResult DensityMapper::refineModelDensity(std::vector<Datum> & workingSet, PofRData * pofrData,
                                                           densrefine::Options options, const std::string & prefix){
    if (options.spacing <= 0){
        float r = readEdgeRadius(inputFilename);
        if (r <= 0) throw std::invalid_argument("refineModelDensity: lattice spacing not set and no EDGE RADIUS in " + inputFilename);
        options.spacing = 2.0*r;
    }

    std::vector<densrefine::V3> beads; beads.reserve(centered_coordinates.size());
    for (auto & v : centered_coordinates) beads.push_back({v.x, v.y, v.z});

    std::vector<double> q, iofq, sigma;
    for (auto & d : workingSet){ q.push_back(d.getQ()); iofq.push_back(d.getI()); sigma.push_back(d.getSigma()); }
    // dmax for the Shannon number: from the P(r) file, else the largest distance in the aligned set plus a bead
    double dmax = (pofrData != nullptr) ? pofrData->getDmax() : 0.0;
    if (dmax <= 0){
        double r2 = 0; for (auto & v : centered_coordinates) r2 = std::max(r2, (double)(v.x*v.x + v.y*v.y + v.z*v.z));
        dmax = 2.0*std::sqrt(r2) + options.spacing;
    }
    densrefine::Data data = densrefine::makeData(q, iofq, sigma, densrefine::rebinForShannon(q, dmax, options.pointsPerShannon));

    if (options.snapshots > 0 && options.snapshotPrefix.empty()) options.snapshotPrefix = prefix;
    densrefine::Refiner refiner(beads, data, options);

    if (pofrData != nullptr){
        const unsigned int bins = pofrData->getShannonBins();
        std::vector<double> prob(bins);
        for (unsigned int b = 0; b < bins; b++) prob[b] = pofrData->getProbabilityPerBin(b);
        refiner.setPrTarget(prob, [pofrData](double d){ return (int)pofrData->convertToBin((float)d); });
    }

    refiner.run();
    refiner.write(prefix);

    const densrefine::LambdaResult & best = refiner.chosen();
    if (options.verbose){
        std::printf("  refineModelDensity: lambda %.1f | chi2 work %.3f free %.3f | DW %.2f | P(r) KL %.2e | Rg %.1f  (start chi2 %.3f / %.3f) => %s_map.xplor\n",
                    best.lambda, best.chiW, best.chiF, best.dw, best.kl, best.rg, refiner.startChiW(), refiner.startChiF(), prefix.c_str());
    }
    return best;
}
