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
#include <boost/regex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include "Objective.h"

using namespace boost::filesystem;

Objective::Objective(){

}

/**
 * returns a pointer to Data Object
 * holds the memory address to element in datasets vector
 *
 */
DataBase * Objective::getMainDataset() {

    DataBase * preturnMe = datasets.begin()->second;

    int phaseCount = 0;

    if (datasets.size() > 1){
        // iterate over each dataset
        // determine total number of phases in each dataset
        // identify dataset with largest number of phases
        for(auto & iterator : datasets) {
            // getData object
            DataBase * temp = iterator.second;
//            if (temp.getTotalPhases() > phaseCount){
//                preturnMe = &(iterator->second);
//                phaseCount = temp.getTotalPhases();
//            }
        }

    }

    return preturnMe;
}


/**
 *
 */
//void Objective::createWorkingSets(PointSetModel & model) {
//
//    int totalDatasets = datasets.size();
//
//    // generates different sets of q-values per dataset
//    for(int i=0; i<totalDatasets; i++){
//        this->getDataObjectByIndex(i)->creatingWorkingSet(model);
//    }
//
//}


/**
 * Create dataset object using both Intensity (Reciprocal Space) and Real-Space datasets
 * Modeling is based on real-space data.
 *
 * Why keep Reciprocal?
 *
 * @param iofqfile
 * @param pofrfile
 */
void Objective::addDataObject(std::string datFile) {
    // load file and create Data object
    // get name of iofqfile, strip away
    path p = datFile;
    boost::regex slashes("/|\\\\");

    std::string identifierKey;
    std::vector<std::string> tempLine;

    // getFiletype of datFile

    try {
        if (exists(p)){
            //split name if forward or backslash is present
            if(boost::regex_search(datFile, slashes)){

                boost::split(tempLine, datFile, boost::is_any_of("/|\\\\"), boost::token_compress_on);

                //int last = tempLine.size()-1;
                identifierKey = tempLine[tempLine.size()-1];

            } else {
                identifierKey = datFile;
            }

            if (SASTOOLS_UTILS_H::validatePofRFile(datFile)){
                datasets.emplace(identifierKey, new PofRData(datFile, false));
            } else {
                datasets.emplace(identifierKey, new IofQData(datFile, false));
            }

            keys.push_back(identifierKey);
        }

    } catch (const filesystem_error& ex) {
        std::cout << ex.what() << '\n';
    }

    /*
    cout  <<  "  root_name()----------: " << p.root_name() << '\n';
    cout  <<  "  root_directory()-----: " << p.root_directory() << '\n';
    cout  <<  "  root_path()----------: " << p.root_path() << '\n';
    cout  <<  "  relative_path()------: " << p.relative_path() << '\n';
    cout  <<  "  parent_path()--------: " << p.parent_path() << '\n';
    cout  <<  "  filename()-----------: " << p.filename() << '\n';
    */
}
