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


#include <sastools/include/utils.h>
#include <cfloat>
#include <sastools/include/PDBModel.h>
#include "PointSetModel.h"

/**
 * Create the Universe
 * @param searchSpace
 * @param beadradius
 */
PointSetModel::PointSetModel(float searchSpace, float beadradius) : radius_of_universe(searchSpace*0.5f), bead_radius(beadradius), cutOffNeighbor(2.001f*bead_radius) {

//    radial_limit = std::ceil((radius_of_universe + radius_of_universe*0.1011f)/bead_radius)*bead_radius;
    radial_limit = std::ceil((radius_of_universe)/bead_radius)*bead_radius;


    inv_bead_radius = 1.0f/bead_radius;
    bead_volume = (float)(4.0f/3.0f*bead_radius*bead_radius*bead_radius*M_PI); // 4/3*PI*r^3

    this->createUniverse(true);
    this->createDistancesAndConvertToSphericalCoordinates();
}

PointSetModel::PointSetModel(float beadradius, float xaxis, float yaxis, float zaxis) :  bead_radius(beadradius), cutOffNeighbor(2.001f*bead_radius),
                                                                         xaxis(xaxis),
                                                                         yaxis(yaxis),
                                                                         zaxis(zaxis) {

    char* text = new char[7];
    std::sprintf(text,"%.3f", beadradius);
    logger("HCP lattice spacing", text);

    inv_bead_radius = 1.0f/bead_radius;
    bead_volume = (float)(4.0f/3.0f*bead_radius*bead_radius*bead_radius*M_PI); // 4/3*PI*r^3

    if ((this->xaxis*inv_bead_radius < 4 || this->yaxis*inv_bead_radius < 4 || this->zaxis*inv_bead_radius < 4 )){
        throw std::invalid_argument(" AXIS too SMALL ");
    }

    useCylindricalSearchSpace = false;
    this->createUniverse(false);
    this->createDistancesAndConvertToSphericalCoordinates();
    delete[]text;
    //this->writeSubModelToFile(0, beads.size(), bead_indices, "universe");
}


/**
 * Create Rectangular Shaped Universe
 * @param beadradius
 * @param xaxis
 * @param yaxis
 * @param zaxis
 * @param sym
 */
PointSetModel::PointSetModel(float beadradius, float xaxis, float yaxis, float zaxis, std::string sym) : PointSetModel(beadradius, xaxis, yaxis, zaxis) {
    this->symmetry = sym;
    setSymmetryParameters(sym);
}

/*
 * Use for refining map
 * Create a box that contains aligned models (asymmetric units) from many runs
 * Translate box to center of coordinates (asymmetric unit is not centered at 0,0,0)
 * prune the box to be shape of asymmetric unit (wedge)
 */
PointSetModel::PointSetModel(float beadradius, float xmin, float ymin, float zmin, float xmax, float ymax, float zmax, std::string sym) :
        bead_radius(beadradius),
        cutOffNeighbor(2.001f*bead_radius),
        xaxis(xmax-xmin),
        yaxis(ymax-ymin),
        zaxis(zmax-zmin) {

    inv_bead_radius = 1.0f/bead_radius;
    bead_volume = (float)(4.0f/3.0f*bead_radius*bead_radius*bead_radius*M_PI); // 4/3*PI*r^3

    if ((this->xaxis*inv_bead_radius < 4 || this->yaxis*inv_bead_radius < 4 || this->zaxis*inv_bead_radius < 4 )){
        throw std::invalid_argument(" AXIS too SMALL ");
    }

    useCylindricalSearchSpace = false;
    this->symmetry = sym;
    setSymmetryParameters(sym);

    this->createUniverse(false); // create boxed universe

    vector3 translateTo = vector3(0.5f*(xmax+xmin), 0.5f*(ymax+ymin), 0.5f*(zmax+zmin));
    // translate the box
    for(auto & bead : beads){ // beads are centered at zero
        bead.translateTo(translateTo);
    }

    this->pruneUniverseToAsymmetricUnit(); // prune to asymmetric unit
    this->createDistancesAndConvertToSphericalCoordinates();
    //this->writeSubModelToFile(0, beads.size(), bead_indices, "kdeuniverse");
}

/*
 * use the mask to limit the search space
 *
 * refining is to take an existing model do an additional low-temp (low acceptance rate) search
 * so we are looking for small adjustments around the input model
 */
PointSetModel::PointSetModel(std::string maskfile, float searchSpace, float beadradius, std::string sym, bool maskit) : radius_of_universe(searchSpace*0.5f), bead_radius(beadradius), cutOffNeighbor(2.001f*bead_radius) {

    radial_limit = std::ceil((radius_of_universe + radius_of_universe*0.1011f)/bead_radius)*bead_radius;

    inv_bead_radius = 1.0f/bead_radius;
    bead_volume = (float)(4.0f/3.0f*bead_radius*bead_radius*bead_radius*M_PI); // 4/3*PI*r^3

    this->symmetry = sym;
    setSymmetryParameters(sym);

    if (maskit){
        logger("MASKING TO", maskfile);
        this->createUniverseFromMaskSym(maskfile);
    } else {
        this->createUniverse(true);
    }

    logger("","PRUNING TO ASSYMETRIC UNIT");
    this->pruneUniverseToAsymmetricUnit();

    this->createDistancesAndConvertToSphericalCoordinates();
    this->writeSubModelToFile(0, beads.size(), bead_indices, "asym_universe");
}

/**
 * Cylindrical search space with symmetry
 * @param beadradius
 * @param height
 * @param radius
 * @param sym
 */
PointSetModel::PointSetModel(float beadradius, float height, float radius, std::string sym) : bead_radius(beadradius), cutOffNeighbor(2.001f*bead_radius),yaxis(radius),
                                                                              zaxis(height)   {
    this->radial_limit = yaxis;
    this->xaxis = yaxis;
    this->symmetry = sym;
    setSymmetryParameters(sym);

    inv_bead_radius = 1.0f/bead_radius;
    bead_volume = (float)(4.0f/3.0f*bead_radius*bead_radius*bead_radius*M_PI); // 4/3*PI*r^3

    useCylindricalSearchSpace = true;

    this->createUniverse(false);
    this->createDistancesAndConvertToSphericalCoordinates();
    this->writeSubModelToFile(0, beads.size(), bead_indices, "cylinderuniverse");
}

/*
 *
 */
PointSetModel::PointSetModel(float searchSpace, float beadradius, std::string sym) : radius_of_universe(searchSpace*0.5f), bead_radius(beadradius), cutOffNeighbor(2.001f*bead_radius) {

    radial_limit = std::ceil((radius_of_universe + radius_of_universe*0.07f)/bead_radius)*bead_radius;
    inv_bead_radius = 1.0f/bead_radius;
    bead_volume = (float)(4.0f/3.0f*bead_radius*bead_radius*bead_radius*M_PI); // 4/3*PI*r^3

//    this->createRectangularUniverse(true);
    this->createUniverse(true);

    this->symmetry = sym;
    setSymmetryParameters(sym);

    this->pruneUniverseToAsymmetricUnit();

    this->createDistancesAndConvertToSphericalCoordinates();
    this->writeSubModelToFile(0, beads.size(), bead_indices, "rectangularUniverse");

}

PointSetModel::PointSetModel(float beadradius, float height, float radius) : bead_radius(beadradius), cutOffNeighbor(2.001f*bead_radius),yaxis(radius),
                                                             zaxis(height)   {

    this->radial_limit = yaxis;
    this->xaxis = yaxis;

    inv_bead_radius = 1.0f/bead_radius;
    bead_volume = (float)(4.0f/3.0f*bead_radius*bead_radius*bead_radius*M_PI); // 4/3*PI*r^3

    useCylindricalSearchSpace = true;

    this->createUniverse(false);
    this->createDistancesAndConvertToSphericalCoordinates();
    //this->writeSubModelToFile(0, beads.size(), bead_indices, "universe");
}

/**
 * Constructor for only systems with helical symmetry
 * @param helicalfile
 * @param beadradius
 */
PointSetModel::PointSetModel(std::string helicalfile, float beadradius){

    this->symmetry = "H";
    this->useCylindricalSearchSpace = true;

    this->bead_radius = beadradius;
    inv_bead_radius = 1.0f/bead_radius;
    bead_volume = (float)(4.0f/3.0f*bead_radius*bead_radius*bead_radius*M_PI); // 4/3*PI*r^3
    cutOffNeighbor = (2.001f*bead_radius);

    std::ifstream data (helicalfile.c_str());

    // screw symmetry is unknown along long dimesion
    if (data.is_open()) {

        boost::regex pitchFormat("(pitch)|(PITCH)", boost::regex::icase); // zaxis
        boost::regex widthFormat("(radius)|(RADIUS)", boost::regex::icase); // zaxis
        boost::regex riseFormat("RISE", boost::regex::icase); // rise per subunit
        boost::regex thetaFormat("theta", boost::regex::icase); // rise per subunit
        boost::regex volumeFormat("(VOLUME)|(volume)", boost::regex::icase);
        boost::regex pointGroupSymmetry("(C|D)[0-9]+", boost::regex::icase);

        //helicalVolume =
        // dmax divided by rise => total subUnits
        std::string line;

        while(!data.eof()) {
            getline(data, line);
            boost::algorithm::trim(line);
            std::vector<std::string> contents;
            boost::split(contents, line, boost::is_any_of("\t  "), boost::token_compress_on);
            std::string first = std::string(contents[0]);

            if (boost::regex_search(first, pointGroupSymmetry)){
                this->symmetry = contents[1];
            } else if (boost::regex_search(first, pitchFormat)){
                this->zaxis = std::stof(contents[1]);
            } else if (boost::regex_search(first, widthFormat)){
                this->yaxis = std::stof(contents[1]);
                this->xaxis = std::stof(contents[1]);
                std::cout << line << std::endl;
            } else if (boost::regex_search(first, riseFormat)){
                this->risePerSubUnit = std::stof(contents[1]);
            } else if (boost::regex_search(first, thetaFormat)){
                this->rotationPerSubunit = std::stof(contents[1]);
            } else if (boost::regex_search(first, volumeFormat)){
                this->helicalVolume = std::stof(contents[1]);
            }

        }

        std::cout << " HEIGHT " << this->zaxis << std::endl;
        std::cout << " XWIDTH " << this->yaxis << std::endl;
        std::cout << " YWIDTH " << this->xaxis << std::endl;
        std::cout << "    PGS " << this->symmetry << std::endl;

        setSymmetryParameters(this->symmetry);
    }
    data.close();

    this->createUniverse(false);
    this->createDistancesAndConvertToSphericalCoordinates();

    std::string name="cylinder_search_space";
    std::vector<unsigned int> bead_indices(this->number_of_beads); // large vector ~1000's
    //std::clock_t start;
    // c-style is slightly faster for large vector sizes
    // start = std::clock();
    unsigned int * ptr = (this->number_of_beads != 0) ? &bead_indices.front() : nullptr;
    for(unsigned int i = 0; i < this->number_of_beads; i++) {
        ptr[i] = i;
    }

    this->writeModelToFile(this->number_of_beads, bead_indices, name, number_of_beads);

}

PointSetModel::PointSetModel(std::string maskfile, float searchSpace, float beadradius) : radius_of_universe(searchSpace*0.5f), bead_radius(beadradius), cutOffNeighbor(2.001f*bead_radius) {

    radial_limit = std::ceil(radius_of_universe/bead_radius)*bead_radius;
//    radial_limit = std::ceil((radius_of_universe + radius_of_universe*0.1011f)/bead_radius)*bead_radius;

    inv_bead_radius = 1.0f/bead_radius;
    bead_volume = (float)(4.0f/3.0f*bead_radius*bead_radius*bead_radius*M_PI); // 4/3*PI*r^3

    this->createUniverseFromMask(maskfile);
    this->createDistancesAndConvertToSphericalCoordinates();
    this->writeSubModelToFile(0, beads.size(), bead_indices, "masked_universe");
}

void PointSetModel::setSymmetryParameters(std::string sym) {

    symmetry = sym;
    symmetryGroup = symmetry[0];
    boost::to_upper(symmetryGroup);

    logger("SYMMETRY GROUP", symmetryGroup);

    if (symmetry.substr(0,1) == "C" || symmetry.substr(0,1) == "c" ){ // rotation group
        symmetryIndex = (unsigned int) atoi(sym.erase(0,1).c_str());
        numberOfSubUnits = symmetryIndex; // number of identical subunits
        // neighboring subunit is the first index > 0
        neighboringSubunits.push_back(1);
    } else if (symmetry.substr(0,1) == "D" || symmetry.substr(0,1) == "d" ) { // dihedral group
        symmetryIndex = (unsigned int) atoi(sym.erase(0,1).c_str());
        numberOfSubUnits = 2*symmetryIndex; // number of identical subunits
        /*
         * with respect to the base subunit (index 0), neighboring subunits are the
         * symmetry related subunits that are touching the base subunit
         */
        neighboringSubunits.push_back(1);
        neighboringSubunits.push_back(symmetryIndex);
    } else if (symmetry.substr(0,1) == "O" || symmetry.substr(0,1) == "o" ) { // octahedral
        symmetryIndex = 4;
        numberOfSubUnits = 8; // number of identical subunits
        neighboringSubunits.push_back(1);
        neighboringSubunits.push_back(symmetryIndex);
    } else if (symmetry.substr(0,1) == "T" || symmetry.substr(0,1) == "t" ) { // tetrahedral
        symmetryIndex = 1;
        numberOfSubUnits = 4; // number of identical subunits
        neighboringSubunits.push_back(1);
    } else if (symmetry.substr(0,1) == "I" || symmetry.substr(0,1) == "i" ) { // icosahedral

    } else if (symmetry.substr(0,1) == "H" || symmetry.substr(0,1) == "h" ) { // helical
        // model the persistance length and propagate
        symmetryIndex = 1;
        neighboringSubunits.push_back(1);

    } else if (symmetry.substr(0,1) == "X" || symmetry.substr(0,1) == "x" ) { // helical
        // no symmetry, find X number of identical subunits
        // subunits must be touching
        symmetryIndex = (unsigned int) atoi(sym.erase(0,1).c_str());

        if (symmetryIndex < 2){
            exit(0);
        }

        numberOfSubUnits = symmetryIndex; // number of identical subunits
        neighboringSubunits.push_back(1);
    }

    logger("SYMMETRY SET", symmetry);
    logger("SYMMETRY INDEX", std::to_string(symmetryIndex));
    logger("TOTAL SUBUNITS", std::to_string(numberOfSubUnits));

    static const boost::regex bad_words("X", boost::regex::icase);
    // octahedral or icosahedral symmetries would be too large for this approach
    // calculations should be made on the fly
    if ( !(boost::regex_search(symmetry.substr(0,1), bad_words))){
        this->createSubUnitAngles();
    }
}


void PointSetModel::createSubUnitAngles(){

    subUnitAnglesCos.resize(numberOfSubUnits);
    subUnitAnglesSin.resize(numberOfSubUnits);

    subUnitAnglesCos[0]=1.0f;
    subUnitAnglesSin[0]=0.0f;

    for(unsigned int i=1; i<numberOfSubUnits; i++){
        double angle = M_PI*2.0d*i/(double)symmetryIndex;
        subUnitAnglesCos[i]=(float)std::cos(angle);
        subUnitAnglesSin[i]=(float)std::sin(angle);
    }
}


void PointSetModel::createUniverse(bool useSphericalModel) {
    beads.clear();
    std::cout << " --------------------------------------------------- " << std::endl;
    if (useSphericalModel){
        logger("=> CREATING UNIVERSE","SPHERICAL SEARCH SPACE");
        const auto z_prime_index = (int) (2*std::ceil(radial_limit*inv_bead_radius*0.5))+3;
        const auto x_prime_index = (int) (2*std::ceil(radial_limit*inv_bead_radius*0.5))+3;
        const auto y_prime_index = (int) (2*std::ceil(radial_limit*inv_bead_radius*0.5))+3;

        const vector3 a1 = vector3(2*bead_radius, 0, 0);
        const vector3 a2 = vector3(-bead_radius, bead_radius*sqrt(3.0f), 0);
        const vector3 a3 = vector3(0, 0, (float)(4.0f*sqrt(6)/3.0f*bead_radius)); // seems fine, spreads out

        // basis 1 is (0,0,0)
        const vector3 basis2 = a1*(2.0f/3.0f) + a2*(1.0f/3.0f) + a3/2.0f;

        for (int n3 = -z_prime_index; n3<z_prime_index; n3++){
            for (int n2 = -y_prime_index; n2 <= y_prime_index; n2++) {
                for (int n1 = -x_prime_index; n1 <= x_prime_index; n1++) {

                    const vector3 temp = (a1 * n1 + a2 * n2 + a3 * n3);

                    if ((temp.length()) <= radial_limit){
                        beads.emplace_back(Bead(temp.x, temp.y, temp.z, 1));
                    }

                    vector3 temp2 = temp + basis2;

                    if ((temp2.length()) <= radial_limit) {
                        beads.emplace_back(Bead(temp2.x, temp2.y, temp2.z, 1));
                    }
                }
            }
        }

    } else {

        const auto x_prime_index = (int) (2*std::ceil(xaxis*inv_bead_radius*0.5))+10;

        if (useCylindricalSearchSpace){
            logger("=> CREATING UNIVERSE","CYLINDER SEARCH SPACE");
            // klimit is height of cylinder
            const auto y_prime_index = x_prime_index;

            const vector3 a1 = vector3(2*bead_radius, 0, 0);
            const vector3 a2 = vector3(-bead_radius, bead_radius*sqrt(3.0f), 0);
            const vector3 a3 = vector3(0, 0, (float)(4.0f*sqrt(6)/3.0f*bead_radius)); // seems fine, spreads out

            const auto z_prime_index = (int) std::ceil(zaxis/a3.length()*0.5);

            // basis 1 is (0,0,0)
            const vector3 basis2 = a1*(2.0/3.0) + a2*(1.0f/3.0) + a3/2.0;

            for (int n3 = -z_prime_index; n3<z_prime_index; n3++){

                for (int n2 = -y_prime_index; n2 <= y_prime_index; n2++) {
                    for (int n1 = -x_prime_index; n1 <= x_prime_index; n1++) {
                        const vector3 temp = (a1 * n1 + a2 * n2 + a3 * n3);
                        if (std::sqrt(temp.x*temp.x + temp.y*temp.y) < radial_limit){
                            beads.emplace_back(Bead(temp.x, temp.y, temp.z, 1));
                        }

                        vector3 temp2 = temp + basis2;
                        if (std::sqrt(temp2.x*temp2.x + temp2.y*temp2.y) < radial_limit){
                            beads.emplace_back(Bead(temp2.x, temp2.y, temp2.z, 1));
                        }
                    }
                }
            }

        } else {
            logger("=> CREATING UNIVERSE","BOXED SEARCH SPACE");
            const auto z_prime_index = (int) (std::ceil(zaxis*inv_bead_radius*0.5) + 10);
            const auto y_prime_index = (int) (std::ceil(yaxis*inv_bead_radius*0.5) + 10);

            const vector3 a1 = vector3(2*bead_radius, 0, 0);
            const vector3 a2 = vector3(-bead_radius, bead_radius*sqrt(3.0f), 0);
            const vector3 a3 = vector3(0, 0, (float)(4.0f*sqrt(6)/3.0f*bead_radius)); // seems fine, spreads out

            // basis 1 is (0,0,0)
            const vector3 basis2 = a1*(2.0f/3.0f) + a2*(1.0f/3.0f) + a3/2.0f;
            float xlimit = xaxis*0.5f;
            float ylimit = yaxis*0.5f;
            float zlimit = zaxis*0.5f;

            for (int n3 = -z_prime_index; n3<z_prime_index; n3++){
                for (int n2 = -y_prime_index; n2 <= y_prime_index; n2++) {
                    for (int n1 = -x_prime_index; n1 <= x_prime_index; n1++) {
                        const vector3 temp = (a1 * n1 + a2 * n2 + a3 * n3);
                        if (std::abs(temp.x) < xlimit && std::abs(temp.y) < ylimit && std::abs(temp.z) < zlimit){
                            beads.emplace_back(Bead(temp.x, temp.y, temp.z, 1));
                        }

                        vector3 temp2 = temp + basis2;
                        if (std::abs(temp2.x) < xlimit && std::abs(temp2.y) < ylimit && std::abs(temp2.z) < zlimit) {
                            beads.emplace_back(Bead(temp2.x, temp2.y, temp2.z, 1));
                        }
                    }
                }
            }
        }
    }

    number_of_beads = (unsigned int) beads.size();
    logger("TOTAL LATTICE POINTS",std::to_string(number_of_beads));
}

/*
 * Instead of using a universe with a geometric shape, use one from a mask
 * Use if you make a model at 0.2 and then scale to 0.4 or higher resolution
 */
void PointSetModel::createUniverseFromMask(std::string filename) {
    beads.clear();
    std::cout << " --------------------------------------------------- " << std::endl;

    PDBModel tempmodel (filename, true, false);
    unsigned int totalCoords = tempmodel.getTotalCoordinates();
    std::set<unsigned int> keepers;

    if(!tempmodel.getEdgeRadiusStatus()){
        throw std::invalid_argument("*** ERROR => MASK FILE MUST INCLUDE EDGE RADIUS SPECIFIED AS: REMARK 265                   EDGE RADIUS : 2.5 ");
    }

    float edge_radius = tempmodel.getEdgeRadius() + bead_radius;
    logger("=> CREATING UNIVERSE","MASKED");

    const auto z_prime_index = (int) (2*std::ceil(radial_limit*inv_bead_radius*0.5))+10;
    const auto x_prime_index = (int) (2*std::ceil(radial_limit*inv_bead_radius*0.5))+10;
    const auto y_prime_index = (int) (2*std::ceil(radial_limit*inv_bead_radius*0.5))+10;

    const vector3 a1 = vector3(2*bead_radius, 0, 0);
    const vector3 a2 = vector3(-bead_radius, bead_radius*sqrt(3.0f), 0);
    const vector3 a3 = vector3(0, 0, (float)(4.0f*sqrt(6)/3.0f*bead_radius)); // seems fine, spreads out

    const vector3 basis2 = a1*(2.0/3.0) + a2*(1/3.0) + a3/2.0;

    float part_x = 0.0f, part_y = 0.0f, part_z = 0.0f, part_length = 0.0f;
    for (int n3 = -z_prime_index; n3<z_prime_index; n3++){
        for (int n2 = -y_prime_index; n2 <= y_prime_index; n2++) {
            for (int n1 = -x_prime_index; n1 <= x_prime_index; n1++) {

                const vector3 temp = (a1 * n1 + a2 * n2 + a3 * n3);

                if ((temp.length()) <= radial_limit){
                    for(unsigned int i = 0; i<totalCoords; i++) {
                        vector3 maskVec(*(tempmodel.getCenteredXVec() + i), *(tempmodel.getCenteredYVec() + i),
                                        *(tempmodel.getCenteredZVec() + i));

                        float len = (temp - maskVec).length();
                        if (len <= edge_radius || len <= bead_radius){
                            beads.emplace_back(Bead(temp.x, temp.y, temp.z, 1));
                            part_x += temp.x;
                            part_y += temp.y;
                            part_z += temp.z;
                            part_length += temp.length()*temp.length();
                            break;
                        }
                    }
                }

                vector3 temp2 = temp + basis2;

                if ((temp2.length()) <= radial_limit){
                    for(unsigned int i = 0; i<totalCoords; i++) {
                        vector3 maskVec(*(tempmodel.getCenteredXVec() + i), *(tempmodel.getCenteredYVec() + i),
                                        *(tempmodel.getCenteredZVec() + i));
                        float len = (temp2 - maskVec).length();
                        if (len <= edge_radius || len <= bead_radius){
                            beads.emplace_back(Bead(temp2.x, temp2.y, temp2.z, 1));
                            part_x += temp2.x;
                            part_y += temp2.y;
                            part_z += temp2.z;
                            part_length += temp2.length()*temp2.length();
                            break;
                        }
                    }
                }
            }
        }
    }


    // rg calculations within sphere
    float totalInModel = beads.size();
    std::cout << "center " << (part_x)/totalInModel << " " << (part_y)/totalInModel << " " << (part_z)/totalInModel << std::endl;
    const vector3 part_mean = vector3((part_x)/totalInModel, (part_y)/totalInModel, (part_z)/totalInModel);
    float part_rg = std::sqrt(-part_mean.length()*part_mean.length() + part_length/totalInModel);
    std::cout << "particle :: " << part_mean.length() << " " << part_rg << std::endl;

    unsigned int countIt = 0;
    float comp_x = 0.0f, comp_y = 0.0f, comp_z = 0.0f, comp_length=0.0;
    for (int n3 = -z_prime_index; n3<z_prime_index; n3++){
        for (int n2 = -y_prime_index; n2 <= y_prime_index; n2++) {
            for (int n1 = -x_prime_index; n1 <= x_prime_index; n1++) {

                const vector3 temp = (a1 * n1 + a2 * n2 + a3 * n3);

                if ((temp.length()) <= radial_limit){
                    bool useIt = true;
                    for(unsigned int i = 0; i<totalCoords; i++) {
                        vector3 maskVec(*(tempmodel.getCenteredXVec() + i), *(tempmodel.getCenteredYVec() + i),
                                        *(tempmodel.getCenteredZVec() + i));

                        float len = (temp - maskVec).length();
                        if (len <= edge_radius || len <= bead_radius){
                            useIt = false;
                            break;
                        }
                    }

                    if (useIt){
                        comp_x += temp.x;
                        comp_y += temp.y;
                        comp_z += temp.z;
                        comp_length += temp.length()*temp.length();
                        countIt++;
                    }
                }

                vector3 temp2 = temp + basis2;

                if ((temp2.length()) <= radial_limit){
                    bool useIt = true;
                    for(unsigned int i = 0; i<totalCoords; i++) {
                        vector3 maskVec(*(tempmodel.getCenteredXVec() + i), *(tempmodel.getCenteredYVec() + i),
                                        *(tempmodel.getCenteredZVec() + i));
                        float len = (temp2 - maskVec).length();
                        if (len <= edge_radius || len <= bead_radius){
                            useIt = false;
                            break;
                        }
                    }

                    if (useIt){
                        comp_x += temp2.x;
                        comp_y += temp2.y;
                        comp_z += temp2.z;
                        comp_length += temp2.length()*temp2.length();
                        countIt++;
                    }
                }
            }
        }
    }

    totalInModel = (float)countIt;
    std::cout << "center 2 " << (comp_x)/totalInModel << " " << (comp_y)/totalInModel << " " << (comp_z)/totalInModel << std::endl;
    const vector3 comp_mean = vector3((comp_x)/totalInModel, (comp_y)/totalInModel, (comp_z)/totalInModel);
    float comp_rg = std::sqrt(-comp_mean.length()*comp_mean.length() + comp_length/totalInModel);
    std::cout << "particle :: " << comp_mean.length() << " " << comp_rg << " r_pp " << (comp_mean - part_mean).length() << std::endl;

    number_of_beads = (unsigned int) beads.size();
    logger("TOTAL LATTICE POINTS",std::to_string(number_of_beads));

    float totalOf = number_of_beads + totalInModel;

    float rpp = (comp_mean - part_mean).length();
    float rgsphere = totalInModel/totalOf*comp_rg*comp_rg + number_of_beads/totalOf*part_rg*part_rg + totalInModel/totalOf*number_of_beads/totalOf*rpp*rpp;
    std::cout << " rg sphere " << std::sqrt(rgsphere) << " :: " << std::sqrt((3.0/5.0*radial_limit*radial_limit)) << std::endl;


    exit(0);
}

/*
 * Instead of using a universe with a geometric shape, use one from a mask
 * Use if you make a model at 0.2 and then scale to 0.4 or higher resolution
 */
void PointSetModel::createUniverseFromMaskSym(std::string filename) {
    beads.clear();
    std::cout << " --------------------------------------------------- " << std::endl;

    PDBModel tempmodel (filename, true, false);
    unsigned int totalCoords = tempmodel.getTotalCoordinates();
    std::set<unsigned int> keepers;

    if(!tempmodel.getEdgeRadiusStatus()){
        throw std::invalid_argument("*** MASK FILE MUST INCLUDE EDGE RADIUS SPECIFIED AS: \nREMARK 265                   EDGE RADIUS : 2.5 ");
    }

    float edge_radius = tempmodel.getEdgeRadius() + bead_radius;
    logger("=> CREATING UNIVERSE","MASKED");

    const auto z_prime_index = (int) (2*std::ceil(radial_limit*inv_bead_radius*0.5))+10;
    const auto x_prime_index = (int) (2*std::ceil(radial_limit*inv_bead_radius*0.5))+10;
    const auto y_prime_index = (int) (2*std::ceil(radial_limit*inv_bead_radius*0.5))+10;

    const vector3 a1 = vector3(2*bead_radius, 0, 0);
    const vector3 a2 = vector3(-bead_radius, bead_radius*sqrt(3.0f), 0);
    const vector3 a3 = vector3(0, 0, (float)(4.0f*sqrt(6)/3.0f*bead_radius)); // seems fine, spreads out

    const vector3 basis2 = a1*(2.0/3.0) + a2*(1/3.0) + a3/2.0;

    for (int n3 = -z_prime_index; n3<z_prime_index; n3++){
        for (int n2 = -y_prime_index; n2 <= y_prime_index; n2++) {
            for (int n1 = -x_prime_index; n1 <= x_prime_index; n1++) {

                const vector3 temp = (a1 * n1 + a2 * n2 + a3 * n3);

                if ((temp.length()) <= radial_limit){
                    for(unsigned int i = 0; i<totalCoords; i++) {
                        vector3 maskVec(tempmodel.getX()[i], tempmodel.getY()[i],tempmodel.getZ()[i]);

                        float len = (temp - maskVec).length();
                        if (len <= edge_radius ){
                            beads.emplace_back(Bead(temp.x, temp.y, temp.z, 1));
                            break;
                        }
                    }
                }

                vector3 temp2 = temp + basis2;

                if ((temp2.length()) <= radial_limit){
                    for(unsigned int i = 0; i<totalCoords; i++) {
                        vector3 maskVec(tempmodel.getX()[i], tempmodel.getY()[i],tempmodel.getZ()[i]);
                        float len = (temp2 - maskVec).length();
                        if (len <= edge_radius ){
                            beads.emplace_back(Bead(temp2.x, temp2.y, temp2.z, 1));
                            break;
                        }
                    }
                }
            }
        }
    }

    number_of_beads = (unsigned int) beads.size();
    logger("TOTAL LATTICE POINTS",std::to_string(number_of_beads));
}

/**
 * for each bead in the bead universe, calculate distances and convert to spherical coordinates
 * pre-populate a neighbors list for each bead based on a distance cutoff
 */
void PointSetModel::createDistancesAndConvertToSphericalCoordinates(){

    vector3 diff;

    totalDistances = ( ((unsigned long int)number_of_beads*(number_of_beads-1)) / 2);

    logger("TOTAL DISTANCES",std::to_string(totalDistances));
    logger("MAX MEMORY",std::to_string(bead_indices.max_size()));
    logger("MAX SHORT INDEX",std::to_string(std::numeric_limits<unsigned short int>::max()));
    logger("MAX INDEX",std::to_string(std::numeric_limits<int>::max()));
    logger("MAX UNSIGNED INDEX",std::to_string(std::numeric_limits<unsigned int>::max()));
    logger("MAX UL INDEX",std::to_string(std::numeric_limits<unsigned long int>::max()));

    // total beads must be less than std::numeric_limits<unsigned int>::max()
    bead_indices.resize(number_of_beads);
    std::cout << "*******************         RESIZING NEIGHBORS VECTOR      *******************" << std::endl;
    neighbors.resize(sizeOfNeighborhood*number_of_beads);
    neighborLimit =number_of_beads + 1; // positions with no neighbor will have this limit
    logger("NEIGHBORS SIZE",std::to_string(neighbors.size()));
    startingSetVector.resize(number_of_beads);
    std::fill(neighbors.begin(), neighbors.end(), neighborLimit);

    try{

        if (totalDistances > std::numeric_limits<unsigned int>::max() || (sizeOfNeighborhood*number_of_beads > std::numeric_limits<int>::max())){
//        if (true){

            // distances and bins can't not be precalculated
            std::cout << " SLOW MODE DISTANCES AND BINS CAN NOT BE PRECALCULATED " << std::endl;

            float ratio;
            float invCutOff = 1.0f/cutOffNeighbor;

            std::cout << "*******************           POPULATING NEIGHBORS         *******************" << std::endl;
            for (unsigned int n=0; n < number_of_beads; n++) {

                const vector3 * pCurrentVec = &(beads[n].getVec());
                // populate distance matrix as 1-D
                // very slow, for each bead, calculate distances to every other bead and determine neighbors
                unsigned int count=0;
                for(unsigned int m=0; m < n; m++){
                    diff = (*pCurrentVec) - (&(beads[m]))->getVec();
                    ratio = diff.length()*invCutOff;
                    if (diff.length() <= cutOffNeighbor || (std::abs(1.0f - ratio) < 100*FLT_EPSILON) ){
                        neighbors[sizeOfNeighborhood*n + count] = m;
                        count++;
                    }
                }

                for(unsigned int m=(n+1); m < number_of_beads; m++){
                    diff = (*pCurrentVec) - (&(beads[m]))->getVec();
                    ratio = diff.length()*invCutOff;
                    if (diff.length() <= cutOffNeighbor || (std::abs(1.0f - ratio) < 100*FLT_EPSILON)){
                        neighbors[sizeOfNeighborhood*n + count] = m;
                        count++;
                    }
                }

                if (count > 12){
                    throw std::invalid_argument("NEIGHBOR COUNT TOO HIGH ");
                }

                bead_indices[n] = n;
            }

            std::cout << "*******************            FINISHED NEIGHBORS          *******************" << std::endl;

            //throw std::invalid_argument( "PHYSICAL SYSTEM IS TOO SMALL < NOT ENOUGH MEMEORY => REDUCE RESOLUTION: \n");
            // use directMethohd
            useDirectMethod=true;

        } else {
            std::cout << "*******************                                        *******************" << std::endl;
            std::cout << "*******************        RESIZING DISTANCES VECTOR       *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            //distances.resize(totalDistances);
//            std::vector<float> distances(totalDistances); // linear array (n*(n-1)/2)
//            std::cout << "            DISTANCES SIZE : " << distances.size() << std::endl;
            std::cout << "*******************           RESIZING BINS VECTOR         *******************" << std::endl;
            logger("BINS MAX SIZE",std::to_string(bins.max_size()));
            bins.resize(totalDistances); // this will hold pre-calculated bin assignments
            logger("BINS SIZE",std::to_string(bins.size()));

            std::cout << "*******************           POPULATING NEIGHBORS         *******************" << std::endl;
            /*
             * populate neighbors list
             * distances and neighbors should be made thread safe
             * count should be atomic
             */
            unsigned int * pB = bead_indices.data();

            for (unsigned int n=0; n < number_of_beads; n++) {
                pB[n] = n;
                const vector3 * pCurrentVec = &(beads[n].getVec());
                // populate distance matrix as 1-D
                // very slow, for each bead, calculate distances to every other bead and determine neighbors
                unsigned int count=0;
                for(unsigned int m=0; m < n; m++){
                    diff = (*pCurrentVec) - (&(beads[m]))->getVec();
                    if (diff.length() < cutOffNeighbor){
                        neighbors[sizeOfNeighborhood*n + count] = m;
                        count++;
                    }
                }

                for(unsigned int m=(n+1); m < number_of_beads; m++){
                    diff = (*pCurrentVec) - (&(beads[m]))->getVec();
                    if (diff.length() < cutOffNeighbor){
                        neighbors[sizeOfNeighborhood*n + count] = m;
                        count++;
                    }
                }
            }

            std::cout << "*******************            FINISHED NEIGHBORS          *******************" << std::endl;
        }

    } catch (std::exception &err) {
        std::cerr<<"Caught "<<err.what()<< std::endl;
        std::cerr<<"Type "<<typeid(err).name()<< std::endl;
        exit(0);
    }
}


/**
 * indices of lattice points that comprise neighborhood of lattice point at index
 */
std::vector<unsigned int>::iterator PointSetModel::getPointerToNeighborhood(unsigned int index){
    return (this->neighbors.begin() + this->sizeOfNeighborhood*index);
}


/**
 * Go over all the lattice points in the Universe and convert each distance pair to a bin
 * Also determine max bin.
 * There will be n*(n-1)/2 distances or bins
 *
 * @param pData
 * @return
 */
unsigned short int PointSetModel::populateBins(PofRData * pData) {

    logger("CALCULATING", "POPULATING BINS");
    //unsigned short int * pBin = &bins[0];
    const Bead * pBead = beads.data();
    unsigned short int * pBin = bins.data();
    unsigned short int maxbin = 0;

    //Bead * currentbead;
    unsigned long int discount=0;
    for (unsigned int n=0; n < number_of_beads; n++) {

        const vector3 & pVec1 = (pBead+n)->getVec();
//        const vector3 & pVec1 = beads[n].getVec();

        unsigned int next = n+1;
        // populate distance matrix as 1-D
        for(unsigned int m=next; m < number_of_beads; m++){
//            *(pBin+discount) = pData->convertToBin((pVec1 - (&(beads[m]))->getVec()).length());
            *(pBin+discount) = pData->convertToBin((pVec1 - ((pBead+m))->getVec()).length());
            if (*(pBin+discount) > maxbin){
                maxbin = *(pBin+discount);
            }
            discount++;
        }

        if (pBin != bins.data()){
            throw std::invalid_argument(" CATASTROPHIC POINTER ERROR from PointSetModel::populateBins ");
            exit(0);
        }
    }

//    std::vector<unsigned int> values(maxbin);
//
//    unsigned int totalcounts = 0;
//    for(unsigned int i=0; i<(maxbin+1); i++){ // distribution of contacts
//
//        unsigned int count = 0;
//        for(auto bin : bins){
//            if (bin == i){
//                ++count;
//            }
//        }
//        totalcounts += count;
//        values[i] = count;
//        std::cout << " " << i << " " << count << std::endl;
//    }
//
//    std::cout << " total " << totalcounts << std::endl;
//    float invCounts = 1.0/(float)totalcounts;
//    for(unsigned int i=0; i<(maxbin+1); i++){ // distribution of contacts
//
//        float value = values[i]*invCounts;
//        float rvalue = pData->getBinWidth()*i + 0.5*pData->getBinWidth();
//        float calc = 3.0*std::pow(rvalue, 2)/std::pow(radial_limit,3) - 9.0/4.0*std::pow(rvalue, 3)/std::pow(radial_limit,4) + 3.0/16.0*std::pow(rvalue, 5)/std::pow(radial_limit,6);
//        std::cout << " " << rvalue << " " << value << " " << calc  << std::endl;
//    }
//    std::cout << " CONTACTS " << contacts << std::endl;
    logger("MAX BIN", std::to_string(maxbin));
    return maxbin;
}

void PointSetModel::pruneUniverseToAsymmetricUnit(){
    /*
     * set angular range for
     */
    //const vector3 normX = vector3(1,0,0);
    double angle = M_PI*2.0d/(double)symmetryIndex;
    double cosLimit = std::cos(angle*0.5d);

    logger("ROT ANGLE", SASTOOLS_UTILS_H::formatNumber(angle*180/M_PI,2)); // half angle centered on X-axis
    logger("COSINE LIMIT", SASTOOLS_UTILS_H::formatNumber(cosLimit, 3));

    std::vector<vector3> tempbeads;

    if (symmetryGroup=="D"){

        for(auto & bead : beads){
            // calculate angle between bead and XY-plane
            if (bead.getX() >= 0 && bead.getZ() >= 0){

                const vector3 * pVec = &bead.getVec();
                vector3 tVec = vector3(pVec->x, pVec->y, 0);
                double test = pVec->x/tVec.length();

                if (test >= cosLimit){
                    tempbeads.emplace_back(vector3(*pVec));
                }

            }
        }
    } else { // rotational symmetry such as C
        for(auto & bead : beads){
            // calculate angle between bead and XY-plane
            if (bead.getX() >= 0){
                const vector3 * pVec = &bead.getVec();
                vector3 tVec = vector3(pVec->x, pVec->y, 0);
                double test = pVec->x/tVec.length();

                if (test >= cosLimit){
                    tempbeads.emplace_back(vector3(*pVec));
                }
            }
        }
    }

    logger("TOTAL ASYMMETRIC UNIT ",std::to_string(tempbeads.size()));
    beads.clear();
    beads.reserve(tempbeads.size());
    for(auto & vec : tempbeads){
        beads.emplace_back(Bead(vec.x, vec.y, vec.z, 1));
    }

    number_of_beads = (unsigned int) beads.size();
    logger("TOTAL LATTICE POINTS",std::to_string(number_of_beads));
}

void PointSetModel::writeModelToFile(const unsigned int workingLimit, std::vector<unsigned int> &selectedBeads, std::string & nameOf, unsigned int steps){
    //const char *outputFileName;
    nameOf = nameOf+"_" + std::to_string(steps) + ".pdb";
    FILE * pFile = fopen(nameOf.c_str(), "w");

    Bead * currentBead;
    unsigned int residue_index;

    fprintf(pFile, "REMARK 265                   EDGE RADIUS : %.3f\n", this->bead_radius);
    for (unsigned int i=0; i < workingLimit; i++){
        currentBead = this->getBead(selectedBeads[i]);
        residue_index = (i + 1);
        SASTOOLS_UTILS_H::printAtomLine(pFile, i+1, "A", residue_index, currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    }

    fclose(pFile);
}

void PointSetModel::writeModelToFileFlipped(const unsigned int workingLimit, std::vector<unsigned int> &selectedBeads, std::string & nameOf, unsigned int steps){
    //const char *outputFileName;
    nameOf = nameOf+"_" + std::to_string(steps) + ".pdb";
    FILE * pFile = fopen(nameOf.c_str(), "w");

    Bead * currentBead;
    unsigned int residue_index;

    fprintf(pFile, "REMARK 265                   EDGE RADIUS : %.3f\n", this->bead_radius);
    for (unsigned int i=0; i < workingLimit; i++){
        currentBead = this->getBead(selectedBeads[i]);
        residue_index = (i + 1);
        SASTOOLS_UTILS_H::printAtomLine(pFile, i+1, "A", residue_index, -currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    }

    fclose(pFile);
}


std::string PointSetModel::writeBasicSymModelToFile(float dkl, unsigned int workingLimitS, std::vector<unsigned int> &beads, std::string name) {

    char alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" ;
    std::vector<char> alphabet( alpha, alpha+sizeof(alpha)-1 ) ;

    //FILE * pFile;

    //const char *outputFileName;
    std::string newName = name + ".pdb";
    const char * outputFileName = newName.c_str() ;
    FILE * pFile = fopen(outputFileName, "w");

//    std::string temp = createHeader(dkl, annealedObject, pData, totalSteps, workingLimitS, volume, averageContacts);
//    fprintf(pFile, temp.c_str());

    unsigned int residue_index;

    unsigned int totalCoordinates = numberOfSubUnits*workingLimitS;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;
    Bead * tempBead;
    // create sym partners and add to Pr
    // calculate P(r) for subunit
    for (unsigned int i=0; i<workingLimitS; i++){
        tempBead = this->getBead(beads[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    unsigned int count = workingLimitS;
    for (unsigned int s=1; s<numberOfSubUnits; s++){  // create sym related subunits and add to coordinates vector
        this->transformCoordinatesBySymmetry(s, workingLimitS, count, coordinates);
    }

    std::string chain;
    unsigned int index=0;
    for (unsigned int s=0; s<numberOfSubUnits; s++){  // create sym related subunits and add to coordinates vector

        chain = alphabet[s];

        for (unsigned int i=0; i<workingLimitS; i++){
            // convert coordinate to
            residue_index = i+1;
            //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %8.3f%8.3f%8.3f  1.00 100.00\n", "ATOM", i+1, "CA", "ALA", chain.c_str(), residue_index.c_str(), coordinates[index].x, coordinates[index].y, coordinates[index].z );
            printAtomLine(pFile, i+1, chain, residue_index, coordinates[index].x, coordinates[index].y, coordinates[index].z );
            index++;
        }
        fprintf(pFile, "TER\n");
    }
    fprintf(pFile, "END\n");
    fclose(pFile);
    return newName;
}


vector3 PointSetModel::transformVectorBySymmetry(unsigned int subunitIndex, const vector3 & vec){
    auto angle = (float)(M_PI*2.0f*subunitIndex/(float)symmetryIndex);
    float cosine = cos(angle);
    float sine = sin(angle);
    float x, y, z;

    if (subunitIndex < symmetryIndex){ // rotations for C and D symmetry groups
        x = (vec).x;
        y = (vec).y;
        z = (vec).z;
    } else { // reflections for the D_group
        // rotate on x-axis first by 180 and then by
        // rotating on x-axis (x,y,z) -> (x,-y,-z)
        x = vec.x;
        y = -vec.y;
        z = -vec.z;
    }

    return {cosine*x - sine*y, sine*x + cosine*y, z};
}

/**
 * transform subunit contained within coordinates limited to workingLimit
 */
void PointSetModel::transformCoordinatesBySymmetry(const unsigned int subunitIndex, const unsigned int workingLimit, unsigned int &startCount, std::vector<vector3> &coordinates){
    vector3 * transformed, * nextVec;
    // if subunit index > sym_index, implies D_n, otherwise, its all rotations

    auto angle = (float)(M_PI*2.0f*subunitIndex/(float)symmetryIndex);
    float cosine = cos(angle);
    float sine = sin(angle);
    float x, y, z;

    if (subunitIndex < symmetryIndex){ // rotations for C and D symmetry groups

        for(unsigned int i=0; i<workingLimit; i++){
            transformed = &coordinates[i];
            x = (*transformed).x;
            y = (*transformed).y;
            z = (*transformed).z;

            nextVec = &coordinates[startCount];
            (*nextVec).x = cosine*x - sine*y;
            (*nextVec).y = sine*x + cosine*y;
            (*nextVec).z = z;

            //cout << startCount << " " << angle << " SUBUNIT " << subunitIndex << " " << (*nextVec).x << " " << (*nextVec).y << " " << (*nextVec).z << endl;
            startCount++;
        }

    } else { // reflections for the D_group
        // rotate on x-axis first by 180 and then by
        // rotating on x-axis (x,y,z) -> (x,-y,-z)
        for(unsigned int i=0; i<workingLimit; i++){
            transformed = &coordinates[i];
            x = (*transformed).x;
            y = -(*transformed).y;
            z = -(*transformed).z;

            nextVec = &coordinates[startCount];
//            (*nextVec).x = cosine*x + sine*y;
//            (*nextVec).y = sine*x - cosine*y;
            (*nextVec).x = cosine*x - sine*y;
            (*nextVec).y = sine*x + cosine*y;
            (*nextVec).z = z;
            //cout << startCount << " " << angle << " SUBUNIT " << subunitIndex << " " << (*nextVec).x << " " << (*nextVec).y << " " << (*nextVec).z << endl;
            startCount++;
        }
    }
}


/*
 * find the largest distance within the search space
 * for helical symmetry, search space is largest within cylinder
 *
 */
unsigned short int PointSetModel::getMaxBin(PofRData * pData) {

    unsigned short int maxbin = 0, tempBin;

    Bead * currentbead;
    for (unsigned int n=0; n < number_of_beads; n++) {

        currentbead = &(beads[n]);
        const vector3 & tempVec = currentbead->getVec();

        unsigned int next = n+1;
        // populate distance matrix as 1-D
        for(unsigned int m=next; m < number_of_beads; m++){
            tempBin = pData->convertToBin((tempVec - (&(beads[m]))->getVec()).length());
            if (tempBin > maxbin){
                maxbin = tempBin;
            }
        }
    }

    return maxbin;
}

/**
 * write subset of selected bead to file
 */
void PointSetModel::writeSubModelToFile(unsigned int startIndex, unsigned int workingLimit, std::vector<unsigned int> &selectedBeads, std::string nameOf){

    nameOf = nameOf + ".pdb";
    //const char * outputFileName = nameOf.c_str();
    FILE * pFile;
    pFile = fopen(nameOf.c_str(), "w");

    unsigned int residue_index;
    std::vector<std::string> chains(25);
    chains[0] = "A";
    chains[1] = "B";
    chains[2] = "C";
    chains[3] = "D";
    chains[4] = "E";
    chains[5] = "F";
    chains[6] = "G";
    chains[7] = "H";
    chains[8] = "I";
    chains[9] = "K";
    chains[10] = "L";
    chains[11] = "M";
    chains[12] = "N";
    chains[13] = "O";
    chains[14] = "P";
    chains[15] = "Q";
    chains[16] = "R";
    chains[17] = "S";
    chains[18] = "T";
    chains[19] = "U";
    chains[20] = "V";
    chains[21] = "W";
    chains[22] = "X";
    chains[23] = "Y";
    chains[24] = "J";

    std::string chainID = chains[0];
    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    unsigned int residueCnt = 1;

    for (unsigned int i=startIndex; i<workingLimit; i++){
        Bead * currentBead = this->getBead(selectedBeads[i]);
        residue_index = residueCnt;

        if (i%9999 == 0 && i > startIndex){ // reinitialize counter
            residueCnt = 1;

            auto locale = (unsigned int)std::floor((double)i/9999.0);
            locale = (locale >= 25) ? 25-locale : locale;
            if (locale > 0){
                chainID = chains[locale];
            }
        }

        printAtomLine(pFile, residueCnt, chainID, residue_index, currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        residueCnt++;
    }

    fclose(pFile);
}


std::string PointSetModel::writeSymModelToFile(float dkl, unsigned int workingLimitS, std::vector<unsigned int> &beads, std::vector<unsigned int> &pofrModel, std::string name, Anneal *annealedObject, PofRData *pData, unsigned int totalSteps, float volume, float averageContacts) {

    char alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" ;
    std::vector<char> alphabet( alpha, alpha+sizeof(alpha)-1 ) ;

    //FILE * pFile;

    //const char *outputFileName;
    std::string newName = name + ".pdb";
    const char * outputFileName = newName.c_str() ;
    FILE * pFile = fopen(outputFileName, "w");

    std::string temp = createHeader(dkl, annealedObject, pData, totalSteps, workingLimitS, volume, averageContacts);
    fprintf(pFile, temp.c_str());

    unsigned int residue_index;

    unsigned int totalCoordinates = numberOfSubUnits*workingLimitS;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;
    Bead * tempBead;
    // create sym partners and add to Pr
    // calculate P(r) for subunit
    for (unsigned int i=0; i<workingLimitS; i++){
        tempBead = this->getBead(beads[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    unsigned int count = workingLimitS;
    for (unsigned int s=1; s<numberOfSubUnits; s++){  // create sym related subunits and add to coordinates vector
        this->transformCoordinatesBySymmetry(s, workingLimitS, count, coordinates);
    }

    std::string chain;
    unsigned int index=0;
    for (unsigned int s=0; s<numberOfSubUnits; s++){  // create sym related subunits and add to coordinates vector

        chain = alphabet[s];

        for (unsigned int i=0; i<workingLimitS; i++){
            // convert coordinate to
            residue_index = i+1;
            //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %8.3f%8.3f%8.3f  1.00 100.00\n", "ATOM", i+1, "CA", "ALA", chain.c_str(), residue_index.c_str(), coordinates[index].x, coordinates[index].y, coordinates[index].z );
            printAtomLine(pFile, i+1, chain, residue_index, coordinates[index].x, coordinates[index].y, coordinates[index].z );
            index++;
        }
        fprintf(pFile, "TER\n");
    }
    fprintf(pFile, "END\n");
    fclose(pFile);
    return newName;
}

void PointSetModel::writeSetToFile(std::set<unsigned int> &selectedBeads, std::string & nameOf){
    //const char *outputFileName;
    nameOf = nameOf+"_set.pdb";
    FILE * pFile = fopen(nameOf.c_str(), "w");

    Bead * currentBead;

    std::vector<std::string> chains(13);
    chains[0] = "A";
    chains[1] = "B";
    chains[2] = "C";
    chains[3] = "D";
    chains[4] = "E";
    chains[5] = "F";
    chains[6] = "G";
    chains[7] = "H";
    chains[8] = "I";
    chains[9] = "K";
    chains[10] = "L";
    chains[11] = "M";
    chains[12] = "N";

    std::string chainID = chains[0];

    unsigned int count = 1;
    unsigned int residue_index = 1;
    for(auto ind : selectedBeads){
        currentBead = this->getBead(ind);


        if (residue_index%9999 == 0){ // reinitialize counter
            residue_index = 1;

            auto locale = (unsigned int)std::floor((double)count/9999.0);

            if (locale > 0){
                chainID = chains[locale];
            }
        }

        printAtomLine(pFile, count, chainID, residue_index, currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        residue_index++;
        count++;
    }

    fclose(pFile);
}


/**
 * center the input set of indices, reassigning the bead indices
 *
 */
void PointSetModel::centerLatticeModel(unsigned int *workingLimit, std::vector<unsigned int> & indices, std::set<unsigned int> & hullPts){

    Bead * pBead;// = &beads[indices[0]];

    const unsigned int activeLimit = *workingLimit;

    // find center of a filled object
    float avex=0, avey=0, avez=0;
    for(auto it : hullPts){
        pBead = &beads[it];
        avex += pBead->getVec().x;
        avey += pBead->getVec().y;
        avez += pBead->getVec().z;
    }

    float invTotal = 1.0f/(float)hullPts.size();

    avex *= invTotal;
    avey *= invTotal;
    avez *= invTotal;
    const vector3 totalVec = vector3(avex, avey, avez);

    // which bead is closest to centering position

    // which bead (when translated) in universe is closest to center?
    float closest=(beads[0].getVec() - totalVec).length();
    unsigned int closestIndex = 0;
    std::vector<unsigned int> indicesNotFound;
    std::vector<unsigned int> indicesNotNotFound;

    //
    // geometric center, not mass center, which bead is closest to average position and also on the basis plane
    //
    // less than or equal to
    // lower < distance <= upper
    // bin = 0 is the distance between two beads
    // due to floating point errors, sometimes distance will be slightly larger than 1
    // need to subtract 1 to make it zero
    //

    for(unsigned int i=1; i< number_of_beads; i++) {
        const vector3 *pVec = &beads[i].getVec();
        float ratiox = pVec->x/bead_radius;
        float flooredX = floor(ratiox);
        float diffx = std::abs(ratiox-flooredX);
        float ratioy = pVec->y/(bead_radius*sqrt(3.0f));
        float flooredY = floor(ratioy);
        float diffy = std::abs(ratioy-flooredY);
        float ratioz = pVec->z/(4.0f*sqrt(6)/3.0f*bead_radius);
        float flooredZ = floor(ratioz);
        float diffz = std::abs(ratioz-flooredZ);

        if (diffx <= 100*FLT_EPSILON && diffy <= 100*FLT_EPSILON && diffz <= 100*FLT_EPSILON){
            const vector3 newVec = *pVec - totalVec; // translate the bead
            if (newVec.length() < closest) {
                closest = newVec.length();
                closestIndex = i;
            }
        }
    }


    const vector3 * centeringVec = &beads[closestIndex].getVec();

    std::vector<unsigned int> temp_bead_indices(number_of_beads); // make a copy of indices to re-arrange
    unsigned int searchLimit = 0;
    for (unsigned int j=searchLimit; j<number_of_beads; j++){
        temp_bead_indices[j] = j;
    }

    float distance;

    for(unsigned int i=0; i < activeLimit; i++){
        const vector3 newVec = beads[(indices[i])].getVec() - *centeringVec; // translate the bead
        bool foundIt = false;
        float mindistance = 0.51*1.7312f*bead_radius; // beads should be within each radii
        //float mindistance = 2.0017312f*bead_radius; // beads should be within each radii
        // determine which bead is closest to this new vector

        unsigned int swapindex=0;
        for (unsigned int j=searchLimit; j < number_of_beads; j++){

            distance = (newVec - beads[temp_bead_indices[j]].getVec()).length();

            if (distance < mindistance){ // find the bead that is closest to translated vector
                mindistance = distance;
                swapindex=j;
                foundIt = true;
            }
        }

        if (foundIt){
            std::iter_swap( temp_bead_indices.begin() + searchLimit,  temp_bead_indices.begin() + swapindex);
            searchLimit++;
        } else { // if not found, which is closest?
            // if not found find the neighborhood
            indicesNotFound.push_back(indices[i]);
        }
    }


    if (!indicesNotFound.empty()){ // place any remaining beads
        // check if any bead is within two*beadradius and it contacts a previously selected bead
        for(auto it = indicesNotFound.begin(); it != indicesNotFound.end(); ++it){

            const vector3 newVec = beads[*it].getVec() - *centeringVec; // translate the bead
            bool foundIt = false;
            float mindistance = 4.0f*bead_radius; // beads should be within each radii
            // determine which bead is closest to this new vector
            unsigned int swapindex=0;
            for (unsigned int j=searchLimit; j < number_of_beads; j++){

                const vector3 * pVec = &beads[temp_bead_indices[j]].getVec();
                distance = (newVec - *pVec).length();

                if (distance < mindistance){ // find the bead that is closest to translated vector
                    // is this bead in contact with previously selected bead?
                    for(unsigned int k=0; k<searchLimit; k++) {
                        float contactdistance = (*pVec - beads[temp_bead_indices[k]].getVec()).length();
                        if (contactdistance < cutOffNeighbor) {
                            mindistance = distance;
                            swapindex = j;
                            foundIt = true;
                        }
                    }
                }
            }

            if (foundIt){ // make thhe swap
                std::iter_swap( temp_bead_indices.begin() + searchLimit,  temp_bead_indices.begin() + swapindex);
                searchLimit++;
            } else {
                indicesNotNotFound.push_back(*it);
            }
        }
    }

    std::copy(temp_bead_indices.begin(), temp_bead_indices.end(), indices.begin());
    std::sort(indices.begin(), indices.begin() + searchLimit);

    std::cout << "Centering vec Bead Index "<< closestIndex << std::endl;
    beads[closestIndex].printCoordinates();

    if (searchLimit != *workingLimit){
        std::cout << "LIMITS DO NOT MATCH" << std::endl;
//        exit(0);
    }

//    std::string test = "missing";
//    unsigned int totalnotfound = indicesNotFound.size();
//    this->writeModelToFile(totalnotfound, indicesNotFound, test, 0);
    //test = "notnot";
    //this->writeModelToFile(indicesNotNotFound.size(), indicesNotNotFound, test, 0);
    *workingLimit = searchLimit;
    //this->writeModelToFile(*workingLimit, temp_bead_indices, "centered_model", 0);
}




std::string PointSetModel::writeModelToFile2(float dkl, unsigned int workingNumber, std::vector<unsigned int> &selectedBeads, std::vector<unsigned int> &pofrModel, std::string nameOf, Anneal *annealedObject, PofRData *pData, unsigned int steps, float volume, float averageContacts){
    FILE * pFile;

    const char *outputFileName;
    nameOf = nameOf + ".pdb";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    unsigned int residue_index;
    std::string temp = createHeader(dkl, annealedObject, pData, steps, workingNumber, volume, averageContacts);
    fprintf(pFile, temp.c_str());

    // Add P(r) distributions
    // final D_kl, volume and energy
    float totalCounts = 0.0;
    float prob;
    int totalm = pofrModel.size();

    // normalization constant of model Pr
    // treats each value as discrete (i.e. not integrating via trapezoid)
    for (int i=0; i<totalm; i++){
        totalCounts += pofrModel[i];
    }

    unsigned int shannon_bins = annealedObject->gettotalBinsDerivedFromData();
    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
    float binwidth = pData->getBinWidth(), r_value;
    for (unsigned int i=0; i < shannon_bins; i++){
        prob = pData->getProbabilityPerBin(i);  // bounded by experimental Shannon Number
        r_value = i*binwidth + 0.5f*binwidth;
        fprintf(pFile, "REMARK 265  %8.3f  %.4E  %.4E  %-.4E\n", r_value, prob, (pofrModel[i]/totalCounts), (prob * log(prob/pofrModel[i]*totalCounts)));
    }

//    fprintf(pFile, annealedObject->getConnectivityTableText().c_str());

    fprintf(pFile, "REMARK 265\n");
    // write coordinates
    for (unsigned int i=0; i<workingNumber; i++){
        Bead * currentBead = this->getBead(selectedBeads[i]);
        //residue_index = std::to_string(selectedBeads[i]);
        residue_index = i + 1;
        //residue_index = std::to_string(selectedBeads[i]);
        //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        //fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", i+1," CA ", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        printAtomLine(pFile, i+1, "A", residue_index, currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    }

    fprintf(pFile, "END\n");
    fclose(pFile);
    return nameOf;
}


std::string PointSetModel::createHeader(float dkl, Anneal * annealedObject, PofRData *pData, unsigned int totalSteps, unsigned int workingNumber, float volume, float averageContacts){

    std::string datafilename = pData->getFilename();
    char buffer[80];
    // int cstring;

    std::string tempHeader = "REMARK 265\n";
    tempHeader += "REMARK 265 EXPERIMENTAL DETAILS\n";
    tempHeader += "REMARK 265\n";

    tempHeader += "REMARK 265 EXPERIMENT TYPE : X-RAY SOLUTION SCATTERING\n";
    tempHeader += "REMARK 265 DATA ACQUISITION\n";
    tempHeader += "REMARK 265              RADIATION SOURCE : BENDING MAGNET\n";
    tempHeader += "REMARK 265             SYNCHROTRON (Y/N) : Y\n";
    tempHeader += "REMARK 265                      BEAMLINE : B21 DIAMOND LIGHT SOURCE (example)\n";
    tempHeader += "REMARK 265          TEMPERATURE (KELVIN) : 298 (example)\n";
    tempHeader += "REMARK 265                            PH : 6.9 (example)\n";
    tempHeader += "REMARK 265              BUFFER COMPONENT : 20 mM HEPES (example)\n";
    tempHeader += "REMARK 265              BUFFER COMPONENT : 150 mM KCl (example)\n";
    tempHeader += "REMARK 265              BUFFER COMPONENT : 2 mM TCEP (example)\n";
    tempHeader += "REMARK 265              BUFFER COMPONENT : 5 mM KNitrate (example)\n";
    tempHeader += "REMARK 265\n";
    tempHeader += "REMARK 265       DATA REDUCTION SOFTWARE : SCATTER (v3.x)\n";
    tempHeader += "REMARK 265               SOFTWARE AUTHOR : RP RAMBO\n";
    tempHeader += "REMARK 265        DATA ANALYSIS SOFTWARE : SCATTER (v3.x)\n";
    tempHeader += "REMARK 265               SOFTWARE AUTHOR : RP RAMBO\n";
    tempHeader += "REMARK 265\n";
    tempHeader += "REMARK 265 MODELING METHOD\n";
    tempHeader += "REMARK 265                      SOFTWARE : IKETAMA (eKay-TaMa) v 1.9\n";
    tempHeader += "REMARK 265                        METHOD : Cross-Entropy Minimization\n";
    tempHeader += "REMARK 265                        TARGET : REAL-SPACE P(r)-distribution\n";
    tempHeader += "REMARK 265                        AUTHOR : RP RAMBO\n";
    tempHeader += "REMARK 265 \n";

    std::snprintf(buffer, 80, "REMARK 265  EXPERIMENTAL REAL SPACE FILE : %s\n", datafilename.c_str());
    tempHeader.append(buffer);

    tempHeader += "REMARK 265\n";
    tempHeader += "REMARK 265 SEARCH AND REFINEMENT CONSTRAINTS\n";

    std::snprintf(buffer, 80, "REMARK 265                   EDGE RADIUS : %.3f\n", this->bead_radius);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265      AVG PTS HIGH TEMP SEARCH : %.1f\n", this->beadAverage);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265                         STDEV : %.1f\n", this->beadStDev);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265             SEARCH SPACE DMAX : %.1f Angstrom\n", 2*radial_limit);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265        EXPERIMENTAL P(r) DMAX : %.1f Angstrom\n", pData->getDmax());
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265        EMPIRICAL POROD VOLUME : %i Angstrom^3\n", (int)pData->getVolume());
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265   SEARCH SPACE SHANNON NUMBER : %i \n", annealedObject->getMaxBin());
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265   EXPERIMENTAL SHANNON NUMBER : %i\n", annealedObject->gettotalBinsDerivedFromData());
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265                      SYMMETRY : %s\n", symmetry.c_str());
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265             TEMP RANGE (HIGH) : %.4E \n", annealedObject->getHighTempStartForCooling());
    tempHeader.append(buffer);
//    cstring = snprintf(buffer, 80, "REMARK 265              TEMP RANGE (LOW) : %.4E\n", annealedObject->getLowTempStop());
//    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265       TOTAL TEMPERATURE STEPS : %i\n", totalSteps);
    tempHeader.append(buffer);
//    cstring = snprintf(buffer, 80, "REMARK 265        TEMP STEP SCALE FACTOR : %.4f\n", annealedObject->getStepFactor());
//    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265           ALPHA (COMPACTNESS) : %.4E\n", annealedObject->getAlpha());
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265         LAMBDA (CONNECTIVITY) : %.4E\n", annealedObject->getLambda());
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265              MU (COMPACTNESS) : %.4E\n", annealedObject->getMu());
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265     BETA (LATTICE VIOLATIONS) : %.4E\n", annealedObject->getBeta());
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265                 ETA (CONTACT) : %.4E\n", annealedObject->getEta());
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265            PERCENT ADD REMOVE : %.2f\n", annealedObject->getPercentAddRemove());
    tempHeader.append(buffer);

    tempHeader += "REMARK 265\n";

    tempHeader += "REMARK 265 REFINED VALUES\n";
    std::snprintf(buffer, 80, "REMARK 265 TOTAL LATTICE POINTS IN MODEL : %i \n", workingNumber);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265          TOTAL LATTICE VOLUME : %.0f Angstrom^3\n", workingNumber*bead_volume);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265               CVX HULL VOLUME : %i Angstrom^3\n", (int)volume);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265     AVERAGE CONTACTS PER BEAD : %.2f \n", averageContacts);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265   KULLBACK LEIBLER DIVERGENCE : %.5E\n", dkl);
    tempHeader.append(buffer);

    tempHeader += "REMARK 265\nREMARK 265\tFollowing table describes model-data agreement\nREMARK 265\t<r> is the average r-value per Shannon bin\n";
    tempHeader += "REMARK 265\nREMARK 265      <r>\t   P(r)_OBS   P(r)_MODEL   CROSS-ENTROPY\n";

//cout << tempHeader << endl;
    return tempHeader;
}


void PointSetModel::transformCoordinatesBySymmetryPreCalc(const unsigned int subunitIndex, const unsigned int workingLimit, unsigned int &startCount, std::vector<vector3> &coordinates){

    vector3 * transformed, * nextVec;
    // if subunit index > sym_index, implies D_n, otherwise, its all rotations
    float * cosine = &subUnitAnglesCos[subunitIndex];
    float * sine = &subUnitAnglesSin[subunitIndex];
    // float x, y, z;

    if (subunitIndex < symmetryIndex){ // rotations for C and D symmetry groups

        for(unsigned int i=0; i<workingLimit; i++){
            transformed = &coordinates[i];
//            x = (*transformed).x;
//            y = (*transformed).y;
//            z = (*transformed).z;

            nextVec = &coordinates[startCount];
            (*nextVec).x = *cosine*(*transformed).x - *sine*(*transformed).y;
            (*nextVec).y = *sine*(*transformed).x + *cosine*(*transformed).y;
            (*nextVec).z = (*transformed).z;

            //cout << startCount << " " << angle << " SUBUNIT " << subunitIndex << " " << (*nextVec).x << " " << (*nextVec).y << " " << (*nextVec).z << endl;
            startCount++;
        }

    } else { // reflections for the D_group
        // rotate on x-axis first by 180 and then by
        // rotating on x-axis (x,y,z) -> (x,-y,-z)
        for(unsigned int i=0; i<workingLimit; i++){
            transformed = &coordinates[i];
//            x = (*transformed).x;
//            y = -(*transformed).y;
//            z = -(*transformed).z;

            nextVec = &coordinates[startCount];
//            (*nextVec).x = cosine*x + sine*y;
//            (*nextVec).y = sine*x - cosine*y;
//            (*nextVec).x = *cosine*(*transformed).x - *sine*(-(*transformed).y);
//            (*nextVec).y = *sine*(*transformed).x + *cosine*(-(*transformed).y);
//            (*nextVec).z = -(*transformed).z;
            (*nextVec).x = *cosine*(*transformed).x - (*sine*(*transformed).y);
            (*nextVec).y = -(*sine*(*transformed).x) - (*cosine*(*transformed).y);
            (*nextVec).z = -((*transformed).z);

            //cout << startCount << " " << angle << " SUBUNIT " << subunitIndex << " " << (*nextVec).x << " " << (*nextVec).y << " " << (*nextVec).z << endl;
            startCount++;
        }
    }


    if (false){ //icosahedral

        if (subunitIndex == 2){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = 0.30902f*(*transformed).x - 0.80902f*(*transformed).y + 0.5f*(*transformed).z;
                (*nextVec).y = 0.80902f*(*transformed).x - 0.5f*(*transformed).y + 0.30902f*(*transformed).z;
                (*nextVec).z = -0.5f*(*transformed).x + 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 3){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.80902f*(*transformed).x - 0.5f*(*transformed).y + 0.30902f*(*transformed).z;
                (*nextVec).y = 0.5f*(*transformed).x - 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;
                (*nextVec).z = -0.30902f*(*transformed).x + 0.80902f*(*transformed).y + 0.5f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 4){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.80902f*(*transformed).x + 0.5f*(*transformed).y - 0.30902f*(*transformed).z;
                (*nextVec).y = -0.5f*(*transformed).x - 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;
                (*nextVec).z =  0.30902f*(*transformed).x + 0.80902f*(*transformed).y + 0.5f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 5){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = 0.30902f*(*transformed).x + 0.80902f*(*transformed).y - 0.5f*(*transformed).z;
                (*nextVec).y = -0.80902f*(*transformed).x + 0.5f*(*transformed).y + 0.30902f*(*transformed).z;
                (*nextVec).z =  0.5f*(*transformed).x + 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 6){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -1.0f*(*transformed).x;
                (*nextVec).y = -1.0f*(*transformed).y;
                (*nextVec).z =  1.0f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 7){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.30902f*(*transformed).x - 0.80902f*(*transformed).y + 0.5f*(*transformed).z;
                (*nextVec).y =  0.80902f*(*transformed).x - 0.5f*(*transformed).y - 0.30902f*(*transformed).z;
                (*nextVec).z =  0.5f*(*transformed).x + 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 8){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = 0.30902f*(*transformed).x + 0.80902f*(*transformed).y + 0.5f*(*transformed).z;
                (*nextVec).y = -0.80902f*(*transformed).x + 0.5f*(*transformed).y - 0.30902f*(*transformed).z;
                (*nextVec).z =  -0.5f*(*transformed).x - 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 9){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.30902f*(*transformed).x + 0.80902f*(*transformed).y + 0.5f*(*transformed).z;
                (*nextVec).y = -0.80902f*(*transformed).x - 0.5f*(*transformed).y + 0.30902f*(*transformed).z;
                (*nextVec).z =  0.5f*(*transformed).x - 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 10){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.5000f*(*transformed).x + 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;
                (*nextVec).y =  0.30902f*(*transformed).x - 0.80902f*(*transformed).y + 0.5000f*(*transformed).z;
                (*nextVec).z =  0.80902f*(*transformed).x + 0.5000f*(*transformed).y + 0.30902f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 11){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.5000f*(*transformed).x - 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;
                (*nextVec).y = -0.30902f*(*transformed).x + 0.80902f*(*transformed).y + 0.5000f*(*transformed).z;
                (*nextVec).z = -0.80902f*(*transformed).x - 0.5000f*(*transformed).y + 0.30902f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 12){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.80902f*(*transformed).x + 0.5000f*(*transformed).y + 0.30902f*(*transformed).z;
                (*nextVec).y = -0.5000f*(*transformed).x + 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;
                (*nextVec).z =  0.30902f*(*transformed).x - 0.80902f*(*transformed).y + 0.5000f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 13){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.80902f*(*transformed).x - 0.5000f*(*transformed).y - 0.30902f*(*transformed).z;
                (*nextVec).y =  0.5000f*(*transformed).x + 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;
                (*nextVec).z = -0.30902f*(*transformed).x - 0.80902f*(*transformed).y + 0.5000f*(*transformed).z;

                startCount++;
            }
        }  else if (subunitIndex == 14){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.5000f*(*transformed).x + 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;
                (*nextVec).y =  0.30902f*(*transformed).x + 0.80902f*(*transformed).y + 0.5000f*(*transformed).z;
                (*nextVec).z =  0.80902f*(*transformed).x - 0.5000f*(*transformed).y + 0.30902f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 15){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.5000f*(*transformed).x - 0.30902f*(*transformed).y - 0.80902f*(*transformed).z;
                (*nextVec).y = -0.30902f*(*transformed).x - 0.80902f*(*transformed).y + 0.5000f*(*transformed).z;
                (*nextVec).z = -0.80902f*(*transformed).x + 0.5000f*(*transformed).y + 0.30902f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 16){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.30902f*(*transformed).x - 0.80902f*(*transformed).y - 0.5000f*(*transformed).z;
                (*nextVec).y =  0.80902f*(*transformed).x - 0.5000f*(*transformed).y + 0.30902f*(*transformed).z;
                (*nextVec).z = -0.5000f*(*transformed).x - 0.30902f*(*transformed).y + 0.80902f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 17){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.30902*(*transformed).x - 0.80902*(*transformed).y - 0.5000f*(*transformed).z;
                (*nextVec).y =  0.80902*(*transformed).x + 0.5000f*(*transformed).y - 0.30902*(*transformed).z;
                (*nextVec).z =  0.5000f*(*transformed).x - 0.30902*(*transformed).y + 0.80902*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 18){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.30902*(*transformed).x + 0.80902*(*transformed).y - 0.5000f*(*transformed).z;
                (*nextVec).y = -0.80902*(*transformed).x - 0.5000f*(*transformed).y - 0.30902*(*transformed).z;
                (*nextVec).z = -0.5000f*(*transformed).x + 0.30902*(*transformed).y + 0.80902*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 19){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.5000f*(*transformed).x - 0.30902*(*transformed).y - 0.80902*(*transformed).z;
                (*nextVec).y = -0.30902*(*transformed).x + 0.80902*(*transformed).y - 0.5000f*(*transformed).z;
                (*nextVec).z =  0.80902*(*transformed).x + 0.5000f*(*transformed).y + 0.30902*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 20){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.5000f*(*transformed).x + 0.30902*(*transformed).y - 0.80902*(*transformed).z;
                (*nextVec).y =  0.30902*(*transformed).x - 0.80902*(*transformed).y - 0.5000f*(*transformed).z;
                (*nextVec).z = -0.80902*(*transformed).x - 0.5000f*(*transformed).y + 0.30902*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 21){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.80902*(*transformed).x - 0.5000f*(*transformed).y - 0.30902*(*transformed).z;
                (*nextVec).y =  0.5000f*(*transformed).x - 0.30902*(*transformed).y - 0.80902*(*transformed).z;
                (*nextVec).z =  0.30902*(*transformed).x - 0.80902*(*transformed).y + 0.5000f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 22){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.80902*(*transformed).x + 0.5000f*(*transformed).y - 0.30902*(*transformed).z;
                (*nextVec).y = -0.5000f*(*transformed).x + 0.30902*(*transformed).y - 0.80902*(*transformed).z;
                (*nextVec).z = -0.30902*(*transformed).x + 0.80902*(*transformed).y + 0.5000f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 23){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.80902*(*transformed).x - 0.5000f*(*transformed).y + 0.30902*(*transformed).z;
                (*nextVec).y =  0.5000f*(*transformed).x + 0.30902*(*transformed).y - 0.80902*(*transformed).z;
                (*nextVec).z =  0.30902*(*transformed).x + 0.80902*(*transformed).y + 0.5000f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 24){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.80902*(*transformed).x + 0.5000f*(*transformed).y + 0.30902*(*transformed).z;
                (*nextVec).y = -0.5000f*(*transformed).x - 0.30902*(*transformed).y - 0.80902*(*transformed).z;
                (*nextVec).z = -0.30902*(*transformed).x - 0.80902*(*transformed).y + 0.5000f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 25){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.5000f*(*transformed).x - 0.30902*(*transformed).y + 0.80902*(*transformed).z;
                (*nextVec).y = -0.30902*(*transformed).x - 0.80902*(*transformed).y - 0.5000f*(*transformed).z;
                (*nextVec).z =  0.80902*(*transformed).x - 0.5000f*(*transformed).y + 0.30902*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 26){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.5000f*(*transformed).x + 0.30902*(*transformed).y + 0.80902*(*transformed).z;
                (*nextVec).y =  0.30902*(*transformed).x + 0.80902*(*transformed).y - 0.5000f*(*transformed).z;
                (*nextVec).z = -0.80902*(*transformed).x + 0.5000f*(*transformed).y + 0.30902*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 27){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  1.0f*(*transformed).z;
                (*nextVec).y =  1.0f*(*transformed).x;
                (*nextVec).z =  1.0f*(*transformed).y;

                startCount++;
            }
        } else if (subunitIndex == 28){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  1.0f*(*transformed).y;
                (*nextVec).y =  1.0f*(*transformed).z;
                (*nextVec).z =  1.0f*(*transformed).x;

                startCount++;
            }
        } else if (subunitIndex == 29){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -1.0f*(*transformed).y;
                (*nextVec).y =  1.0f*(*transformed).z;
                (*nextVec).z = -1.0f*(*transformed).x;

                startCount++;
            }
        } else if (subunitIndex == 30){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -1.0f*(*transformed).z;
                (*nextVec).y = -1.0f*(*transformed).x;
                (*nextVec).z =  1.0f*(*transformed).y;

                startCount++;
            }
        } else if (subunitIndex == 31){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -1.0f*(*transformed).y;
                (*nextVec).y = -1.0f*(*transformed).z;
                (*nextVec).z =  1.0f*(*transformed).x;

                startCount++;
            }
        } else if (subunitIndex == 32){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  1.0f*(*transformed).y;
                (*nextVec).y = -1.0f*(*transformed).z;
                (*nextVec).z = -1.0f*(*transformed).x;

                startCount++;
            }
        } else if (subunitIndex == 33){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x = -0.80902f*(*transformed).x - 0.5000f*(*transformed).y + 0.30902f*(*transformed).z;
                (*nextVec).y = -0.5000f*(*transformed).x + 0.30902f*(*transformed).y - 0.80902f*(*transformed).z;
                (*nextVec).z =  0.30902f*(*transformed).x - 0.80902f*(*transformed).y - 0.5000f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 34){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.80902f*(*transformed).x - 0.5000f*(*transformed).y - 0.30902f*(*transformed).z;
                (*nextVec).y = -0.5000f*(*transformed).x - 0.30902f*(*transformed).y - 0.80902f*(*transformed).z;
                (*nextVec).z =  0.30902f*(*transformed).x + 0.80902f*(*transformed).y - 0.5000f*(*transformed).z;

                startCount++;
            }
        } else if (subunitIndex == 35){
            for(unsigned int i=0; i<workingLimit; i++){
                transformed = &coordinates[i];

                nextVec = &coordinates[startCount];
                (*nextVec).x =  0.80902f*(*transformed).x - 0.5000f*(*transformed).y - 0.30902f*(*transformed).z;
                (*nextVec).y = -0.5000f*(*transformed).x - 0.30902f*(*transformed).y - 0.80902f*(*transformed).z;
                (*nextVec).z =  0.30902f*(*transformed).x + 0.80902f*(*transformed).y - 0.5000f*(*transformed).z;

                startCount++;
            }
        }
    }
}

unsigned int PointSetModel::getNumberOfSubUnits() const {
    return numberOfSubUnits;
}


/**
 * returns sorted indices of lattice positions that match input PDB model
 * the output model will likely contain more beads as we are finding all beads within
 * Lattice is based on input pr data file.
 */
void PointSetModel::createSeedFromPDBSym(std::string filename, PofRData * pData, unsigned int totalBins, std::vector<double> * pdbPr){

    PDBModel pdbModel(filename, true, false); // coordinates are centered
    baseWorkingLimit = pdbModel.getTotalCoordinates();

    if(!pdbModel.getEdgeRadiusStatus()){
        throw std::invalid_argument("*** ERROR => PDB FILE MUST INCLUDE EDGE RADIUS SPECIFIED AS: REMARK 265                   EDGE RADIUS : 2.5 ");
    }
    float edge_radius = pdbModel.getEdgeRadius() + bead_radius;
//    float edge_radius = bead_radius;

    const unsigned int totalBeads = this->getTotalNumberOfBeadsInUniverse();

    float diffx, diffy, diffz, bx, by, bz;
    const unsigned int totalAtomsInPDBModel = pdbModel.getTotalCoordinates();

    float beadradius = this->getBeadRadius()*1.0003f;

    std::vector<unsigned int> indicesToCheck(totalAtomsInPDBModel);

    unsigned int * ptr = (totalAtomsInPDBModel != 0) ? &indicesToCheck.front() : nullptr;
    for(unsigned int i = 0; i < totalAtomsInPDBModel; i++) {
        ptr[i] = i;
    }

    logger("RUNNING", "Model::createSeedFromPDBSym");
    logger("","CONVERTING TO LATTICE MODEL");
    logger("TOTAL PTS IN INPUT", std::to_string(totalAtomsInPDBModel));
    logger("TOTAL PTS IN ASYMMETRIC UNIT", std::to_string(totalBeads));
    logger("","...please wait");

    std::vector<unsigned int> temp_indices(bead_indices.begin(), bead_indices.end());  // contains true model of PDB model converted to lattice model
    std::set<unsigned int>uniqueIndices;

    const auto centeringVec = pdbModel.getCenteringVector();

    for(unsigned int i=0; i < totalBeads && totalAtomsInPDBModel>0; i++){ // iterate over each bead in Universe

        Bead * currentBead = this->getBead(i);

        float fromCenter = (currentBead->getVec() - *centeringVec).length();

        for (unsigned int t=0; t < totalAtomsInPDBModel; t++){

            vector3 modelPt (pdbModel.getX()[t], pdbModel.getY()[t], pdbModel.getZ()[t]);

            float len = (currentBead->getVec() - modelPt).length();
            float diffLengthFromCenter = fromCenter - (modelPt - *centeringVec).length();

            if ((diffLengthFromCenter < beadradius || diffLengthFromCenter < pdbModel.getEdgeRadius()) && len < edge_radius){
                uniqueIndices.insert(i);
                break;
            }
        }
    }

    // remove beads that are isolated
    for(auto uni = uniqueIndices.begin(); uni != uniqueIndices.end(); ++uni){
        auto it = this->getPointerToNeighborhood(*uni);
        unsigned int neighborContacts = 0;
        // go through each member of the neighborhood
        for (unsigned int i=0; i< sizeOfNeighborhood; i++){
            unsigned int neighbor = *(it+i);
            if ((neighbor < neighborLimit ) && uniqueIndices.find(neighbor) != uniqueIndices.end()){
                neighborContacts += 1;
            } else if (neighbor == neighborLimit) {
                break;
            }
        }
        if (neighborContacts==0){
            uniqueIndices.erase(uni);
            logger("ISOLATED PT", std::to_string(*uni));
        }
    }
    seed_indices.resize(uniqueIndices.size());
    std::copy(uniqueIndices.begin(), uniqueIndices.end(), seed_indices.begin());


    // make PDB P(r) for refining beadmodel (doesn't matter if centered)
    double sum=0.0;
    for (unsigned int i=0; i<totalAtomsInPDBModel; i++){

        bx = *(pdbModel.getCenteredXVec()+i);
        by = *(pdbModel.getCenteredYVec()+i);
        bz = *(pdbModel.getCenteredZVec()+i);

        unsigned int next = i+1;
        for (unsigned int j=next; j<totalAtomsInPDBModel; j++){
            diffx = bx - *(pdbModel.getCenteredXVec()+j);
            diffy = by - *(pdbModel.getCenteredYVec()+j);
            diffz = bz - *(pdbModel.getCenteredZVec()+j);
            //float dis = sqrt((diffx*diffx + diffy*diffy + diffz*diffz));
            // convert to bin
            (*pdbPr)[pData->convertToBin(sqrt((diffx*diffx + diffy*diffy + diffz*diffz)))]++;
            sum += 1.0;
        }
    }

    //normalize
    double invSum = 1.0/sum;
    for (unsigned int i=0; i<totalBins; i++){
        //std::cout << i << " " << (*pdbPr)[i] << std::endl;
        (*pdbPr)[i] *= invSum;
    }
    // plot should integrate to 1

    std::sort(seed_indices.begin(), seed_indices.end());
    total_in_seed = seed_indices.size();
    //pModel->printSelectedBeads(0, totalKept, *keptBeads);
    //pdbModel.writeCenteredCoordinatesToFile("inputPDBSubUnitCentered_SYM");
    this->writeSubModelToFile(0, total_in_seed, seed_indices, "latticeSubUnit_SYM");
}



/**
 * returns sorted indices of lattice positions that match input PDB model
 */
void PointSetModel::createSeedFromPDB(std::string filename, PofRData * pData, unsigned int totalBins, std::vector<double> * pdbPr){

    PDBModel pdbModel(filename, true, true); // coordinates are centered

    const unsigned int totalBeads = this->getTotalNumberOfBeadsInUniverse();

    float diffx, diffy, diffz, bx, by, bz;
    const unsigned int totalAtoms = pdbModel.getTotalCoordinates();

    const float beadradius = (this->getBeadRadius()*std::sqrt(3.0f/2.0f));

    std::vector<unsigned int> indicesToCheck(totalAtoms);

    unsigned int * const ptr = (totalAtoms != 0) ? indicesToCheck.data() : nullptr;
    for(unsigned int i = 0; i < totalAtoms; i++) {
        ptr[i] = i;
    }

    /*
     * calculate surface accessible volume of the PDB model
     */
    std::vector<float> weights(totalAtoms);
    const float delta_r = 1.4f;
    std::vector<Eigen::Vector3f> coordinates(totalAtoms);
    for (unsigned int i=0; i < totalAtoms; i++){
        bx = *(pdbModel.getCenteredXVec() + ptr[i]);
        by = *(pdbModel.getCenteredYVec() + ptr[i]);
        bz = *(pdbModel.getCenteredZVec() + ptr[i]);

        weights[i] = pdbModel.getAtomicRadius(i) + delta_r;
        coordinates[i] = Eigen::Vector3f(bx,by,bz);
    }

    POWERSASA::PowerSasa<float, Eigen::Vector3f> *ps =
            new POWERSASA::PowerSasa<float, Eigen::Vector3f>(coordinates, weights, 1, 1, 1, 1);

    ps->calc_sasa_all();

    sasa_volume_start = 0.0f;
    sasa_surface_area = 0.0f;
    for (unsigned int i = 0; i < totalAtoms; ++i) {
        sasa_surface_area += ps->getSasa()[i];
        sasa_volume_start += ps->getVol()[i];
    }

    surface_to_volume = ((float)M_PI*(36.0f*sasa_volume_start*sasa_volume_start)/(sasa_surface_area*sasa_surface_area*sasa_surface_area));

    delete ps;
    logger("PDB SASA VOLUME", formatNumber(sasa_volume_start,1));

    std::cout << "** Model::createSeedFromPDB CONVERT TO LATTICE MODEL => please wait " << std::endl;
    std::set<unsigned int> uniqueIndices;

    for(unsigned int i=0; i < totalBeads && totalAtoms>0; i++){ // iterate over each bead in Universe

        Bead * currentBead = this->getBead(i);
        bx = currentBead->getX();
        by = currentBead->getY();
        bz = currentBead->getZ();

        float tempclosest, closest=1.1f*beadradius;
        bool foundIt=false;
        unsigned int found=0;
        for (unsigned int t=0; t < totalAtoms; t++){

            diffx = bx - *(pdbModel.getCenteredXVec() + ptr[t]);
            diffy = by - *(pdbModel.getCenteredYVec() + ptr[t]);
            diffz = bz - *(pdbModel.getCenteredZVec() + ptr[t]);
            tempclosest = std::sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
            if (tempclosest <= closest && uniqueIndices.find(i) == uniqueIndices.end()){
                closest = tempclosest;
                //uniqueIndices.insert(i);
                found = i;
                foundIt = true;
                //break;
            }
        }

        if (foundIt){
            uniqueIndices.insert(found);
        }
    }


    for(auto uni = uniqueIndices.begin(); uni != uniqueIndices.end(); ++uni){
        auto it = this->getPointerToNeighborhood(*uni);
        unsigned int neighborContacts = 0;
        // go through each member of the neighborhood
        for (unsigned int i=0; i< sizeOfNeighborhood; i++){
            unsigned int neighbor = *(it+i);
            if ((neighbor < neighborLimit ) && uniqueIndices.find(neighbor) != uniqueIndices.end()){
                neighborContacts += 1;
            } else if (neighbor == neighborLimit) {
                break;
            }
        }
        if (neighborContacts==0){ // remove any lattice points that are isolated
            uniqueIndices.erase(uni);
        }
    }

    seed_indices.resize(uniqueIndices.size());
    std::copy(uniqueIndices.begin(), uniqueIndices.end(), seed_indices.begin());
    logger("TOTAL UNIQUE INDICES", formatNumber((unsigned int)seed_indices.size()));
    logger("TOTAL IN UNIVERSE", formatNumber(totalBeads));
    // make PDB P(r) for refining beadmodel
    double sum=0.0;
    for (unsigned int i=0; i<totalAtoms; i++){

        bx = *(pdbModel.getCenteredXVec()+i);
        by = *(pdbModel.getCenteredYVec()+i);
        bz = *(pdbModel.getCenteredZVec()+i);

        unsigned int next = i+1;
        for (unsigned int j=next; j<totalAtoms; j++){
            diffx = bx - *(pdbModel.getCenteredXVec()+j);
            diffy = by - *(pdbModel.getCenteredYVec()+j);
            diffz = bz - *(pdbModel.getCenteredZVec()+j);
            float dis = sqrt((diffx*diffx + diffy*diffy + diffz*diffz));
            // convert to bin
            (*pdbPr)[pData->convertToBin(dis)]++;
            sum += 1.0;
        }
    }
    std::cout << " A " << std::endl;
    //normalize
    double invSum = 1.0d/sum;
    for (unsigned int i=0; i<totalBins; i++){
        //std::cout << i << " " << (*pdbPr)[i] << std::endl;
        (*pdbPr)[i] *= invSum;
    }

    std::cout << " AA " << std::endl;
    // plot should integrate to 1
    // prune beads that are single
    EulerTour eulerTour(seed_indices.begin(), seed_indices.size(), this);
    if (eulerTour.getNumberOfComponents() > 1){ // add neighbors

    }
    std::cout << " AAA " << std::endl;
    std::sort(seed_indices.begin(), seed_indices.end());

    total_in_seed = seed_indices.size();
    //pModel->printSelectedBeads(0, totalKept, *keptBeads);
    std::cout << " AAAA " << std::endl;
    pdbModel.writeCenteredCoordinatesToFile("centeredPDBSeed");

    this->writeSubModelToFile(0, total_in_seed, seed_indices, "centeredLatticeSeed");
}

std::string PointSetModel::writeModelToFileBare(float dkl, unsigned int workingNumber, std::vector<unsigned int> &selectedBeads, std::vector<unsigned int> &pofrModel, std::string nameOf, Anneal *annealedObject, unsigned int steps, float volume, float averageContacts){

    const char *outputFileName;
    nameOf = nameOf + ".pdb";
    outputFileName = nameOf.c_str() ;
    FILE * pFile = fopen(outputFileName, "w");

    unsigned int residue_index;

    // Add P(r) distributions
    // final D_kl, volume and energy
    float totalCounts = 0.0;

    auto totalm = (unsigned int)pofrModel.size();

    // normalization constant of model Pr
    // treats each value as discrete (i.e. not integrating via trapezoid)
    for (unsigned int i=0; i<totalm; i++){
        totalCounts += pofrModel[i];
    }

//    unsigned int shannon_bins = annealedObject->gettotalBinsDerivedFromData();
//    for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
//    fprintf(pFile, annealedObject->getConnectivityTableText().c_str());
    fprintf(pFile, "REMARK 265\n");
    // write coordinates
    for (unsigned int i=0; i<workingNumber; i++){
        Bead * currentBead = this->getBead(selectedBeads[i]);
        //residue_index = std::to_string(selectedBeads[i]);
        residue_index = i + 1;
        //residue_index = std::to_string(selectedBeads[i]);
        //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        //fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", i+1," CA ", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        printAtomLine(pFile, i+1, "A", residue_index, currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    }

    fprintf(pFile, "END\n");
    fclose(pFile);
    return nameOf;
}