/*
 * genarateDistributions.cpp
 *
 *  Created on: 23 mars 2014
 *      Author: Olivier Coulaud
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Utils/FGlobal.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FPoint.hpp"
#include "Files/FGenerateDistribution.hpp"
#include "Files/FFmaGenericLoader.hpp"
#include "Files/FExportWriter.hpp"

#include "Utils/FParameterNames.hpp"

/**
 * \file
 *
 * \brief Generates points (non)uniformly distributed on a given geometry
 *
 * The goal of this driver is to generate uniform or non uniform points on the
 * following geometries
 *
 *   - Uniform : cube, cuboid, sphere, prolate,
 *   - Non uniform : ellipsoid, prolate
 *
 *  You can set two kind of physical values depending of your problem. By
 *   default all values are between 0 and 1.  If you select the argument -charge
 *   (see bellow) the values are between -1 and 1.  The arguments available are
 *
 * <b> General arguments:</b>
 * \param -help (-h)      to see the parameters available in this driver
 * \param -N     The number of points in the distribution (default 20000)
 * \param -fout name: generic name for files (with extension) and save data with
 *                   following format in name.fma or name.bfma in -bin is set"
 * \param -fvisuout Filename for the visu file (vtk, vtp, cvs or cosmo). vtp is
 *                  the default
 * \param -extraLength value extra length to add to the boxWidth (default 0.0)
 * <b> Geometry arguments:</b>
 * \param -unitCube uniform distribution in unit cube
 * \param -cube uniform distribution in a cube
 *     \arg -size LX:LY:LZ - default value for R is 1.0:1.0:2.0
 * \param -unitSphere uniform distribution on unit sphere
 * \param -sphere uniform distribution on sphere of radius given by
 *     \arg -radius R - default value for R is 2.0
 * \param -ball uniform distribution in ball of radius given by
 *     \arg -radius R - default value for R is 2.0
 * \param -ellipsoid non uniform distribution on an ellipsoid of aspect ratio
 *                   given by
 *     \arg -size a:b:c with a, b and c > 0
 * \param -prolate ellipsoid with aspect ratio a:a:c given by
 *     \arg -size a:a:c with c > a > 0
 * \param -plummer (Highly non uniform) plummer distribution (astrophysics)
 *     \arg -radius R - default value 10.0"
 *
 *
 * <b> Physical values argument:</b>
 * \param -charge generate physical values between -1 and 1 otherwise generate between 0 and 1
 * \param -zeromean  the average of the physical values is zero
 *
 *
 * <b> examples</b>
 *
 *    generateDistributions -prolate -size 2:2:4   -N 20000 -fout prolate
 *
 *  or
 *
 *    generateDistributions -cuboid 2:2:4 -N 100000 -fout cuboid.bfma  -fvisuout cuboid.vtp -charge  -zeromean
 *
 */


namespace Param {
    const FParameterNames Ellipsoid
    = {{"-ellipsoid"}, "non uniform distribution on an ellipsoid of aspect ratio given by -size a:b:c with a, b and c > 0"};
    const FParameterNames UnitCube
    = {{"-unitCube"}, "uniform distribution on unit cube"};
    const FParameterNames Cube
    = {{"-cuboid"}, "uniform distribution on rectangular cuboid of size -size a:b:c - default values are 1.0:1.0:2.0 "};
    const FParameterNames UnitSphere
    = {{"-unitSphere"}, "uniform distribution on unit sphere"};
    const FParameterNames Ball
    = {{"-ball"}, "uniform distribution in a ball of radius given by -radius R - default value for R is 2.0"};
    const FParameterNames Sphere
    = {{"-sphere"}, "uniform distribution on sphere of radius given by -radius R - default value for R is 2.0"};
    const FParameterNames Prolate
    = {{"-prolate"}, "ellipsoid with aspect ratio a:a:c given by -size a:a:c with c > a > 0"};
    const FParameterNames Plummer
    = {{"-plummer"}, "(Highly non uniform) plummer distribution (astrophysics) -radius R - default value 10.0"};
    const FParameterNames Size
    = {{"-size"}, "Size of the geometry a:b:c - default values are 1.0:1.0:2.0"};
    const FParameterNames Radius
    = {{"-radius"}, "used to specified the radius of the sphere and the plummer distribution or R - default value for R is 2.0"};
    const FParameterNames Charge
    = {{"-charge"}, "generate physical values between -1 and 1 otherwise generate between 0 and 1"};
    const FParameterNames ZM
    = {{"-zeromean"}, "the average of the physical values is zero"};
    const FParameterNames EL
    = {{"-extraLength"}, "-extraLength value extra length to add to the boxWidth"};
}

#define getParamV(name, default)                                \
    FParameters::getValue(argc,argv,(name).options,(default))

#define getParamS(name, default)                                \
    FParameters::getStr(argc,argv,(name).options,(default))

int main(int argc, char ** argv){

    FHelpDescribeAndExit(
        argc, argv,
        ">> Driver to generate N points (non)uniformly distributed on a given geometry.\n"
        "Options  \n"
        "   -help       to see the parameters    ",
        FParameterDefinitions::OutputFile, FParameterDefinitions::NbParticles,
        FParameterDefinitions::OutputVisuFile,
        Param::UnitCube, Param::Cube,      Param::UnitSphere, Param::Sphere,
        Param::Radius,   Param::Ellipsoid, Param::Prolate,    Param::Plummer,
        Param::Ball,
        Param::Charge,   Param::ZM,        Param::EL,    Param::Size
        );

    using FReal = double;

    const FSize NbPoints = getParamV(FParameterDefinitions::NbParticles, FSize(20000));

    FReal extraRadius = 0.000;
    FReal BoxWith = 0.0;
    FPoint<FReal> Centre(0.0, 0.0,0.0);

    // Allocate particle array
    FReal * particles;
    particles = new FReal[4*NbPoints] ;
    memset(particles, 0, 4*NbPoints*sizeof(FReal));
    FmaRWParticle<FReal, 4, 4>* ppart = (FmaRWParticle<FReal, 4, 4>*)(&particles[0]);

    // Generate physical values
    FReal sum = 0;
    FReal a = 1.0;
    FReal b = 0.0;
    if(FParameters::existParameter(argc, argv, "-charge")){
        a = 2.0; b = -1.0;
    }

    for(int i = 0, j = 3 ; i< NbPoints; ++i, j+=4){
        particles[j] = a * getRandom<FReal>() + b;
        sum += particles[j] ;
    }

    if(FParameters::existParameter(argc, argv, "-zeromean")){
        FReal rm = FReal(sum) / FReal(NbPoints) ;
        sum -= static_cast<FReal>(NbPoints) * rm;
        for(int i = 0, j = 3 ; i< NbPoints; ++i, j+=4){
            particles[j] -= rm ;
        }
    }

    std::cout << "Physical value sum: " << sum
              << " mean: " << sum / FReal(NbPoints)
              << std::endl;

    // Read arguments
    // Radius
    const FReal Radius = getParamV(Param::Radius,  2.0);
    // Aspect ratio
    std::string aspectRatio = getParamS(Param::Size, "1:1:2");
    std::replace(aspectRatio.begin(), aspectRatio.end(), ':', ' ');
    FReal A, B, C;
    std::stringstream(aspectRatio) >> A >> B >> C;

    // Point  generation
    if(FParameters::existParameter(argc, argv, "-unitCube")) {
        unifRandomPointsInCube<FReal>(NbPoints, 1, 1, 1, particles);
        Centre.setPosition(0.5,0.5,0.5);
        BoxWith = 1.0;
        std::cout << "Unit cube "<< std::endl;
    }
    else if(FParameters::existParameter(argc, argv, "-ball")) {
        unifRandomPointsInBall<FReal>(NbPoints, Radius, particles);
        BoxWith = 2.0 * Radius;
        std::cout << "Ball radius: " << Radius << std::endl;
    }
    else if(FParameters::existParameter(argc, argv, "-cuboid")) {
        unifRandomPointsInCube(NbPoints, A, B, C, particles);
        BoxWith = FMath::Max(A, FMath::Max(B,C));
        FReal halfBW = BoxWith * 0.5;
        Centre.setPosition(halfBW, halfBW, halfBW);
        std::cout << "Cuboid: "<< A << ":" << B << ":" << C << std::endl;
    }
    else if(FParameters::existParameter(argc, argv, "-unitSphere")) {
        unifRandomPointsOnSphere<FReal>(NbPoints, 1.0, particles);
        BoxWith = 2.0;
    }
    else if(FParameters::existParameter(argc, argv, "-sphere")) {
        unifRandomPointsOnSphere(NbPoints, Radius, particles);
        BoxWith = 2.0 * Radius;
        std::cout << "Sphere radius: " << Radius << std::endl;
    }
    else if(FParameters::existParameter(argc, argv, "-prolate")) {
        if(A != B){
            std::cerr << " A != B in prolate ellipsoid. Your aspect ratio: "
                      << aspectRatio << std::endl;
        }
        std::cout << "Prolate A: " << A << " B: " << B << " C: " << C << std::endl;
        unifRandomPointsOnProlate(NbPoints, A, C, particles);
        BoxWith = 2.0 * C;
    }
    else if(FParameters::existParameter(argc, argv, "-hyperpara")) {
        unifRandomPointsOnHyperPara(NbPoints, A, B, C, particles);
        BoxWith = 2.0 * FMath::Max(A, FMath::Max(B, C));
        std::cout << "Hyperpara "<< A << ":"<< B<<":"<<C<<std::endl;
        std::cout << "BoxWith: " << BoxWith << std::endl;

    }
    else if(FParameters::existParameter(argc, argv, "-ellipsoid")){
        nonunifRandomPointsOnElipsoid(NbPoints, A, B, C, particles);
        BoxWith =  2.0 * FMath::Max(A, FMath::Max(B, C));
        std::cout << "Ellipsoid " << A << ":" << B << ":" << C << std::endl;
    }
    else if(FParameters::existParameter(argc, argv, "-plummer")){
        unifRandomPlummer(NbPoints, Radius, particles);
        BoxWith = 2.0 * Radius;
        std::cout << "Plummer radius: " << Radius << std::endl;
    }
    else {
        std::cout << "Bad geometry option"<< std::endl;
        exit(-1);
    }

    /////////////////////////////////////////////////////////////////////////
    //                                           Save data
    /////////////////////////////////////////////////////////////////////////
    //  Generate FMA file for FFmaGenericLoader<FReal> Loader
    if(FParameters::existParameter(argc, argv, "-extraLength")){
        extraRadius = FParameters::getValue(argc, argv, "-extraLength",  0.0);
        BoxWith += 2 * extraRadius;
    }
    const std::string name(getParamS(FParameterDefinitions::OutputFile, "unifPointDist"));
    std::cout << "Write "<< NbPoints <<" particles to '" << name << "'" << std::endl;
    FFmaGenericWriter<FReal> writer(name);
    writer.writeHeader(Centre, BoxWith, NbPoints, *ppart);
    writer.writeArrayOfParticles(ppart, NbPoints);
    std::cout << "End of writing" <<std::endl;

    //  Generate  file for visualization
//    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputVisuFile.options)){
//        std::string outfilename(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output.vtp"));
//        driverExportData(outfilename, particles , NbPoints,loader.getNbRecordPerline() );
//    }
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputVisuFile.options)) {
        std::string visufile(FParameters::getStr(argc, argv, FParameterDefinitions::OutputVisuFile.options, "output.vtp"));
        driverExportData(visufile, particles , NbPoints);
    }
    //
    delete [] particles;
}
