/*
 * massFit.cpp
 *
 *  Created on: 13 Sep 2017
 *      Author: jkiesele
 */




#include "../include/extractor.h"
#include "../include/globals.h"
#include "TString.h"
#include <string>
#include <iostream>
#include "../include/plotMaker.h"
#include "../include/fileReader.h"

int main(int argc, char* argv[]){



	//fileReader::debug=true;
	if(argc<2){
		return -1;
	}

	std::string predfile = argv[1];

	std::cout << "PDF " << predfile <<std::endl;
	extractor tmp_ex;

	std::cout << "reading prediction"<<std::endl;
	tmp_ex.readPrediction(predfile);


	std::cout << "mtop pred scan: ";
	for(float mtop=166.5;mtop<176;mtop+=2){
		std::cout << mtop<<":"<<tmp_ex.getPred().eval(mtop,globals::default_alphas,0,0)<<" ";
	}


	extractResult resmtop=tmp_ex.extractMtop();
	std::cout << "scale Gauss: "<< resmtop.central << " + " << resmtop.errup << " - " << resmtop.errdown << std::endl;
	tmp_ex.fixScaleVar(extractor::prior_fixed0);
	resmtop=tmp_ex.extractMtop();
	std::cout << "scale 0: "<< resmtop.central << " + " << resmtop.errup << " - " << resmtop.errdown << std::endl;
	tmp_ex.fixScaleVar(extractor::prior_fixedm1);
	resmtop=tmp_ex.extractMtop();
	std::cout << "scale -1: "<< resmtop.central << " + " << resmtop.errup << " - " << resmtop.errdown << std::endl;
	tmp_ex.fixScaleVar(extractor::prior_fixed1);
	resmtop=tmp_ex.extractMtop();
	std::cout << "scale 1: "<< resmtop.central << " + " << resmtop.errup << " - " << resmtop.errdown << std::endl;

	return 0;

}
