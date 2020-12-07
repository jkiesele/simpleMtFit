/*
 * extractor.h
 *
 *  Created on: 14 Sep 2017
 *      Author: jkiesele
 */

#ifndef TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_EXTRACTOR_H_
#define TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_EXTRACTOR_H_


#include "prediction.h"
#include "simpleFitter.h"
#include <string>
#include <algorithm>
#include <vector>
#include "TGraphAsymmErrors.h"

class extractResult{
public:
    double central;
    double errup;
    double errdown;
    TString pdfname;
    std::vector<std::pair<double,double> > scan;

    TGraph * createScanGraph()const;
};

class extractor{
public:

	extractor():scaleprior_(prior_gauss), pdfscale_(1.){}

    prediction * readPrediction(const std::string& file);

    extractResult extractMtop()const;

    const prediction& getPred()const{return pred;}

    enum priors{prior_gauss,prior_box,prior_free,prior_fixed1,prior_fixedm1,prior_fixed0};

    void scalePDFUnc(double scale=1){pdfscale_=scale;}

    void fixScaleVar(priors p){
    	scaleprior_=p;
    }

private:
    double toBeMinimized(const double* pars)const;
    extractResult extractPriv(bool ismtop);

    prediction pred;

    //temps for the fit
    enum paranumbers{para_mtop=0,para_alphas=1,para_scale=2,para_pdf=3};
    priors scaleprior_;
    std::vector<priors> priors_;
    double pdfscale_;

};


#endif /* TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_EXTRACTOR_H_ */
