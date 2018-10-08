/*
 * plotMaker.h
 *
 *  Created on: 15 May 2018
 *      Author: jkiesele
 */

#ifndef TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_PLOTMAKER_H_
#define TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_PLOTMAKER_H_

#include "prediction.h"
#include "extractor.h"

class plotMaker{
public:
    plotMaker():cmsandatlas_(false){}
    ~plotMaker(){}

    void setExtractor(const extractor&m){extractor_=m;}

    void addPrediction(const prediction&p7,
            const prediction&p8,
            const extractResult & mtopres,
            const extractResult & alphasres,
            const TString& name){
        predictions7_.push_back(p7);
        predictions8_.push_back(p8);
        results_mtop.push_back(mtopres);
        results_alphas.push_back(alphasres);
        names.push_back(name);
    }


    void makePlotAlphaS()const;
    void makePlotMtop()const;

    void plotCMSandATLAS(bool set){
        cmsandatlas_=set;
    }

private:


    extractor extractor_;
    std::vector<prediction> predictions7_,predictions8_;
    std::vector<extractResult> results_mtop,results_alphas;
    std::vector<TString> names;

    bool cmsandatlas_;

};


#endif /* TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_PLOTMAKER_H_ */
