/*
 * extractor.cpp
 *
 *  Created on: 14 Sep 2017
 *      Author: jkiesele
 */

#include "../include/extractor.h"
#include "../include/simpleFitter.h"
#include <stdexcept>
#include "../include/globals.h"


TGraph * extractResult::createScanGraph()const{
	return 0;
}




prediction * extractor::readPrediction(const std::string& file){
    pred.readFromFile(file);
    return &pred;
}

extractResult extractor::extractPriv(bool ismtop){
	simpleFitter fitter;
	size_t freepara=0;
	for(size_t i=0;i<priors_.size();i++){
		if(priors_.at(i)==prior_free){
			freepara=i;
			break;
		}
	}
    std::vector<double> startparas={
            globals::default_mtop,
            globals::default_alphas,
            0,0};
    fitter.setParameters(startparas);
    ROOT::Math::Functor func(this,&extractor::toBeMinimized,startparas.size());
    fitter.setMinFunction(func);
  //  fitter.setAsMinosParameter(freepara,true);
    //set parameter limits for limted paras FIXME

    fitter.fit();

    if(!fitter.wasSuccess())
        throw std::runtime_error("fit not successful");


    extractResult res;

    res.central =fitter.getParameter(freepara);
    res.errup   =fitter.getParameterErrUp()->at(freepara);
    res.errdown =fitter.getParameterErrDown()->at(freepara);

    //make a scan in 10 points between
    //FIXME
    // use ismtop
    std::vector<float> scanpoints;
    float lower=globals::default_mtop_low;
    float higher=globals::default_mtop_high;
    if(!ismtop){
        lower=globals::default_alphas_low;
        higher=globals::default_alphas_high;
    }



    for(float m=lower;m<=higher+0.0001;m+= (higher-lower)/50){

        break;
        std::cout << m << std::endl;


        fitter.setMinFunction(func);
        fitter.setParameter(freepara,m);
        fitter.setStrategy(1);//enough for min chi2
        fitter.setParameterFixed(freepara,true);
        fitter.fit();
        if(!fitter.wasSuccess()){
            continue;
            throw std::runtime_error("fit not successful");
        }
        float minchi2=fitter.getChi2Min();

        res.scan.push_back(std::pair<float,float>(m,minchi2));
    }

    return res;
}

extractResult extractor::extractMtop()const{
    extractor cpex=*this;
    cpex.priors_.clear();
    // mtop, alphas, scale, pdf
    cpex.priors_.push_back(prior_free);
    cpex.priors_.push_back(prior_gauss);
    cpex.priors_.push_back(scaleprior_);
    cpex.priors_.push_back(prior_gauss);

    return cpex.extractPriv(true);
}



double extractor::toBeMinimized(const double* pars)const{

    //the actual likelihood
    //mtop,alphas,scale,pdf
    const double& mtop=pars[0];
    const double& alphas=pars[1];
    double scale=pars[2];
    const double& pdf=pars[3];

    if(priors_.at(2) == prior_fixed1)
    	scale=1;
    else if(priors_.at(2) == prior_fixedm1)
    	scale=-1;
    else if(priors_.at(2) == prior_fixed0)
    	scale=0;

    double delta7=globals::measured_xsec-pred.eval(mtop,alphas,scale,0);
    double err7=globals::measured_errdown;
    if(delta7<0)
    	err7=globals::measured_errup;
    double pdferr=pred.getRelPdfUp();
    if(delta7<0)
    	pdferr=pred.getRelPdfDown();
    pdferr*=pred.eval(mtop,alphas,scale,0);

    pdferr*=pdferr;
    err7*=err7;

    double out= (delta7*delta7)/(pdferr+err7);


    for(size_t i=0;i<priors_.size();i++){
        if(priors_.at(i) == prior_gauss){
            if(i==para_mtop){
                double delta=pars[i]-globals::default_mtop;
                delta/=globals::default_mtop_err;
                out+=delta*delta;
            }
            else if(i==para_alphas){
                double delta=pars[i]-globals::default_alphas;
                delta/=globals::default_alphas_err;
                out+=delta*delta;
            }
            else
                out+=pars[i]*pars[i];
        }
        else if (priors_.at(i) == prior_fixed0 || priors_.at(i) == prior_fixed1 || priors_.at(i) == prior_fixedm1){
        	out+=pars[i]*pars[i];
        }
    }
    return out;
}







