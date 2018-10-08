/*
 * prediction.h
 *
 *  Created on: 12 Sep 2017
 *      Author: jkiesele
 */

#ifndef TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_PREDICTION_H_
#define TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_PREDICTION_H_

#include <string>
#include <vector>
#include "TGraphAsymmErrors.h"

/**
 *
 * pseudo-2D dependence on mtop and alphas
 * both are rather independent so pseudo-2D is enough
 *
 */
class prediction{
public:
    prediction():linestyle(1),linewidth(2),linecolor(kBlack),
        rel_scaleup_(0),rel_scaledown_(0),rel_pdfup_(0),rel_pdfdown_(0),nominal_(0){}

    /*
     * mtop and alpha_s are real values while pdf and scale variations are:
     * -1 .. 1, with 0 being the nominal
     */
    double eval(const double& mtop,const double& alphaas,
            const double& scale, const double& pdf)const;

    double getEnvelope(const double& mtop,const double& alphas, bool isMtop, bool isUp)const;
    /*
     * File format:
     *
     * # comment
     *
     * [nominal]
     * <nominal> <scale up> <scale down> <pdf up> <pdf down>
     * [end nominal]
     *
     * [mtop]
     * <top mass> <xsec>
     * ...
     * [end mtop]
     *
     * [alpha_s]
     * <alpha_s> <xsec>
     * ...
     * [end alpha_s]
     * -----
     * the values for mtop are given for alpha_s=0.118, the values for alpha_s at mtop=172.5
     */
    void readFromFile(const std::string& filename);

    TGraphAsymmErrors* createGraphMtop(int color, int style=1, int width=2)const;
    TGraphAsymmErrors* createGraphAlphas(int color, int style=1, int width=2)const;

    TGraphAsymmErrors* createGraphMtop()const;
    TGraphAsymmErrors* createGraphAlphas()const;

    int linestyle,linewidth,linecolor;
    TString name;

    bool setup()const{return nominal_;}


private:

    double getDeltaAlphaS(const double& alpha_s,const double* paras)const;
    double getDeltaMtop(const double& mt,const double* paras)const;
    double evalparas(const double& x,const double& xnom,const double* paras)const;

    double minFunction(const double* pars)const;
    void fitLine(double* setpars);

    TGraphAsymmErrors*  createGraph(bool mtop)const;

    static const size_t nparas_=4;

    double alphas_paras_[nparas_];
    double mtop_paras_[nparas_];

    double rel_scaleup_,rel_scaledown_,rel_pdfup_,rel_pdfdown_;

    double nominal_;//for mtop=172.5, alphas=0.118


    //intermediate for fit
    std::vector<std::pair<double,double> > fitpoints_;
    double fitnominal_;



};


#endif /* TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_PREDICTION_H_ */
