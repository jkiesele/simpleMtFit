/*
 * prediction.cpp
 *
 *  Created on: 12 Sep 2017
 *      Author: jkiesele
 */




#include "../include/prediction.h"
#include "../include/fileReader.h"
#include "../include/simpleFitter.h"
#include "../include/globals.h"
#include <stdexcept>
#include <algorithm>

//defaults


double prediction::getDeltaAlphaS(const double& alpha_s,const double* paras)const{
    return  evalparas(alpha_s,globals::default_alphas,paras);
}
double prediction::getDeltaMtop(const double& mt,const double* paras)const{
    return  evalparas(mt,globals::default_mtop,paras);
}
double prediction::evalparas(const double& xin,const double& xnom,const double* paras)const{
    double x=xin-xnom;
    double ret=0;
    double xmult=1;
    for(int i=0;i<nparas_;i++){
        ret+=paras[i]*xmult;
        xmult*=x;
    }
    return ret;
}

double prediction::eval(const double& mtop,const double& alphas,
        const double& scale, const double& pdf)const{

    double central=nominal_+getDeltaMtop(mtop,mtop_paras_)+getDeltaAlphaS(alphas,alphas_paras_);

    double deltascale=scale*rel_scaleup_*central;
    if(scale<0)
        deltascale=-scale*rel_scaledown_*central;


    double deltapdf=pdf*rel_pdfup_*central;
    if(pdf<0)
        deltapdf=-pdf*rel_pdfdown_*central;


    return central+deltascale+deltapdf;
}

double prediction::getEnvelope(const double& mtop,const double& alphas, bool isMtop, bool isUp)const{

    double nom=eval(mtop,alphas,0,0);

    double scaleup=eval(mtop,alphas,1,0)-nom;
    double scaledown=eval(mtop,alphas,-1,0)-nom;
    double pdfup=eval(mtop,alphas,0,1)-nom;
    double pdfdown=eval(mtop,alphas,0,-1)-nom;

    double maxpdf=std::max(pdfup,pdfdown);
    double minpdf=std::min(pdfup,pdfdown);
    double maxscale=std::max(scaleup,scaledown);
    double minscale=std::min(scaleup,scaledown);

    double mtop_as_var_up=0;
    double mtop_as_var_down=0;

    if(isMtop){
        mtop_as_var_up = getDeltaAlphaS(alphas+globals::default_alphas_err,alphas_paras_);
        mtop_as_var_down = getDeltaAlphaS(alphas-globals::default_alphas_err,alphas_paras_);
    }
    else{
        mtop_as_var_up = getDeltaMtop(mtop-globals::default_mtop_err,mtop_paras_);
        mtop_as_var_down = getDeltaMtop(mtop+globals::default_mtop_err,mtop_paras_);
    }
    if(isUp){
        return std::sqrt(maxpdf*maxpdf + maxscale*maxscale + mtop_as_var_up*mtop_as_var_up);
    }
    else{
        return std::sqrt(minpdf*minpdf + minscale*minscale + mtop_as_var_down*mtop_as_var_down);
    }

}

void prediction::readFromFile(const std::string& filename){
    fileReader fr;
    fr.setComment("#");

    fr.setStartMarker("[global]");
    fr.setEndMarker("[end global]");
    fr.setDelimiter(" ");

    fr.readFile(filename);
    name=fr.getValue<TString>("name");

    fr.setStartMarker("[nominal]");
    fr.setEndMarker("[end nominal]");


    fr.readFile(filename);

    nominal_=0;
    for(size_t l=0;l<fr.nLines();l++){
        if(fr.nEntries(l) > 4){
            nominal_=       fr.getData<double>(l,0);
            rel_scaleup_=   fr.getData<double>(l,1);
            rel_scaledown_= fr.getData<double>(l,2);
            rel_pdfup_=     fr.getData<double>(l,3);
            rel_pdfdown_=   fr.getData<double>(l,4);
            break;
        }
    }
    if(!nominal_){
        throw std::runtime_error("prediction::readFromFile: no nominal entry found (check file)");
    }

    fr.setStartMarker("[mtop]");
    fr.setEndMarker("[end mtop]");
    fr.readFile(filename);

    fitpoints_.clear();
    for(size_t l=0;l<fr.nLines();l++){
        if(fr.nEntries(l)<2)continue;
        fitpoints_.push_back(std::pair<double,double>(fr.getData<double>(l,0), fr.getData<double>(l,1) ) );
    }
    fitnominal_=globals::default_mtop;

    //fit the parameters:
    fitLine(mtop_paras_);

    fr.setStartMarker("[alpha_s]");
    fr.setEndMarker("[end alpha_s]");
    fr.readFile(filename);

    fitpoints_.clear();
    for(size_t l=0;l<fr.nLines();l++){
        if(fr.nEntries(l)<2)continue;
        fitpoints_.push_back(std::pair<double,double>(fr.getData<double>(l,0), fr.getData<double>(l,1) ) );
    }
    fitnominal_=globals::default_alphas;

    fitLine(alphas_paras_);

    fitpoints_.clear();

}

TGraphAsymmErrors*  prediction::createGraph(bool mtop)const{
    const int npoints=50;

    double x[npoints];
    double xerr[npoints];
    double y[npoints];
    double yerrup[npoints];
    double yerrdown[npoints];


    float rangemin=globals::default_alphas_low;
    float rangemax=globals::default_alphas_high;
    if(mtop){
      rangemin=globals::default_mtop_low;
      rangemax=globals::default_mtop_high;
    }

    const float step=(rangemax-rangemin)/(float)(npoints-1);

    for(int i=0;i<npoints;i++){
        x[i]=rangemin+step*(float)i;
        xerr[i]=0;
        if(mtop){
            y[i]=eval(x[i],globals::default_alphas,0,0);
            yerrup[i]  = getEnvelope(x[i],globals::default_alphas,true,true);
            yerrdown[i]= getEnvelope(x[i],globals::default_alphas,true,false);
        }
        else{
            y[i]=eval(globals::default_mtop,x[i],0,0);
            yerrup[i]  = getEnvelope(globals::default_mtop,x[i],false,true);
            yerrdown[i]= getEnvelope(globals::default_mtop,x[i],false,false);
        }


    }
    TGraphAsymmErrors * g = new TGraphAsymmErrors(npoints,x,y,xerr,xerr,yerrdown,yerrup);
    g->SetName(name);
    g->SetTitle(name);
    return g;
}

TGraphAsymmErrors* prediction::createGraphMtop(int lcolor, int lstyle, int lwidth)const{
    TGraphAsymmErrors* g= createGraph(true);
    g->SetLineColor(lcolor);
    g->SetLineWidth(lwidth);
    g->SetLineStyle(lstyle);
    return g;
}

TGraphAsymmErrors* prediction::createGraphAlphas(int lcolor, int lstyle, int lwidth)const{
    TGraphAsymmErrors* g=  createGraph(false);
    g->SetLineColor(lcolor);
    g->SetLineWidth(lwidth);
    g->SetLineStyle(lstyle);
    return g;
}

TGraphAsymmErrors* prediction::createGraphMtop()const{
    return createGraphMtop(linecolor, linestyle, linewidth);
}
TGraphAsymmErrors* prediction::createGraphAlphas()const{
    return createGraphAlphas(linecolor, linestyle, linewidth);
}

void prediction::fitLine(double* setpars){
    simpleFitter fitter;
    std::vector<double> startparas(nparas_,0);
    fitter.setParameters(startparas);
    ROOT::Math::Functor func(this,&prediction::minFunction,nparas_);
    fitter.setMinFunction(func);
    fitter.fit();
    for(size_t i=0;i<nparas_;i++){
        setpars[i]=fitter.getParameters()->at(i);
    }
}

double prediction::minFunction(const double* pars)const{

    double out=0;
    for(const auto& p:fitpoints_){
        double delta=evalparas(p.first,fitnominal_,pars)+nominal_-p.second;
        delta/=p.second;
        delta*=delta;
        out+=delta;
    }
    return out;

}




