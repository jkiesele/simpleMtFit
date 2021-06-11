/*
 * massFit.cpp
 *
 *  Created on: 13 Sep 2017
 *      Author: jkiesele
 */




#include "../include/extractor.h"
#include "../include/globals.h"
#include "../include/fileReader.h"
#include "TString.h"
#include <string>
#include <iostream>
#include "../include/plotMaker.h"
#include "../include/fileReader.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TColor.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGAxis.h"
#include "TStyle.h"
#include "TLegendEntry.h"
#include "TPad.h"


class as_extraction{
public:
    as_extraction(){read();}

    void read(std::string filename="data/as_bril.dat"){
        lumi.clear();
        as.clear();
        exp_unc.clear();
        pdf_unc.clear();
        scale_unc.clear();

        fileReader fr;
        fr.setDelimiter(" ");
        fr.setComment("lumi"); //remove the first line
        fr.readFile(filename);

        for(size_t l=0;l<fr.nLines();l++){
            lumi.push_back(fr.getData<double>(l,0));
            double asval = fr.getData<double>(l,1);

            as.push_back(asval);
            double exp = fr.getData<double>(l,2);
            exp_unc.push_back(exp/asval*100.);

            //this is not really correct, but to get a feeling...
            double exppdf = fr.getData<double>(l,3);
            double pdf = std::sqrt(exppdf*exppdf  - exp*exp);

            pdf_unc.push_back(pdf/asval*100.);

            double scale = fr.getData<double>(l,4);
            if(fr.getData<double>(l,5)> scale)
                scale = fr.getData<double>(l,5);

            std::cout << "scale " << scale <<" pdf "<< pdf<< std::endl;

            scale_unc.push_back(scale/asval*100.);

            std::cout << fr.getData<double>(l,0) << " " << asval << " " << exp << " " << exppdf << " " << scale << " " << std::endl;
            std::cout << fr.getData<double>(l,0) << " " << asval << " " << exp << " " << pdf << " " << scale << " "
                    << std::endl << std::endl;
        }
    }

    TGraph *  createGraph(TString name, double scalepdf=1, double scalescale=1, double scaletotal=1){
        std::vector<double> x,y;
        for(size_t i=0;i<lumi.size();i++){
            x.push_back(lumi.at(i));

            double uncsq = exp_unc.at(i)*exp_unc.at(i)
                    + scalepdf*scalepdf*pdf_unc.at(i)*pdf_unc.at(i)
                    + scalescale*scalescale*scale_unc.at(i)*scale_unc.at(i);
            y.push_back(std::sqrt(uncsq)*scaletotal);

        }
        TGraph * g = new TGraph(x.size(), &x.at(0), &y.at(0));
        g->SetTitle(name);
        g->SetName(name);
        return g;
    }

    TGraph *  createUncGraphExp(TString name)const{
        TGraph * g = new TGraph(exp_unc.size(), &lumi.at(0), &exp_unc.at(0));
        g->SetTitle(name);
        g->SetName(name);
        g->SetLineWidth(2);
        g->SetLineColor(kRed);
        g->GetXaxis()->SetTitle("#Delta(lumi) [%]");
        g->GetYaxis()->SetTitle("#Delta [%]");
        return g;
    }
    TGraph *  createUncGraphPDF(TString name)const{
        TGraph * g = new TGraph(pdf_unc.size(), &lumi.at(0), &pdf_unc.at(0));
        g->SetTitle(name);
        g->SetName(name);
        g->SetLineWidth(2);
        g->SetLineColor(kBlue);
        return g;
    }
    TGraph *  createUncGraphScale(TString name)const{
        TGraph * g = new TGraph(scale_unc.size(), &lumi.at(0), &scale_unc.at(0));
        g->SetTitle(name);
        g->SetName(name);
        g->SetLineWidth(2);
        g->SetLineColor(kViolet);
        return g;
    }

private:
    //all unc relative in percent
    std::vector<double> lumi, as, exp_unc, pdf_unc, scale_unc;

};

void drawSubLegend(TLegend* leg, double x0, double y0, double x1, double x2,
        int linestyleas, int linestylemt){

    //create dummy for second legend
    if(!leg){
        leg=new TLegend(x0,y0,x1,x2);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);}
    else{
        auto e=leg->AddEntry("","","");
        e->SetLabel("");
        e->SetFillStyle(0);
        e->SetTextColorAlpha(0,0);
        e->SetFillColor(0);
        e->SetLineColorAlpha(0,0);
        e->SetLineColor(0);
        e->SetLineWidth(0);
        e->SetDrawOption("p");
        e->SetMarkerSize(0);
        e->SetMarkerColor(0);
        e->SetMarkerColorAlpha(0,0);
    }
    TGraph * g = new TGraph();
    g->SetLineWidth(2);
    g->SetLineStyle(linestyleas);

    leg->AddEntry(g,"#alpha_{S}","l");
    g = new TGraph();
    g->SetLineWidth(2);
    g->SetLineStyle(linestylemt);
    leg->AddEntry(g,"m_{t}","l");

    leg->Draw("same");

}

static void setDefaultAxis(TH1F& axis, double multiplier=1.){
    double padwidth=1/multiplier;

    double labelsize=0.05;
    //axis.GetYaxis()->SetTickLength(0);
    axis.GetYaxis()->SetLabelSize(labelsize*multiplier);
    axis.GetXaxis()->SetLabelSize(labelsize*multiplier);
    axis.GetXaxis()->SetTitleSize(0.05*multiplier);
    axis.GetYaxis()->SetTitleSize(0.05*multiplier);
    //axis.GetXaxis()->SetNdivisions(505);
}

void drawLogos(bool prelim=true){

    TLatex *tex = new TLatex(0.185,0.82,"CMS");
    tex->SetNDC(true);
    tex->SetTextFont(61);
    tex->SetTextSize(0.08);
    tex->SetLineWidth(2);
    tex->Draw();

    TString add="#it{Projection}";
    if(prelim)
        add+=" #it{Preliminary}";
    tex = new TLatex(0.185,0.76,add);
    tex->SetNDC(true);
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->Draw();

    TString lumi = "3000 fb^{-1} (13 TeV)";
    tex = new TLatex(0.9,0.91,lumi);
    tex->SetNDC(true);
    tex->SetTextAlign(31);
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(1);
    tex->SetNDC(true);
    tex->Draw();

    gPad->SetTickx();
    gPad->SetTicky();

}

double sumsq(double a, double b){
    return std::sqrt(a*a+b*b);
}

TGraph*  scanMeasUnc(
        const double meas_unc,
        TString title){

    std::vector<double> lumi_unc;
    const double measured_xsec=globals::measured_xsec;
    for(double i=0.5;i<3.01;i+=0.1)
            lumi_unc.push_back(i);

    std::vector<double>  out;
    for(const auto lumi: lumi_unc){
        double totalunc = sumsq(meas_unc, lumi);
        out.push_back(totalunc);
    }
    auto g=  new TGraph(out.size(),&lumi_unc.at(0),&out.at(0));

    g->SetTitle(title);
    return g;

}

TGraph*  performExtraction(const std::string predfile,
        const double meas_unc,
        const double pdf_scaler,
        const double scale_scaler,
        TString title
){

    std::vector<double> lumi_unc;

    for(double i=0.5;i<3.01;i+=0.1)
        lumi_unc.push_back(i);


    std::vector<double>  out;

    for(const auto lumi: lumi_unc){
        const double measured_xsec=globals::measured_xsec;
        globals::measured_errdown=sumsq(meas_unc*0.01*measured_xsec, lumi*0.01*measured_xsec);
        globals::measured_errup=sumsq(meas_unc*0.01*measured_xsec, lumi*0.01*measured_xsec);

       // std::cout << "globals::measured_errup " << globals::measured_errup << std::endl;

        extractor tmp_ex;
        tmp_ex.scalePDFUnc(pdf_scaler);
        tmp_ex.readPrediction(predfile);
        extractResult resmtop=tmp_ex.extractMtop();


        tmp_ex.fixScaleVar(extractor::prior_fixed0);
        //  std::cout << "scale 0: "<< resmtop.central << " + " << resmtop.errup << " - " << resmtop.errdown << std::endl;
        tmp_ex.fixScaleVar(extractor::prior_fixedm1);
        auto resmtop2=tmp_ex.extractMtop();
        //   std::cout << "scale -1: "<< resmtop2.central << " + " << resmtop2.errup << " - " << resmtop2.errdown << std::endl;
        tmp_ex.fixScaleVar(extractor::prior_fixed1);
        auto resmtop3=tmp_ex.extractMtop();
        //   std::cout << "scale 1: "<< resmtop3.central << " + " << resmtop3.errup << " - " << resmtop3.errdown << std::endl;

        double scale_low =  scale_scaler*fabs(resmtop2.central-resmtop.central);
        double scale_hi =  scale_scaler*fabs(resmtop3.central-resmtop.central);
        double relup = sumsq(resmtop.errup, scale_hi)/resmtop.central*100.;
        double reldown=sumsq(resmtop.errdown, scale_low)/resmtop.central*100.;


        double max=relup;
        if(max<reldown)
            max=reldown;

        out.push_back(max);
    }

    auto g=  new TGraph(out.size(),&lumi_unc.at(0),&out.at(0));

    g->SetTitle(title);
    return g;
}


TGraph * drawDoubleLine(TGraph * g,
        int lineStyle,
        int lineColor){

    TGraph *g2 = (TGraph *)g->Clone();
    g2->SetLineColor(kBlack);
    g2->SetLineWidth(4);
    g2->SetLineStyle(lineStyle);
    g2->Draw("L,same");

    g->SetLineColor(lineColor);
    g->SetLineWidth(3);
    g->SetLineStyle(lineStyle);
    g->Draw("L,same");
    return g;
}

void drawMarker(double ymin,double ymax){
    TLine * l = new TLine(1.,ymin,1,ymax);//ymin+(ymax-ymin)/10.);
    l->SetLineColorAlpha(kBlack,0.05);
    l->SetLineWidth(10);
   // l->SetLineStyle(9);
    l->Draw("same");

}

TH1F makeAxisHist(TGraph* g, double xmin=-1, double xmax=-1){
    if(xmin==xmax)
        return *g->GetHistogram();

    else return TH1F("","",100,xmin,xmax);

}

TCanvas* createCanvas(){
    TCanvas * cv=new TCanvas();
    cv->SetBottomMargin(0.15);
    cv->SetLeftMargin(0.15);
    return cv;
}

int main(int argc, char* argv[]){

    bool do_half_unc=false;

    double reduced=0.01;
    if(do_half_unc){
        reduced=0.5;
    }

    //some graphs to save for BRIL plot
    gStyle->SetOptStat(0);

    TGraph * mxsec_run2, * mxsec_yr,
    * mt_full , * mt_noscale, * mt_nopdf, * mt_nopdf_noscale,
    * as_full , * as_noscale, * as_nopdf, * as_nopdf_noscale;



    globals::init();

    std::vector<int> colours;
    TColor col;


	//fileReader::debug=true;
    std::string predfile="data/ABMP16_13TeV.data";
	if(argc>=2){
	    predfile = argv[1];
	}



	std::cout << "PDF " << predfile <<std::endl;
	extractor tmp_ex;
	tmp_ex.readPrediction(predfile);


	std::vector<double> lumi_unc;

	//this is without lumi. Put numbers here
	std::vector<TGraph * > outgraphs;
	std::vector<double> meas_uncs = {3.1, 1.5}; //percent


	std::vector<TString> names = {
	        "#Delta(#sigma_{t#bar{t}}) = 3.1%",
	        "#Delta(#sigma_{t#bar{t}}) = 1.5%",
	};
    colours={col.GetColor("#5e3c99"),col.GetColor("#e66101")};

    std::vector<int> line_styles={1,1,1,1,1,1,1,1,1}; //1,2,9

	double min_unc=1e6, max_unc=-1;

	for(size_t im=0;im<meas_uncs.size();im++){

	    const auto meas_unc = meas_uncs.at(im);

	    auto g = performExtraction(predfile,
	            meas_unc,
	            1.,
	            1.,
	            names.at(im));

	    outgraphs.push_back(g);
	}

	mxsec_run2 = outgraphs.at(0);
    mxsec_yr = outgraphs.at(1);


	TCanvas * cv = createCanvas();
	TH1F axis= *outgraphs.at(0)->GetHistogram();
    setDefaultAxis(axis);
    axis.Draw("AXIS");
    axis.GetYaxis()->SetRangeUser(1.1, 1.5);
    drawMarker(1.1, 1.5);

    axis.GetYaxis()->SetTitle("#Delta(m_{t}) [%]");
    axis.GetXaxis()->SetTitle("#Delta(luminosity) [%]");

    TLegend * leg=new TLegend(0.57,0.22,0.85,0.44);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

	for(size_t i=0;i<outgraphs.size();i++){

	    auto g = outgraphs.at(i);
	    drawDoubleLine(g, line_styles.at(i),colours.at(i));
        leg->AddEntry(g,g->GetTitle(),"lp");
	}

	leg->Draw("same");
	drawLogos();
	cv->Print("meas_lumi.pdf");


	outgraphs.clear();
	names = {
	            "#Delta(scale)/2",
	            "#Delta(PDF)/2",
                "#Delta(scale)/2, #Delta(PDF)/2",
	    };
    std::vector<std::pair<double,double> > scale_pdf_multis = {{0.5,1.}, {1.,0.5}, {0.5, 0.5}};

	colours={col.GetColor("#5e3c99"),col.GetColor("#e66101"),col.GetColor("#fdb863"),col.GetColor("#b2abd2")};


	for(size_t im=0;im<scale_pdf_multis.size();im++){

	    const auto multi = scale_pdf_multis.at(im);

	    std::cout << multi.first << ", PDF: " << multi.second << std::endl;

	    auto g = performExtraction(predfile,
	            meas_uncs.at(1),
	            multi.second,
                multi.first,
	            names.at(im));

	    outgraphs.push_back(g);
	}


    mt_noscale=outgraphs.at(0);
    mt_nopdf=outgraphs.at(1);
    mt_nopdf_noscale=outgraphs.at(2);


	cv = createCanvas();
    axis= *outgraphs.at(0)->GetHistogram();
    setDefaultAxis(axis);
    axis.Draw("AXIS");
    axis.GetYaxis()->SetRangeUser(0.65, 1.4);
    drawMarker(0.65, 1.4);

    axis.GetYaxis()->SetTitle("#Delta(m_{t}) [%]");
    axis.GetXaxis()->SetTitle("#Delta(luminosity) [%]");

    leg=new TLegend(0.50,0.18,0.88,0.45);
    leg->SetBorderSize(0);

    for(size_t i=0;i<outgraphs.size();i++){
        auto g = outgraphs.at(i);
        drawDoubleLine(g, line_styles.at(i),colours.at(i));
        leg->AddEntry(g,g->GetTitle(),"lp");
    }

    leg->Draw("same");
    drawLogos();
    cv->Print("theo_lumi.pdf");

    names = {
                "#Delta(#sigma_{t#bar{t}}) = 3.1%",
                "#Delta(#sigma_{t#bar{t}}) = 1.5%",
        };
    outgraphs.clear();

    for(size_t im=0;im<meas_uncs.size();im++){

        const auto meas_unc = meas_uncs.at(im);
        auto g = scanMeasUnc(meas_unc, names.at(im));
        outgraphs.push_back(g);
    }
    cv = createCanvas();
    axis= *outgraphs.at(0)->GetHistogram();
    setDefaultAxis(axis);
    axis.Draw("AXIS");
    axis.GetYaxis()->SetRangeUser(1.0, 5.0);
    drawMarker(1.0,5.0);

    axis.GetYaxis()->SetTitle("#Delta_{tot}(#sigma_{t#bar{t}}) [%]");
    axis.GetXaxis()->SetTitle("#Delta(luminosity) [%]");

    leg=new TLegend(0.60,0.2,0.87,0.35);
    leg->SetBorderSize(0);

    for(size_t i=0;i<outgraphs.size();i++){
        auto g = outgraphs.at(i);
        drawDoubleLine(g, line_styles.at(i),colours.at(2*i));
        leg->AddEntry(g,g->GetTitle(),"lp");
    }

    leg->Draw("same");
    drawLogos();
    cv->Print("meas_lumi.pdf");

    //now the alpa_s data from Katerina

    // * as_noscale, as_nlowpdf, as_lowpdf_lowscale

    as_extraction as_ex;

    //just check
    cv = createCanvas();
    auto guexp = as_ex.createUncGraphExp("exp");
    auto gupdf = as_ex.createUncGraphPDF("pdf");
    auto guscale = as_ex.createUncGraphScale("scale");

    leg=new TLegend(0.60,0.7,0.87,0.8);
    leg->AddEntry(guexp   ,guexp->GetName(),"l");
    leg->AddEntry(gupdf   ,gupdf->GetName(),"l");
    leg->AddEntry(guscale ,guscale->GetName(),"l");


    guexp  ->Draw("APl");
    guexp->GetYaxis()->SetRangeUser(0.,2.5);
    drawMarker(0,2.5);
    gupdf  ->Draw("same,l");
    guscale->Draw("same,l");

    leg->Draw("same");

    cv->Print("as_unc.pdf");

    cv = createCanvas();
    double mtmulti=1.;

    as_full = as_ex.createGraph("full",1,1.,1/mtmulti);
    as_noscale = as_ex.createGraph("lowscale",1,reduced,1/mtmulti);
    as_nopdf = as_ex.createGraph("lowpdf",reduced,1.,1/mtmulti);
    as_nopdf_noscale = as_ex.createGraph("lowpdfscale",reduced,reduced,1/mtmulti);







    double legendleft = 0.67;
    double legendright = 0.89;
    double legendlow = 0.3;
    double legendheight=0.5;

    axis= makeAxisHist(0, 0.45, 3.95); //*as_noscale->GetHistogram();
    setDefaultAxis(axis);
    axis.Draw("AXIS");

    double ymin=0.65, ymax=2.3;



    axis.GetYaxis()->SetRangeUser(ymin, ymax);
    drawMarker(ymin,ymax);

    axis.GetYaxis()->SetTitle("#Delta(m_{t}), #Delta(#alpha_{S}) [%]");
    axis.GetXaxis()->SetTitle("#Delta(luminosity) [%]");


    leg=new TLegend(legendleft,legendlow,legendright,legendlow+legendheight);
    leg->SetNColumns(1);
    //leg->SetBorderSize(0);
   // leg->SetFillStyle(0);

    drawDoubleLine(as_full, 1,colours.at(0));

    drawDoubleLine(as_nopdf, 1,colours.at(2));

    drawDoubleLine(as_noscale, 1,colours.at(1));


    drawDoubleLine(as_nopdf_noscale, 1,colours.at(3));

    if(do_half_unc){
        leg->AddEntry(as_full,"Full #Delta(theo)","l");
        leg->AddEntry(as_nopdf,"#Delta(PDF)/2","l");
        leg->AddEntry(as_noscale,"#Delta(scale)/2","l");
        leg->AddEntry(as_nopdf_noscale,"#Delta(theo)/2","l");
    }
    else{
       leg->AddEntry(as_full,"Full #Delta(theo)","l");
       leg->AddEntry(as_nopdf,"No #Delta(PDF)","l");
       leg->AddEntry(as_noscale,"No #Delta(scale)","l");
       leg->AddEntry(as_nopdf_noscale,"No #Delta(theo)","l");

    }
/*
 *
    * mt_full , * mt_noscale, * mt_nopdf, * mt_nopdf_noscale,
    * as_full , * as_noscale, * as_nopdf, * as_nopdf_noscale;
 */


    mt_full =          performExtraction(predfile,meas_uncs.at(1),1.,1., "mt_full");
    mt_noscale =       performExtraction(predfile,meas_uncs.at(1),1.,reduced, "mt_noscale");
    mt_nopdf =         performExtraction(predfile,meas_uncs.at(1),reduced,1., "mt_nopdf");
    mt_nopdf_noscale = performExtraction(predfile,meas_uncs.at(1),reduced,reduced, "mt_nopdf_noscale");


    drawDoubleLine(mt_full,        2,colours.at(0));
    drawDoubleLine(mt_noscale,        2,colours.at(1));
    drawDoubleLine(mt_nopdf,          2,colours.at(2));
    drawDoubleLine(mt_nopdf_noscale, 2,colours.at(3));

    leg->Draw("same");

    drawSubLegend(leg,legendleft,0.2,legendright,0.3,
            1,2);


    drawLogos();
    cv->Print("full_theo_lumi.pdf");


    cv = createCanvas();
    axis= makeAxisHist(0, 0.45, 3.95); //*as_noscale->GetHistogram();
    setDefaultAxis(axis);
    axis.Draw("AXIS");

    axis.GetYaxis()->SetRangeUser(ymin, ymax);
    drawMarker(ymin,ymax);

    axis.GetYaxis()->SetTitle("#Delta(m_{t}), #Delta(#alpha_{S}) [%]");
    axis.GetXaxis()->SetTitle("#Delta(luminosity) [%]");

    leg=new TLegend(legendleft,legendlow+0.05,legendright,legendlow+legendheight-0.05);
    leg->SetNColumns(1);

    drawDoubleLine(as_full, 1,colours.at(0));
    drawDoubleLine(as_nopdf_noscale, 1,colours.at(2));


    leg->AddEntry(as_full,"Full #Delta(theo)","l");
    leg->AddEntry(as_nopdf_noscale,"No #Delta(theo)","l");



    drawDoubleLine(mt_full,        2,colours.at(0));
    drawDoubleLine(mt_nopdf_noscale, 2,colours.at(2));

    leg->Draw("same");

    drawSubLegend(leg,legendleft,0.2,legendright,0.3,
            1,2);

    drawLogos();
    cv->Print("simple_theo_lumi.pdf");

	return 0;

}







