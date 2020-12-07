/*
 * globals.h
 *
 *  Created on: 14 Sep 2017
 *      Author: jkiesele
 */

#ifndef TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_GLOBALS_H_
#define TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_GLOBALS_H_


namespace globals{

extern double measured_xsec;
extern double measured_errdown;
extern double measured_errup;

const double default_mtop=172.;
const double default_alphas=0.118;
const double default_mtop_err=1;
const double default_alphas_err=0.0012;

const double default_mtop_low=165;
const double default_mtop_high=185;
const double default_alphas_low=0.109;
const double default_alphas_high=0.124;


void init();
}


#endif /* TOPLHCWG_CROSSSECTION_SUMMER16_COMBINATION_MASS_ALPHA_S_INCLUDE_GLOBALS_H_ */
