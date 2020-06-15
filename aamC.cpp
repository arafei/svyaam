// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]]

Rcpp::List aamC(arma::dvec & time, arma::dvec & CR, arma::dmat & DM, arma::dvec & WT, arma::dvec & ST, arma::dvec & CL, std::string & MT){

	unsigned int n = time.n_elem;
	unsigned int r = DM.n_cols-1;

	arma::dcube Yt(n, r+1, n, arma::fill::zeros);
	arma::dcube gYt(r+1, r+1, n, arma::fill::zeros);
	arma::dcube Xt(r+1, n, n, arma::fill::zeros);
	arma::dmat Ast(r+1, n, arma::fill::zeros);
	arma::uvec tmp = (time>=time(0));
	Yt.slice(0) = DM % repmat(arma::conv_to<arma::dvec>::from(tmp), 1, r+1);
	//Yt.slice(0).rows( find(time>=time(0)) ) = DM.rows( find(time>=time(0)) );
	gYt.slice(0) = pinv(trans(Yt.slice(0).each_col() % WT) * Yt.slice(0));
	Xt.slice(0) = gYt.slice(0) * trans(Yt.slice(0).each_col() % WT);
	arma::dvec Ik = arma::conv_to<arma::dvec>::from(time==time(0) && CR==1);
	Ast.col(0) = Xt.slice(0) * Ik;
	arma::dcube Ost(r+1, r+1, n, arma::fill::zeros);
	for(unsigned int t=1; t<n; t++){
		tmp = (time>=time(t));
		Yt.slice(t) = DM % repmat(arma::conv_to<arma::vec>::from(tmp), 1, r+1);
		//Yt.slice(t).rows( find(time>=time(t)) ) = DM.rows( find(time>=time(t)) );
		gYt.slice(t) = pinv(trans(Yt.slice(t).each_col() % WT) * Yt.slice(t));
		Xt.slice(t) = gYt.slice(t) * trans(Yt.slice(t).each_col() % WT);
		Ik = arma::conv_to<arma::vec>::from(time==time(t) && CR==1);
		Ast.col(t) = Ast.col(t-1) + Xt.slice(t) * Ik;
	}
	if(MT=="no"){
		Ik = arma::conv_to<arma::vec>::from(time==time(0) && CR==1);
		Ost.slice(0) = (Xt.slice(0).each_row() % Ik.t()) * trans(Xt.slice(0));
		for(unsigned int i=1; i<n; i++){
			Ik = arma::conv_to<arma::vec>::from(time==time(i) && CR==1);
			Ost.slice(i) = Ost.slice(i-1) + (Xt.slice(i).each_row() % Ik.t()) * trans(Xt.slice(i));
		}
	}else if(MT=="TSL"){
		arma::dcube TDit(n, r+1, n, arma::fill::zeros);
		for(unsigned int t=0; t<n; t++){
			for(unsigned int j=0; j<n; j++){
				TDit.slice(t).row(j) = Yt.slice(t).row(j) * gYt.slice(t); 
 	     		}
		}
		arma::dcube TD(n, r+1, n, arma::fill::zeros);
		arma::dcube wTD(n, r+1, n, arma::fill::zeros);
		Ik = arma::conv_to<arma::vec>::from(time==time(0) && CR==1);
		TD.slice(0) = TDit.slice(0) % repmat(Ik - TDit.slice(0) * trans(Yt.slice(0).each_col() % WT) * Ik, 1, r+1);
		wTD.slice(0) = repmat(WT, 1, r+1) % TD.slice(0);

		arma::uword H = max(ST);
		arma::uword L = 0;
		for(unsigned int h=1; h<=H; h++){
			L = max( CL.elem( find(ST==h) ) );
			arma::dmat wtTDpsu(L, r+1, arma::fill::zeros);
			for(unsigned int l=1; l<=L; l++){
				wtTDpsu.row(l-1) = sum( wTD.slice(0).rows( find(ST==h && CL==l) ), 0);
			}
			arma::dmat meanpsui = mean(wtTDpsu, 0);
			arma::dmat tmpcomp = wtTDpsu.each_row() - meanpsui;
			Ost.slice(0) = Ost.slice(0) + trans(tmpcomp) * tmpcomp * L/(L-1.0);
		}

		for(unsigned int i=1; i<n; i++){
			Ik = arma::conv_to<arma::vec>::from(time==time(i) && CR==1);
			TD.slice(i) = TD.slice(i-1) + TDit.slice(i) % repmat(Ik - TDit.slice(i) * trans(Yt.slice(i).each_col() % WT) * Ik, 1, r+1);
			wTD.slice(i) = repmat(WT, 1, r+1) % TD.slice(i);

			//arma::dmat cs(r+1, r+1, arma::fill::zeros);
			for(unsigned int h=1; h<=H; h++){
				L = max( CL.elem( find(ST==h) ) );
				arma::dmat wtTDpsu(L, r+1, arma::fill::zeros);
				for(unsigned int l=1; l<=L; l++){
					wtTDpsu.row(l-1) = sum( wTD.slice(i).rows( find(ST==h && CL==l) ), 0);
				}
				arma::dmat meanpsui = mean(wtTDpsu, 0);
				arma::dmat tmpcomp = wtTDpsu.each_row() - meanpsui;
				Ost.slice(i) = Ost.slice(i) + trans(tmpcomp) * tmpcomp * L/(L-1.0);
			}
		}
		TDit.reset();
		TD.reset();
		wTD.reset();
		//wTDpsu.reset();
		//tmpcomp.reset();
	}else if(MT=="JK1"){
		arma::dvec repwght(n, arma::fill::zeros);
		arma::dcube gYthl(r+1, r+1, n, arma::fill::zeros);
		arma::dcube Xthl(r+1, n, n, arma::fill::zeros);
		arma::dcube Osthl(r+1, r+1, n, arma::fill::zeros);
		arma::dmat Asthl(r+1, n, arma::fill::zeros);
		arma::uword H = max(ST);
		arma::uword L = 0;
		for(unsigned int h=1; h<=H; h++){
			//generating replication weights;
			L = max( CL.elem( find(ST==h) ) );
			for(unsigned int l=1; l<=L; l++){
				repwght = WT;
				repwght.elem( find(ST==h && CL==l) ).zeros();
				repwght.elem( find(ST==h && CL!=l) ) = repwght.elem( find(ST==h && CL!=l) ) * L/(L-1.0);
				// Construct Xt_(hl)(t);
				Asthl.zeros();
				gYthl.slice(0) = pinv(trans(Yt.slice(0).each_col() % repwght) * Yt.slice(0));
				Xthl.slice(0) = gYthl.slice(0) * trans(Yt.slice(0).each_col() % repwght);
				Ik = arma::conv_to<arma::vec>::from(time==time(0) && CR==1);
				Asthl.col(0) = Xthl.slice(0) * Ik;
				Osthl.slice(0) = (Asthl.col(0) - Ast.col(0)) * trans(Asthl.col(0) - Ast.col(0)) * (L-1.0) / L;
				for(unsigned int t=1; t<n; t++){
					gYthl.slice(t) = pinv(trans(Yt.slice(t).each_col() % repwght) * Yt.slice(t));
					Xthl.slice(t) = gYthl.slice(t) * trans(Yt.slice(t).each_col() % repwght);
					Ik = arma::conv_to<arma::vec>::from(time==time(t) && CR==1);
					Asthl.col(t) = Asthl.col(t-1) + Xthl.slice(t) * Ik;
					Osthl.slice(t) = (Asthl.col(t) - Ast.col(t)) * trans(Asthl.col(t) - Ast.col(t)) * (L-1.0) / L;
				}
				Ost = Ost + Osthl;
			}
		}
		repwght.reset();
		gYthl.reset();
		Xthl.reset();
		Asthl.reset();
		Osthl.reset();
	}
	return(List::create(Named("AstX", Ast), Named("OstX", Ost)));
}


//aamC <- cxxfunction(signature(meth="character", tm = "numeric", cen = "numeric", mx = "numeric", wght = "numeric", clr = "numeric", str = "numeric"), includes=inc, src, plugin="RcppArmadillo")

