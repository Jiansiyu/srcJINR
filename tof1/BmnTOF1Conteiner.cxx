/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "BmnTOF1Conteiner.h"

//--------------------------------------------------------------------------------------------------------------------------------------
BmnTOF1Conteiner::BmnTOF1Conteiner()
{
	//Clear();
}
//--------------------------------------------------------------------------------------------------------------------------------------
BmnTOF1Conteiner::BmnTOF1Conteiner( int plane, int strip, double time, double amp, double x_local, double y_local, double z_local, double x_glob, double y_glob, double z_glob)
{
	fPlane = plane;
	fStrip = strip;
	// fTimeL
	// fTimeR
	// fAmpL
	// fAmpR

	fTime = time;
	fAmp = amp;

	fX_local = x_local;
	fY_local = y_local;
	fZ_local = z_local;
	fX_glob = x_glob;
	fY_glob = y_glob;
	fZ_glob = z_glob;
	
	fdL = sqrt( pow(fX_glob,2) + pow(fY_glob,2) + pow(fZ_glob,2) );
}

void BmnTOF1Conteiner::SetParameters ( int plane, int strip, double time, double amp, double x_local, double y_local, double z_local, double x_glob, double y_glob, double z_glob)
{
	fPlane = plane;
	fStrip = strip;
	// fTimeL
	// fTimeR
	// fAmpL
	// fAmpR

	fTime = time;
	fAmp = amp;

	fX_local = x_local;
	fY_local = y_local;
	fZ_local = z_local;
	fX_glob = x_glob;
	fY_glob = y_glob;
	fZ_glob = z_glob;

	fdL = sqrt( pow(fX_glob,2) + pow(fY_glob,2) + pow(fZ_glob,2) );
}
//--------------------------------------------------------------------------------------------------------------------------------------
void	BmnTOF1Conteiner::print( std::ostream& os, const char* comment)const
{
	os<<" [BmnTof1Digit] "; if(nullptr != comment) os<<comment;
	os<<"  detID: "<<fPlane<<", stripID: "<<fStrip<<", Time: "<<fTime<<", Width: "<<fAmp<<std::endl;
}
//--------------------------------------------------------------------------------------------------------------------------------------
void	BmnTOF1Conteiner::Clear()
{
	fPlane = -1;
	fStrip = -1;
	// fTimeL
	// fTimeR
	// fAmpL
	// fAmpR

	fTime = -1;
	fAmp = -1;

	fX_local = -1000;
	fY_local = -1000;
	fZ_local = -1000;
	fX_glob = -1000;
	fY_glob = -1000;
	fZ_glob = -1000;
}
//
ClassImp(BmnTOF1Conteiner)

