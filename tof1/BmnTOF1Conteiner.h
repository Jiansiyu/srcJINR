/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BmnTOF1Conteiner.h
 * Author: mikhailr
 *
 * Created on May 3, 2018, 3:46 PM
 */

#ifndef BMNTOF1CONTEINER_H
#define BMNTOF1CONTEINER_H 1

#include "TObject.h"
#include <iostream>
#include <cmath>
//------------------------------------------------------------------------------------------------------------------------
class BmnTOF1Conteiner : public TObject {

	private: 
		int fPlane;
		int fStrip;
		// fTimeL
		// fTimeR
		// fAmpL
		// fAmpR

		double fTime;
		double fAmp;

		double fX_local;
		double fY_local;
		double fZ_local;
		double fX_glob;
		double fY_glob;
		double fZ_glob;
		double fdL;

	public:
		BmnTOF1Conteiner();
		BmnTOF1Conteiner( int plane, int strip, double time, double amp, double x_local, double y_local, double z_local, double x_glob, double y_glob, double z_glob);
		virtual ~BmnTOF1Conteiner() {};

		int GetStrip()      	const { return fStrip; }
		int GetPlane()		const { return fPlane; }

		double GetTime()       	const { return fTime; }
		double GetAmp()      	const { return fAmp; }
		double GetXLocal()       	const { return fX_local; }
		double GetYLocal()       	const { return fY_local; }
		double GetZLocal()       	const { return fZ_local; }
		double GetXGlobal()       	const { return fX_glob; }
		double GetYGlobal()       	const { return fY_glob; }
		double GetZGlobal()       	const { return fZ_glob; }
		double GetdL()			const { return fdL; }

		void SetParameters ( int plane, int strip, double time, double amp, double x_local, double y_local, double z_local, double x_glob, double y_glob, double z_glob);
		
		void Clear();

		void	print(std::ostream& os = std::cout, const char* comment = nullptr)const;

		ClassDef(BmnTOF1Conteiner, 2);
};


#endif /* BMNTOF1CONTEINER_H */

