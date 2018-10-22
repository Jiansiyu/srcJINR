/*
 * FairMCPointsDraw.cpp
 *
 *  Created on: Apr 17, 2009
 *      Author: stockman
 */

#include "MpdHitPointSetDraw.h"
#include "MpdEventManager.h"

#include "FairHit.h"

#include "TEveManager.h"

#include <iostream>
using namespace std;


TVector3 MpdHitPointSetDraw::GetVector(TObject* obj)
{
    FairHit* hit = (FairHit*) obj;
    cout<<"MpdHitPointSetDraw::GetVector(): "<<hit->GetX()<<" "<<hit->GetY()<<" "<<hit->GetZ()<<endl;

	// SOMEHOW THE COORDINATES GET MESSED UP??

	cout << "\t " << hit->GetZ() << " " << hit->GetX() << " " << hit->GetY() << "\n";

    return TVector3(hit->GetZ(), hit->GetX(), hit->GetY());
}

void MpdHitPointSetDraw::AddEveElementList()
{
    fEventManager->AddEventElement(fq, RecoPointList);
    return;
}

void MpdHitPointSetDraw::RemoveEveElementList()
{
    gEve->RemoveElement(fq, fEventManager->EveRecoPoints);
    return;
}

ClassImp(MpdHitPointSetDraw)
