#ifndef BMNMWPCTRACKTODC_H
#define BMNMWPCTRACKTODC_H

#include "BmnTrack.h"

class BmnMwpcTrackToDC : public BmnTrack {
public:
  BmnMwpcTrackToDC();
  virtual ~BmnMwpcTrackToDC();
private:
  ClassDef(BmnMwpcTrackToDC, 1);
};

#endif
