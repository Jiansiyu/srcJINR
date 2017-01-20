#include <Rtypes.h>
#include "BmnDchTrackFinder.h"

BmnDchTrackFinder::BmnDchTrackFinder() :
fSegmentMatching(kFALSE),
has7DC1(kFALSE),
has7DC2(kFALSE),
nDC1_segments(0),
nDC2_segments(0),
pair_x2(0),
pair_y2(0),
pair_u2(0),
pair_v2(0),
single_xa2(0),
single_ya2(0),
single_ua2(0),
single_va2(0),
single_xb2(0),
single_yb2(0),
single_ub2(0),
single_vb2(0),
pair_x1(0),
pair_y1(0),
pair_u1(0),
pair_v1(0),
single_xa1(0),
single_ya1(0),
single_ua1(0),
single_va1(0),
single_xb1(0),
single_yb1(0),
single_ub1(0),
single_vb1(0) {
    fEventNo = 0;
    tracksDch = "BmnDchTrack";
    InputDigitsBranchName = "DCH";

    prev_wire = -1;
    prev_time = -1;

    Z_dch1 = 529.65;
    Z_dch2 = 727.65;
    Z_dch_mid = 628.65;

    x1_sh = 15.0;
    x2_sh = 12.2;
    y1_sh = 0.2;
    y2_sh = -0.3;

    x1_slope_sh = -0.08;
    y1_slope_sh = 0.051;
    x2_slope_sh = -0.08;
    y2_slope_sh = 0.055;

    scale = 0.5;
}

BmnDchTrackFinder::~BmnDchTrackFinder() {

}

void BmnDchTrackFinder::Exec(Option_t* opt) {
    PrepareArraysToProcessEvent();
    cout << "\n======================== DCH track finder exec started =====================\n" << endl;
    cout << "Event number: " << fEventNo++ << endl;

    //temporary containers
    // Order used: va1, vb1, ua1, ub1, ya1, yb1, xa1, xb1 (dch1, 0 - 7) - va2, vb2, ua2, ub2, ya2, yb2, xa2, xb2 (dch2, 8 - 15)
    const Int_t nDim = 20;
    const Int_t nPlanes = 16; // Total number of planes in both DCHs (0-7, 8-15)

    Double_t times[nPlanes][nDim];
    Double_t wires[nPlanes][nDim];
    Int_t it[nPlanes];
    Bool_t used[nPlanes][nDim];
    for (Int_t iPlanes = 0; iPlanes < nPlanes; iPlanes++) {
        it[iPlanes] = 0;
        for (Int_t iDim = 0; iDim < nDim; iDim++) {
            used[iPlanes][iDim] = kFALSE;
            times[iPlanes][iDim] = 0.;
            wires[iPlanes][iDim] = 0.;
        }
    }

    Bool_t goodEv = kTRUE;
    Bool_t written = kFALSE;

    for (Int_t iDig = 0; iDig < fBmnDchDigitsArray->GetEntriesFast(); ++iDig) {
        BmnDchDigit* digit = (BmnDchDigit*) fBmnDchDigitsArray->UncheckedAt(iDig);

        //skip identical events
        if (!written) {
            written = kTRUE;
            if (digit->GetTime() == prev_time && digit->GetWireNumber() == prev_wire) {
                goodEv = kFALSE;
            } else {
                prev_time = Int_t(digit->GetTime());
                prev_wire = Int_t(digit->GetWireNumber());
            }
        }//!written

        if (!goodEv)
            return;

        // Order used: va1(0), vb1(1), ua1(2), ub1(3), ya1(4), yb1(5), xa1(6), xb1(7) -> 
        //             va2(8), vb2(9), ua2(10), ub2(11), ya2(12), yb2(13), xa2(14), xb2(15)
        Short_t plane = digit->GetPlane();
        Short_t wire = digit->GetWireNumber();
        Double_t time = digit->GetTime();
        Bool_t secondaries = kFALSE;

        for (Int_t sec = 0; sec < it[plane] - 1; sec++)
            if (wire == wires[plane][sec]) {
                secondaries = kTRUE;
                break;
            }

        if (it[plane] == (nDim - 1) || secondaries)
            continue;

        wires[plane][it[plane]] = wire;
        times[plane][it[plane]] = time;
        it[plane]++;
    }

    pair_v1 = Reconstruction(1, "v", pair_v1, it[0], it[1], wires[0], wires[1], times[0], times[1], used[0], used[1], v1_ab, sigm_v1_ab);
    pair_u1 = Reconstruction(1, "u", pair_u1, it[2], it[3], wires[2], wires[3], times[2], times[3], used[2], used[3], u1_ab, sigm_u1_ab);
    pair_y1 = Reconstruction(1, "y", pair_y1, it[4], it[5], wires[4], wires[5], times[4], times[5], used[4], used[5], y1_ab, sigm_y1_ab);
    pair_x1 = Reconstruction(1, "x", pair_x1, it[6], it[7], wires[6], wires[7], times[6], times[7], used[6], used[7], x1_ab, sigm_x1_ab);

    single_va1 = ReconstructionSingle(1, "v", "a", single_va1, it[0], wires[0], times[0], used[0], v1_single, sigm_v1_single);
    single_vb1 = ReconstructionSingle(1, "v", "b", single_vb1, it[1], wires[1], times[1], used[1], v1_single, sigm_v1_single);
    single_ua1 = ReconstructionSingle(1, "u", "a", single_ua1, it[2], wires[2], times[2], used[2], u1_single, sigm_u1_single);
    single_ub1 = ReconstructionSingle(1, "u", "b", single_ub1, it[3], wires[3], times[3], used[3], u1_single, sigm_u1_single);
    single_ya1 = ReconstructionSingle(1, "y", "a", single_ya1, it[4], wires[4], times[4], used[4], y1_single, sigm_y1_single);
    single_yb1 = ReconstructionSingle(1, "y", "b", single_yb1, it[5], wires[5], times[5], used[5], y1_single, sigm_y1_single);
    single_xa1 = ReconstructionSingle(1, "x", "a", single_xa1, it[6], wires[6], times[6], used[6], x1_single, sigm_x1_single);
    single_xb1 = ReconstructionSingle(1, "x", "b", single_xb1, it[7], wires[7], times[7], used[7], x1_single, sigm_x1_single);

    nDC1_segments = BuildUVSegments(1, pair_u1, pair_v1, pair_x1, pair_y1, single_ua1, single_ub1, single_va1, single_vb1,
            x1_ab, y1_ab, u1_ab, v1_ab, sigm_x1_ab, sigm_y1_ab, sigm_u1_ab, sigm_v1_ab, rh_segDC1, rh_sigm_segDC1, u1_single, v1_single, sigm_u1_single, sigm_v1_single);

    nDC1_segments = BuildXYSegments(1, pair_u1, pair_v1, pair_x1, pair_y1, single_xa1, single_xb1, single_ya1, single_yb1,
            x1_ab, y1_ab, u1_ab, v1_ab, sigm_x1_ab, sigm_y1_ab, sigm_u1_ab, sigm_v1_ab, rh_segDC1, rh_sigm_segDC1, x1_single, y1_single, sigm_x1_single, sigm_y1_single);

    pair_v2 = Reconstruction(2, "v", pair_v2, it[8], it[9], wires[8], wires[9], times[8], times[9], used[8], used[9], v2_ab, sigm_v2_ab);
    pair_u2 = Reconstruction(2, "u", pair_u2, it[10], it[11], wires[10], wires[11], times[10], times[11], used[10], used[11], u2_ab, sigm_u2_ab);
    pair_y2 = Reconstruction(2, "y", pair_y2, it[12], it[13], wires[12], wires[13], times[12], times[13], used[12], used[13], y2_ab, sigm_y2_ab);
    pair_x2 = Reconstruction(2, "x", pair_x2, it[14], it[15], wires[14], wires[15], times[14], times[15], used[14], used[15], x2_ab, sigm_x2_ab);

    single_va2 = ReconstructionSingle(2, "v", "a", single_va2, it[8], wires[8], times[8], used[8], v2_single, sigm_v2_single);
    single_vb2 = ReconstructionSingle(2, "v", "b", single_vb2, it[9], wires[9], times[9], used[9], v2_single, sigm_v2_single);
    single_ua2 = ReconstructionSingle(2, "u", "a", single_ua2, it[10], wires[10], times[10], used[10], u2_single, sigm_u2_single);
    single_ub2 = ReconstructionSingle(2, "u", "b", single_ub2, it[11], wires[11], times[11], used[11], u2_single, sigm_u2_single);
    single_ya2 = ReconstructionSingle(2, "y", "a", single_ya2, it[12], wires[12], times[12], used[12], y2_single, sigm_y2_single);
    single_yb2 = ReconstructionSingle(2, "y", "b", single_yb2, it[13], wires[13], times[13], used[13], y2_single, sigm_y2_single);
    single_xa2 = ReconstructionSingle(2, "x", "a", single_xa2, it[14], wires[14], times[14], used[14], x2_single, sigm_x2_single);
    single_xb2 = ReconstructionSingle(2, "x", "b", single_xb2, it[15], wires[15], times[15], used[15], x2_single, sigm_x2_single);

    // Build segments
    nDC2_segments = BuildUVSegments(2, pair_u2, pair_v2, pair_x2, pair_y2, single_ua2, single_ub2, single_va2, single_vb2,
            x2_ab, y2_ab, u2_ab, v2_ab, sigm_x2_ab, sigm_y2_ab, sigm_u2_ab, sigm_v2_ab, rh_segDC2, rh_sigm_segDC2, u2_single, v2_single, sigm_u2_single, sigm_v2_single);

    nDC2_segments = BuildXYSegments(2, pair_u2, pair_v2, pair_x2, pair_y2, single_xa2, single_xb2, single_ya2, single_yb2,
            x2_ab, y2_ab, u2_ab, v2_ab, sigm_x2_ab, sigm_y2_ab, sigm_u2_ab, sigm_v2_ab, rh_segDC2, rh_sigm_segDC2, x2_single, y2_single, sigm_x2_single, sigm_y2_single);

    // Common procedures over dch1 and dch2 
    // Fit found segments
    FitDchSegments(1, size_segDC1, rh_segDC1, rh_sigm_segDC1, par_ab1, chi2_DC1, xDC1_glob, yDC1_glob); // Dch1
    FitDchSegments(2, size_segDC2, rh_segDC2, rh_sigm_segDC2, par_ab2, chi2_DC2, xDC2_glob, yDC2_glob); // Dch2

    // Leave only longest and best chi2 segments   
    SelectLongestAndBestSegments(1, size_segDC1, rh_segDC1, chi2_DC1); // Dch1
    SelectLongestAndBestSegments(2, size_segDC2, rh_segDC2, chi2_DC2); // Dch2

    // Fill local segments z,x,y global coords; x-slope; y-slope; Chi2; to be continued...
    CreateDchTrack(1, chi2_DC1, par_ab1, size_segDC1); // Dch1
    CreateDchTrack(2, chi2_DC2, par_ab2, size_segDC2); // Dch2

    // Try to match the reconstructed segments from the two chambers
    // Not used in this version
    if (!fSegmentMatching) {
        cout << "\n======================== DCH track finder exec finished ========================" << endl;
        return;
    }


    cout << "\n======================== DCH track finder exec finished ========================" << endl;
}

Int_t BmnDchTrackFinder::BuildXYSegments(Int_t dchID,
        Int_t pairU, Int_t pairV, Int_t pairX, Int_t pairY, Int_t single_xa, Int_t single_xb, Int_t single_ya, Int_t single_yb,
        Float_t** x_ab, Float_t** y_ab, Float_t** u_ab, Float_t** v_ab,
        Float_t** sigm_x_ab, Float_t** sigm_y_ab, Float_t** sigm_u_ab, Float_t** sigm_v_ab,
        Float_t** rh_seg, Float_t** rh_sigm_seg,
        Float_t** x_single, Float_t** y_single, Float_t** sigm_x_single, Float_t** sigm_y_single) {

    Double_t sqrt_2 = sqrt(2.);
    Double_t isqrt_2 = 1. / sqrt_2;

    Int_t nDC_segments = (dchID == 1) ? nDC1_segments : nDC2_segments;
    for (Int_t i = 0; i < pairU; i++) {
        if (nDC_segments > 48)
            break;
        Float_t u_coord = (u_ab[0][i] + u_ab[1][i]) / 2;

        for (Int_t j = 0; j < pairV; j++) {
            if (nDC_segments > 48)
                break;
            Float_t v_coord = (v_ab[0][j] + v_ab[1][j]) / 2;
            Bool_t foundX = kFALSE;
            Float_t x_est = isqrt_2 * (v_coord - u_coord);
            Float_t y_est = isqrt_2 * (u_coord + v_coord);
            if (pairX > 0) {
                Double_t dX_thresh = 1.5;
                for (Int_t k = 0; k < pairX; k++) {
                    Float_t x_coord = (x_ab[0][k] + x_ab[1][k]) / 2;
                    if (nDC_segments > 48)
                        break;
                    if (Abs(x_coord - x_est) > dX_thresh)
                        continue;
                    dX_thresh = Abs(x_coord - x_est);

                    rh_seg[0][nDC_segments] = x_ab[0][k];
                    rh_seg[1][nDC_segments] = x_ab[1][k];
                    rh_seg[4][nDC_segments] = u_ab[0][i];
                    rh_seg[5][nDC_segments] = u_ab[1][i];
                    rh_seg[6][nDC_segments] = v_ab[0][j];
                    rh_seg[7][nDC_segments] = v_ab[1][j];
                    rh_sigm_seg[0][nDC_segments] = sigm_x_ab[0][k];
                    rh_sigm_seg[1][nDC_segments] = sigm_x_ab[1][k];
                    rh_sigm_seg[4][nDC_segments] = sigm_u_ab[0][i];
                    rh_sigm_seg[5][nDC_segments] = sigm_u_ab[1][i];
                    rh_sigm_seg[6][nDC_segments] = sigm_v_ab[0][j];
                    rh_sigm_seg[7][nDC_segments] = sigm_v_ab[1][j];

                    foundX = kTRUE;
                    if (nDC_segments > 48)
                        break;
                }//k
            }//(pair_x2>0)

            Bool_t foundY = kFALSE;
            if (pairY > 0) {
                Double_t dY_thresh = 1.0;
                for (Int_t m = 0; m < pairY; m++) {
                    if (nDC_segments > 48)
                        break;
                    Float_t y_coord = (y_ab[0][m] + y_ab[1][m]) / 2;
                    if (Abs(y_coord - y_est) > dY_thresh)
                        continue;
                    dY_thresh = Abs(y_coord - y_est);
                    foundY = kTRUE;
                    rh_seg[2][nDC_segments] = y_ab[0][m];
                    rh_seg[3][nDC_segments] = y_ab[1][m];
                    rh_seg[4][nDC_segments] = u_ab[0][i];
                    rh_seg[5][nDC_segments] = u_ab[1][i];
                    rh_seg[6][nDC_segments] = v_ab[0][j];
                    rh_seg[7][nDC_segments] = v_ab[1][j];
                    rh_sigm_seg[2][nDC_segments] = sigm_y_ab[0][m];
                    rh_sigm_seg[3][nDC_segments] = sigm_y_ab[1][m];
                    rh_sigm_seg[4][nDC_segments] = sigm_u_ab[0][i];
                    rh_sigm_seg[5][nDC_segments] = sigm_u_ab[1][i];
                    rh_sigm_seg[6][nDC_segments] = sigm_v_ab[0][j];
                    rh_sigm_seg[7][nDC_segments] = sigm_v_ab[1][j];
                    if (!foundX) {
                        Float_t min_a = 999;
                        Float_t min_b = 999;
                        for (Int_t kk = 0; kk < single_xa; kk++) {
                            if (Abs(x_single[1][kk] - x_est) > 1.5)
                                continue; //????? 0.5 needs to be reviewed

                            if (Abs(x_single[0][kk] - x_est) < min_a) {
                                min_a = Abs(x_single[0][kk] - x_est);
                                rh_seg[0][nDC_segments] = x_single[0][kk];
                                rh_sigm_seg[0][nDC_segments] = sigm_x_single[0][kk];
                                foundX = kTRUE;
                            }
                        }//for kk                                                                                                                                                           
                        for (Int_t kk = 0; kk < single_xb; kk++) {

                            if (Abs(x_single[1][kk] - x_est) > 1.5)
                                continue; //????? 0.5 needs to be reviewed                                                                                     
                            if (Abs(x_single[1][kk] - x_est) < min_b) {
                                min_b = Abs(x_single[1][kk] - x_est);
                                rh_seg[1][nDC_segments] = x_single[1][kk];
                                rh_sigm_seg[1][nDC_segments] = sigm_x_single[1][kk];
                                foundX = kTRUE;
                            }
                        }//for kk 
                        if (nDC_segments > 48)
                            break;
                    }//!foundX
                }//m         
                if (foundX && !foundY) {
                    Float_t min_a = 999;
                    Float_t min_b = 999;
                    for (Int_t kk = 0; kk < single_ya; kk++) {
                        if (Abs(y_single[0][kk] - y_est) > 1.5)
                            continue; //????? 0.5 needs to be reviewed                                                                                     
                        if (Abs(y_single[0][kk] - y_est) < min_a) {
                            min_a = Abs(y_single[0][kk] - y_est);
                            rh_seg[2][nDC_segments] = y_single[0][kk];
                            rh_sigm_seg[2][nDC_segments] = sigm_y_single[0][kk];
                            foundY = kTRUE;
                        }
                    }//for kk                                                                                                                                                          
                    for (Int_t kk = 0; kk < single_yb; kk++) {
                        if (Abs(y_single[1][kk] - y_est) > 1.5)
                            continue; //????? 0.5 needs to be reviewed                                                                                     
                        if (Abs(y_single[1][kk] - y_est) < min_b) {
                            min_b = Abs(y_single[1][kk] - y_est);
                            rh_seg[3][nDC_segments] = y_single[1][kk];
                            rh_sigm_seg[3][nDC_segments] = sigm_y_single[1][kk];
                            foundY = kTRUE;
                        }
                    }//for kk 
                }
            }//(pair_y2>0)
            if (foundX || foundY) nDC_segments++;
        }
    }
    return nDC_segments;
}

Int_t BmnDchTrackFinder::BuildUVSegments(Int_t dchID, Int_t pairU, Int_t pairV, Int_t pairX, Int_t pairY, Int_t single_ua, Int_t single_ub, Int_t single_va, Int_t single_vb,
        Float_t** x_ab, Float_t** y_ab, Float_t** u_ab, Float_t** v_ab,
        Float_t** sigm_x_ab, Float_t** sigm_y_ab, Float_t** sigm_u_ab, Float_t** sigm_v_ab,
        Float_t** rh_seg, Float_t** rh_sigm_seg,
        Float_t** u_single, Float_t** v_single, Float_t** sigm_u_single, Float_t** sigm_v_single) {

    Double_t sqrt_2 = sqrt(2.);
    Double_t isqrt_2 = 1. / sqrt_2;

    Int_t nDC_segments = (dchID == 1) ? nDC1_segments : nDC2_segments;

    for (Int_t i = 0; i < pairX; i++) {
        if (nDC_segments > 48)
            break;
        Float_t x_coord = (x_ab[0][i] + x_ab[1][i]) / 2;
        Float_t XU = x_coord;
        Float_t XV = x_coord;

        for (Int_t j = 0; j < pairY; j++) {
            Float_t y_coord = (y_ab[0][j] + y_ab[1][j]) / 2;
            Float_t YU = y_coord;
            Float_t YV = y_coord;
            Bool_t foundU = kFALSE;
            Float_t u_est = isqrt_2 * (YU - XU);
            Float_t v_est = isqrt_2 * (YV + XV);

            Double_t dU_thresh = 1.3;
            for (Int_t k = 0; k < pairU; k++) {
                Float_t u_coord = (u_ab[0][k] + u_ab[1][k]) / 2;

                if (Abs(u_coord - u_est) > dU_thresh)
                    continue;
                dU_thresh = Abs(u_coord - u_est);

                rh_seg[0][nDC_segments] = x_ab[0][i];
                rh_seg[1][nDC_segments] = x_ab[1][i];
                rh_seg[2][nDC_segments] = y_ab[0][j];
                rh_seg[3][nDC_segments] = y_ab[1][j];
                rh_seg[4][nDC_segments] = u_ab[0][k];
                rh_seg[5][nDC_segments] = u_ab[1][k];
                rh_sigm_seg[0][nDC_segments] = sigm_x_ab[0][i];
                rh_sigm_seg[1][nDC_segments] = sigm_x_ab[1][i];
                rh_sigm_seg[2][nDC_segments] = sigm_y_ab[0][j];
                rh_sigm_seg[3][nDC_segments] = sigm_y_ab[1][j];
                rh_sigm_seg[4][nDC_segments] = sigm_u_ab[0][k];
                rh_sigm_seg[5][nDC_segments] = sigm_u_ab[1][k];

                foundU = kTRUE;
                if (nDC_segments > 48)
                    break;
            }

            Bool_t foundV = kFALSE;

            Double_t dV_thresh = 1.3;
            for (Int_t m = 0; m < pairV; m++) {
                if (nDC_segments > 48)
                    break;
                Float_t v_coord = (v_ab[0][m] + v_ab[1][m]) / 2;

                if (Abs(v_coord - v_est) > dV_thresh)
                    continue;
                dV_thresh = Abs(v_coord - v_est);

                foundV = kTRUE;
                rh_seg[0][nDC_segments] = x_ab[0][i];
                rh_seg[1][nDC_segments] = x_ab[1][i];
                rh_seg[2][nDC_segments] = y_ab[0][j];
                rh_seg[3][nDC_segments] = y_ab[1][j];
                rh_seg[6][nDC_segments] = v_ab[0][m];
                rh_seg[7][nDC_segments] = v_ab[1][m];
                rh_sigm_seg[0][nDC_segments] = sigm_x_ab[0][i];
                rh_sigm_seg[1][nDC_segments] = sigm_x_ab[1][i];
                rh_sigm_seg[2][nDC_segments] = sigm_y_ab[0][j];
                rh_sigm_seg[3][nDC_segments] = sigm_y_ab[1][j];
                rh_sigm_seg[6][nDC_segments] = sigm_v_ab[0][m];
                rh_sigm_seg[7][nDC_segments] = sigm_v_ab[1][m];

                if (!foundU) {
                    Float_t min_a = 999;
                    Float_t min_b = 999;
                    for (Int_t kk = 0; kk < single_ua; kk++) {
                        if (Abs(u_single[0][kk] - u_est) > 1.5)
                            continue; //????? 0.5 needs to be reviewed
                        if (Abs(u_single[0][kk] - u_est) < min_a) {
                            min_a = Abs(u_single[0][kk] - u_est);
                            rh_seg[4][nDC_segments] = u_single[0][kk];
                            rh_sigm_seg[4][nDC_segments] = sigm_u_single[0][kk];
                            foundU = kTRUE;
                        }
                    }//for kk
                    for (Int_t kk = 0; kk < single_ub; kk++) {
                        if (Abs(u_single[1][kk] - u_est) > 1.5)
                            continue; //????? 0.5 needs to be reviewed
                        if (Abs(u_single[1][kk] - u_est) < min_b) {
                            min_b = Abs(u_single[1][kk] - u_est);
                            rh_seg[5][nDC_segments] = u_single[1][kk];
                            rh_sigm_seg[5][nDC_segments] = sigm_u_single[1][kk];
                            foundU = kTRUE;
                        }
                    }//for kk
                    if (nDC_segments > 48)
                        break;
                }//!foundU

                if (nDC_segments > 48)
                    break;
            }//m
            //            }//(pair_v2>0)
            if (!foundV && foundU) {
                Float_t min_a = 999;
                Float_t min_b = 999;
                for (Int_t kk = 0; kk < single_va; kk++) {
                    if (Abs(v_single[0][kk] - v_est) > 1.5)
                        continue; //????? 0.5 needs to be reviewed                                                                                     
                    if (Abs(v_single[0][kk] - v_est) < min_a) {
                        min_a = Abs(v_single[0][kk] - v_est);
                        rh_seg[6][nDC_segments] = v_single[0][kk];
                        rh_sigm_seg[6][nDC_segments] = sigm_v_single[0][kk];
                        foundV = kTRUE;
                    }
                }//for kk                                                                                                                                                           
                for (Int_t kk = 0; kk < single_vb; kk++) {
                    if (Abs(v_single[1][kk] - v_est) > 1.5)
                        continue; //????? 0.5 needs to be reviewed                                                                                     
                    if (Abs(v_single[1][kk] - v_est) < min_b) {
                        min_b = Abs(v_single[1][kk] - v_est);
                        rh_seg[7][nDC_segments] = v_single[1][kk];
                        rh_sigm_seg[7][nDC_segments] = sigm_v_single[1][kk];
                        foundV = kTRUE;
                    }
                }
            }
            if (foundV || foundU) nDC_segments++;
        }
    }
    return nDC_segments;
}

void BmnDchTrackFinder::FitDchSegments(Int_t dchID, Int_t* size_seg, Float_t** rh_seg, Float_t** rh_sigm_seg, Float_t** par_ab, Float_t* chi2, Float_t* x_glob, Float_t* y_glob) {
    Int_t nDC_segments = (dchID == 1) ? nDC1_segments : nDC2_segments;
    for (Int_t j = 0; j < nDC_segments; j++) {
        Int_t worst_hit = -1;
        Double_t max_resid = 0;

        Float_t _rh_seg[8];
        Float_t _rh_sigm_seg[8];
        Float_t _par_ab[4];

        for (Int_t i = 0; i < 8; i++)
            if (Abs(rh_seg[i][j] + 999.) > FLT_EPSILON)
                size_seg[j]++;

        for (Int_t rej = 0; rej < 2; rej++) {//allow 2 passes max 8->7 & 7->6
            for (Int_t i = 0; i < 8; i++) {
                _rh_seg[i] = rh_seg[i][j];
                _rh_sigm_seg[i] = rh_sigm_seg[i][j];
            }

            fit_seg(z_loc, _rh_seg, _rh_sigm_seg, _par_ab, -1, -1); //usual fit without skipping any plane
            for (Int_t i = 0; i < 4; i++)
                par_ab[i][j] = _par_ab[i];

            chi2[j] = 0;

            Float_t resid(LDBL_MAX);
            for (Int_t i = 0; i < 8; i++) {
                if (Abs(rh_seg[i][j] + 999.) < FLT_EPSILON)
                    continue;

                resid = CalculateResidual(i, j, rh_seg, par_ab);
                chi2[j] += (resid * resid) / rh_sigm_seg[i][j];
                if (Abs(resid) > max_resid) {
                    worst_hit = i;
                    max_resid = Abs(resid);
                }
            }

            chi2[j] /= (size_seg[j] - 4);

            //if chi2 is big and seg_size = min erase this seg
            if (chi2[j] > 30.)
                if (size_seg[j] == 6) {
                    chi2[j] = 999.;
                    break;
                } else {
                    rh_seg[worst_hit][j] = -999.; //erase worst hit and refit
                    size_seg[j]--;
                    max_resid = 0;
                    continue;
                }
        }

        // Add shifts to slopes and coords
        Float_t x_slope_sh = (dchID == 1) ? x1_slope_sh : x2_slope_sh;
        Float_t y_slope_sh = (dchID == 1) ? y1_slope_sh : y2_slope_sh;
        Float_t x_sh = (dchID == 1) ? x1_sh : x2_sh;
        Float_t y_sh = (dchID == 1) ? y1_sh : y2_sh;

        par_ab[0][j] += x_slope_sh + x_slope_sh * par_ab[0][j] * par_ab[0][j];
        par_ab[2][j] += y_slope_sh + y_slope_sh * par_ab[2][j] * par_ab[2][j];
        par_ab[1][j] += x_sh;
        par_ab[3][j] += y_sh;

        x_glob[j] = par_ab[0][j]*(-99.5) + par_ab[1][j];
        y_glob[j] = par_ab[2][j]*(-99.5) + par_ab[3][j];

        Bool_t has7DC = (dchID == 1) ? has7DC1 : has7DC2;
        if (size_seg[j] > 6)
            has7DC = kTRUE;
    }
}

void BmnDchTrackFinder::CompareDaDb(Float_t d, Float_t& ele) {
    ele = (d < 0.02) ? (0.08 * 0.08) :
            (d >= 0.02 && d < 0.1) ? (0.06 * 0.06) :
            (d >= 0.1 && d < 0.4) ? (0.025 * 0.025) :
            (d >= 0.4 && d < 0.41) ? (0.08 * 0.08) :
            (0.10 * 0.10);
}

void BmnDchTrackFinder::CompareDaDb(Float_t d, Float_t& ele1, Float_t& ele2) {
    ele1 = (d < 0.02) ? (0.08 * 0.08) :
            (d >= 0.02 && d < 0.1) ? (0.06 * 0.06) :
            (d >= 0.1 && d < 0.4) ? (0.025 * 0.025) :
            (d >= 0.4 && d < 0.41) ? (0.08 * 0.08) :
            (0.10 * 0.10);
    ele2 = ele1;
}

void BmnDchTrackFinder::SelectLongestAndBestSegments(Int_t dchID, Int_t* size_seg, Float_t** rh_seg, Float_t* chi2) {
    Int_t nDC_segments = (dchID == 1) ? nDC1_segments : nDC2_segments;
    for (Int_t max_size = 8; max_size > 5; max_size--) 
        for (Int_t it1 = 0; it1 < nDC_segments; it1++) {
            if (size_seg[it1] != max_size)
                continue;
            for (Int_t it2 = 0; it2 < nDC_segments; it2++) {
                if (it2 == it1)
                    continue;
                for (Int_t hit = 0; hit < 4; hit++) 
                    if (rh_seg[2 * hit][it1] == rh_seg[2 * hit][it2] &&
                            rh_seg[2 * hit + 1][it1] == rh_seg[2 * hit + 1][it2] &&
                            (chi2[it1] <= chi2[it2] || size_seg[it1] > size_seg[it2])) {
                        chi2[it2] = 999.;                                                                                                   
                        break;
                    }
            }
        }
}

void BmnDchTrackFinder::CreateDchTrack(Int_t dchID, Float_t* chi2Arr, Float_t** parArr, Int_t* sizeArr) {
    Int_t nDC_segments = (dchID == 1) ? nDC1_segments : nDC2_segments;
    for (Int_t iSegment = 0; iSegment < nDC_segments; iSegment++) {
        if (chi2Arr[iSegment] > 50)
            continue;
        FairTrackParam trackParam;
        Float_t z0 = (dchID == 1) ? Z_dch1 : Z_dch2;
        Float_t x0 = parArr[1][iSegment];
        Float_t y0 = parArr[3][iSegment];
        trackParam.SetPosition(TVector3(-x0, y0, z0));
        trackParam.SetTx(-parArr[0][iSegment]);
        trackParam.SetTy(parArr[2][iSegment]);

        BmnDchTrack* track = new((*fDchTracks)[fDchTracks->GetEntriesFast()]) BmnDchTrack();
        track->SetChi2(chi2Arr[iSegment]);
        track->SetNHits(sizeArr[iSegment]);
        track->SetParamFirst(trackParam);
    }
}

Float_t BmnDchTrackFinder::CalculateResidual(Int_t i, Int_t j, Float_t** rh_seg, Float_t** par_ab) {
    Double_t sqrt_2 = sqrt(2.);
    Double_t isqrt_2 = 1 / sqrt_2;

    return (i < 2) ? rh_seg[i][j] - z_loc[i] * par_ab[0][j] - par_ab[1][j] :
            (i >= 2 && i < 4) ? rh_seg[i][j] - z_loc[i] * par_ab[2][j] - par_ab[3][j] :
            (i >= 4 && i < 6) ? rh_seg[i][j] - isqrt_2 * z_loc[i] * (par_ab[2][j] - par_ab[0][j]) - isqrt_2 * (par_ab[3][j] - par_ab[1][j]) :
            rh_seg[i][j] - isqrt_2 * z_loc[i] * (par_ab[2][j] + par_ab[0][j]) - isqrt_2 * (par_ab[3][j] + par_ab[1][j]);
}

InitStatus BmnDchTrackFinder::Init() {
    cout << endl << "BmnDchTrackFinder::Init()" << endl;
    FairRootManager* ioman = FairRootManager::Instance();

    fBmnDchDigitsArray = (TClonesArray*) ioman->GetObject(InputDigitsBranchName);

    fDchTracks = new TClonesArray(tracksDch.Data());
    ioman->Register(tracksDch.Data(), "DCH", fDchTracks, kTRUE);

    ifstream fin;
    TString dir = getenv("VMCWORKDIR");
    dir += "/input/";
    fin.open((TString(dir + "transfer_func.txt")).Data(), ios::in);
    for (Int_t fi = 0; fi < 16; fi++) {
        fin >> t_dc[0][fi] >> t_dc[1][fi] >> t_dc[2][fi] >> t_dc[3][fi] >> t_dc[4][fi] >>
                pol_par_dc[0][0][fi] >> pol_par_dc[0][1][fi] >> pol_par_dc[0][2][fi] >> pol_par_dc[0][3][fi] >> pol_par_dc[0][4][fi] >>
                pol_par_dc[1][0][fi] >> pol_par_dc[1][1][fi] >> pol_par_dc[1][2][fi] >> pol_par_dc[1][3][fi] >> pol_par_dc[1][4][fi] >>
                pol_par_dc[2][0][fi] >> pol_par_dc[2][1][fi] >> pol_par_dc[2][2][fi] >> pol_par_dc[2][3][fi] >> pol_par_dc[2][4][fi] >>
                scaling[fi];
    }

    fin.close();

    const Int_t N = 2;

    // z local xa->vb (cm) 
    Double_t arr1[4 * N] = {9.3, 8.1, 3.5, 2.3, -2.3, -3.5, -8.1, -9.3};
    for (Int_t iSize = 0; iSize < 4 * N; iSize++)
        z_loc[iSize] = arr1[iSize];

    // z global dc 1 & dc 2 (cm)
    Double_t arr2[8 * N] = {-45.7, -46.9, -51.5, -52.7, -57.3, -58.5, -63.1, -64.3, 64.3, 63.1, 58.5, 57.3, 52.7, 51.5, 46.9, 45.7};
    for (Int_t iSize = 0; iSize < 8 * N; iSize++)
        z_glob[iSize] = arr2[iSize];

    x1_ab = new Float_t*[N];
    y1_ab = new Float_t*[N];
    u1_ab = new Float_t*[N];
    v1_ab = new Float_t*[N];
    sigm_x1_ab = new Float_t*[N];
    sigm_y1_ab = new Float_t*[N];
    sigm_u1_ab = new Float_t*[N];
    sigm_v1_ab = new Float_t*[N];
    x2_ab = new Float_t*[N];
    y2_ab = new Float_t*[N];
    u2_ab = new Float_t*[N];
    v2_ab = new Float_t*[N];
    sigm_x2_ab = new Float_t*[N];
    sigm_y2_ab = new Float_t*[N];
    sigm_u2_ab = new Float_t*[N];
    sigm_v2_ab = new Float_t*[N];
    for (Int_t iDim = 0; iDim < N; iDim++) {
        x1_ab[iDim] = new Float_t[75 * N];
        y1_ab[iDim] = new Float_t[75 * N];
        u1_ab[iDim] = new Float_t[75 * N];
        v1_ab[iDim] = new Float_t[75 * N];
        sigm_x1_ab[iDim] = new Float_t[75 * N];
        sigm_y1_ab[iDim] = new Float_t[75 * N];
        sigm_u1_ab[iDim] = new Float_t[75 * N];
        sigm_v1_ab[iDim] = new Float_t[75 * N];
        x2_ab[iDim] = new Float_t[75 * N];
        y2_ab[iDim] = new Float_t[75 * N];
        u2_ab[iDim] = new Float_t[75 * N];
        v2_ab[iDim] = new Float_t[75 * N];
        sigm_x2_ab[iDim] = new Float_t[75 * N];
        sigm_y2_ab[iDim] = new Float_t[75 * N];
        sigm_u2_ab[iDim] = new Float_t[75 * N];
        sigm_v2_ab[iDim] = new Float_t[75 * N];
    }

    par_ab1 = new Float_t*[2 * N];
    par_ab2 = new Float_t*[2 * N];
    for (Int_t iDim = 0; iDim < 2 * N; iDim++) {
        par_ab1[iDim] = new Float_t[75 * N];
        par_ab2[iDim] = new Float_t[75 * N];
    }

    chi2_DC1 = new Float_t[75 * N];
    chi2_DC2 = new Float_t[75 * N];
    xDC1_glob = new Float_t[75 * N];
    yDC1_glob = new Float_t[75 * N];
    xDC2_glob = new Float_t[75 * N];
    yDC2_glob = new Float_t[75 * N];

    size_segDC1 = new Int_t[75 * N];
    size_segDC2 = new Int_t[75 * N];

    rh_segDC1 = new Float_t*[4 * N];
    rh_segDC2 = new Float_t*[4 * N];
    rh_sigm_segDC1 = new Float_t*[4 * N];
    rh_sigm_segDC2 = new Float_t*[4 * N];
    for (Int_t iDim = 0; iDim < 4 * N; iDim++) {
        rh_segDC1[iDim] = new Float_t[75 * N];
        rh_segDC2[iDim] = new Float_t[75 * N];
        rh_sigm_segDC1[iDim] = new Float_t[75 * N];
        rh_sigm_segDC2[iDim] = new Float_t[75 * N];
    }

    //single hits on ab-plane
    x1_single = new Float_t*[N];
    y1_single = new Float_t*[N];
    u1_single = new Float_t*[N];
    v1_single = new Float_t*[N];
    sigm_x1_single = new Float_t*[N];
    sigm_y1_single = new Float_t*[N];
    sigm_u1_single = new Float_t*[N];
    sigm_v1_single = new Float_t*[N];
    x2_single = new Float_t*[N];
    y2_single = new Float_t*[N];
    u2_single = new Float_t*[N];
    v2_single = new Float_t*[N];
    sigm_x2_single = new Float_t*[N];
    sigm_y2_single = new Float_t*[N];
    sigm_u2_single = new Float_t*[N];
    sigm_v2_single = new Float_t*[N];
    for (Int_t iDim = 0; iDim < N; iDim++) {
        x1_single[iDim] = new Float_t[20 * N];
        y1_single[iDim] = new Float_t[20 * N];
        u1_single[iDim] = new Float_t[20 * N];
        v1_single[iDim] = new Float_t[20 * N];
        sigm_x1_single[iDim] = new Float_t[20 * N];
        sigm_y1_single[iDim] = new Float_t[20 * N];
        sigm_u1_single[iDim] = new Float_t[20 * N];
        sigm_v1_single[iDim] = new Float_t[20 * N];
        x2_single[iDim] = new Float_t[20 * N];
        y2_single[iDim] = new Float_t[20 * N];
        u2_single[iDim] = new Float_t[20 * N];
        v2_single[iDim] = new Float_t[20 * N];
        sigm_x2_single[iDim] = new Float_t[20 * N];
        sigm_y2_single[iDim] = new Float_t[20 * N];
        sigm_u2_single[iDim] = new Float_t[20 * N];
        sigm_v2_single[iDim] = new Float_t[20 * N];
    }
}

void BmnDchTrackFinder::PrepareArraysToProcessEvent() {
    fDchTracks->Clear();
    has7DC1 = kFALSE;
    has7DC2 = kFALSE;
    nDC1_segments = 0;
    nDC2_segments = 0;
    // Array cleaning and initializing 
    for (Int_t iDim1 = 0; iDim1 < 4; iDim1++)
        for (Int_t iDim2 = 0; iDim2 < 150; iDim2++) {
            par_ab1[iDim1][iDim2] = -999.;
            par_ab2[iDim1][iDim2] = -999.;
        }

    for (Int_t iDim = 0; iDim < 150; iDim++) {
        chi2_DC1[iDim] = 50.;
        chi2_DC2[iDim] = 50.;
        size_segDC1[iDim] = 0;
        size_segDC2[iDim] = 0;
        xDC1_glob[iDim] = -999.;
        yDC1_glob[iDim] = -999.;
        xDC2_glob[iDim] = -999.;
        yDC2_glob[iDim] = -999.;
    }

    for (Int_t iDim1 = 0; iDim1 < 8; iDim1++)
        for (Int_t iDim2 = 0; iDim2 < 150; iDim2++) {
            rh_segDC1[iDim1][iDim2] = -999.;
            rh_segDC2[iDim1][iDim2] = -999.;
            rh_sigm_segDC1[iDim1][iDim2] = 1.;
            rh_sigm_segDC2[iDim1][iDim2] = 1.;
        }

    for (Int_t iDim1 = 0; iDim1 < 2; iDim1++) {
        for (Int_t iDim2 = 0; iDim2 < 150; iDim2++) {
            x1_ab[iDim1][iDim2] = -999.;
            y1_ab[iDim1][iDim2] = -999.;
            u1_ab[iDim1][iDim2] = -999.;
            v1_ab[iDim1][iDim2] = -999.;
            sigm_x1_ab[iDim1][iDim2] = 1.;
            sigm_y1_ab[iDim1][iDim2] = 1.;
            sigm_u1_ab[iDim1][iDim2] = 1.;
            sigm_v1_ab[iDim1][iDim2] = 1.;
            x2_ab[iDim1][iDim2] = -999.;
            y2_ab[iDim1][iDim2] = -999.;
            u2_ab[iDim1][iDim2] = -999.;
            v2_ab[iDim1][iDim2] = -999.;
            sigm_x2_ab[iDim1][iDim2] = 1.;
            sigm_y2_ab[iDim1][iDim2] = 1.;
            sigm_u2_ab[iDim1][iDim2] = 1.;
            sigm_v2_ab[iDim1][iDim2] = 1.;
        }
        for (Int_t iDim2 = 0; iDim2 < 40; iDim2++) {
            x1_single[iDim1][iDim2] = -999.;
            y1_single[iDim1][iDim2] = -999.;
            u1_single[iDim1][iDim2] = -999.;
            v1_single[iDim1][iDim2] = -999.;
            sigm_x1_single[iDim1][iDim2] = 1.;
            sigm_y1_single[iDim1][iDim2] = 1.;
            sigm_u1_single[iDim1][iDim2] = 1.;
            sigm_v1_single[iDim1][iDim2] = 1.;
            x2_single[iDim1][iDim2] = -999.;
            y2_single[iDim1][iDim2] = -999.;
            u2_single[iDim1][iDim2] = -999.;
            v2_single[iDim1][iDim2] = -999.;
            sigm_x2_single[iDim1][iDim2] = 1.;
            sigm_y2_single[iDim1][iDim2] = 1.;
            sigm_u2_single[iDim1][iDim2] = 1.;
            sigm_v2_single[iDim1][iDim2] = 1.;
        }
    }
    pair_x2 = 0;
    pair_y2 = 0;
    pair_u2 = 0;
    pair_v2 = 0;
    single_xa2 = 0;
    single_ya2 = 0;
    single_ua2 = 0;
    single_va2 = 0;
    single_xb2 = 0;
    single_yb2 = 0;
    single_ub2 = 0;
    single_vb2 = 0;
    pair_x1 = 0;
    pair_y1 = 0;
    pair_u1 = 0;
    pair_v1 = 0;
    single_xa1 = 0;
    single_ya1 = 0;
    single_ua1 = 0;
    single_va1 = 0;
    single_xb1 = 0;
    single_yb1 = 0;
    single_ub1 = 0;
    single_vb1 = 0;
}

void BmnDchTrackFinder::Finish() {
    for (Int_t i = 0; i < 4; i++) {
        delete [] par_ab1[i];
        delete [] par_ab2[i];
    }

    delete[] par_ab1;
    delete[] par_ab2;

    delete[] chi2_DC1;
    delete[] chi2_DC2;

    delete[] size_segDC1;
    delete[] size_segDC2;

    for (Int_t i = 0; i < 8; i++) {
        delete [] rh_segDC1[i];
        delete [] rh_segDC2[i];
    }

    delete[] rh_segDC1;
    delete[] rh_segDC2;

    delete [] xDC1_glob;
    delete [] xDC2_glob;
    delete [] yDC1_glob;
    delete [] yDC2_glob;

    for (Int_t i = 0; i < 2; i++) {
        delete [] x1_ab[i];
        delete [] y1_ab[i];
        delete [] u1_ab[i];
        delete [] v1_ab[i];
        delete [] sigm_x1_ab[i];
        delete [] sigm_y1_ab[i];
        delete [] sigm_u1_ab[i];
        delete [] sigm_v1_ab[i];
        delete [] x2_ab[i];
        delete [] y2_ab[i];
        delete [] u2_ab[i];
        delete [] v2_ab[i];
        delete [] sigm_x2_ab[i];
        delete [] sigm_y2_ab[i];
        delete [] sigm_u2_ab[i];
        delete [] sigm_v2_ab[i];

        delete [] x1_single[i];
        delete [] y1_single[i];
        delete [] u1_single[i];
        delete [] v1_single[i];
        delete [] sigm_x1_single[i];
        delete [] sigm_y1_single[i];
        delete [] sigm_u1_single[i];
        delete [] sigm_v1_single[i];
        delete [] x2_single[i];
        delete [] y2_single[i];
        delete [] u2_single[i];
        delete [] v2_single[i];
        delete [] sigm_x2_single[i];
        delete [] sigm_y2_single[i];
        delete [] sigm_u2_single[i];
        delete [] sigm_v2_single[i];
    }

    delete [] x1_ab;
    delete [] y1_ab;
    delete [] u1_ab;
    delete [] v1_ab;
    delete [] sigm_x1_ab;
    delete [] sigm_y1_ab;
    delete [] sigm_u1_ab;
    delete [] sigm_v1_ab;
    delete [] x2_ab;
    delete [] y2_ab;
    delete [] u2_ab;
    delete [] v2_ab;
    delete [] sigm_x2_ab;
    delete [] sigm_y2_ab;
    delete [] sigm_u2_ab;
    delete [] sigm_v2_ab;

    delete [] x1_single;
    delete [] y1_single;
    delete [] u1_single;
    delete [] v1_single;
    delete [] sigm_x1_single;
    delete [] sigm_y1_single;
    delete [] sigm_u1_single;
    delete [] sigm_v1_single;
    delete [] x2_single;
    delete [] y2_single;
    delete [] u2_single;
    delete [] v2_single;
    delete [] sigm_x2_single;
    delete [] sigm_y2_single;
    delete [] sigm_u2_single;
    delete [] sigm_v2_single;
}

Int_t BmnDchTrackFinder::Reconstruction(Int_t dchID, TString wire, Int_t pair, Int_t it_a, Int_t it_b,
        Double_t* wirenr_a, Double_t* wirenr_b, Double_t* time_a, Double_t* time_b,
        Bool_t* used_a, Bool_t* used_b,
        Float_t** _ab, Float_t** sigm_ab) {

    const Int_t arrIdxShift = (dchID == 2) ? 8 : 0;
    const Int_t arrIdxStart = (wire == "x") ? 0 : (wire == "y") ? 2 : (wire == "u") ? 4 : 6;

    Float_t a_pm[2], b_pm[2];

    for (Int_t i = 0; i < it_a; ++i)
        for (Int_t j = 0; j < it_b; ++j) {
            if (pair > 48)
                break;
            if ((wirenr_a[i] != wirenr_b[j] && wirenr_a[i] != wirenr_b[j] + 1))
                continue;
            Int_t func_nr_a = -1;
            Int_t func_nr_b = -1;
            for (Int_t t_it = 0; t_it < 4; t_it++)
                if (time_a[i] >= t_dc[t_it][0 + arrIdxStart + arrIdxShift] && time_a[i] < t_dc[t_it + 1][0 + arrIdxStart + arrIdxShift]) {
                    func_nr_a = t_it;
                    break;
                }

            Double_t time = time_a[i];
            Double_t d_a = 0;
            Double_t d_b = 0;

            if (func_nr_a == 1 || func_nr_a == 2) d_a = scale * (pol_par_dc[1][0][0 + arrIdxStart + arrIdxShift] + pol_par_dc[1][1][0 + arrIdxStart + arrIdxShift] * time +
                    pol_par_dc[1][2][0 + arrIdxStart + arrIdxShift] * Power(time, 2) +
                    pol_par_dc[1][3][0 + arrIdxStart + arrIdxShift] * Power(time, 3) +
                    pol_par_dc[1][4][0 + arrIdxStart + arrIdxShift] * Power(time, 4));

            else if (func_nr_a == 0) d_a = 0;

            else if (func_nr_a == 3) d_a = scale * (pol_par_dc[2][0][0 + arrIdxStart + arrIdxShift] + pol_par_dc[2][1][0 + arrIdxStart + arrIdxShift] * time +
                    pol_par_dc[2][2][0 + arrIdxStart + arrIdxShift] * Power(time, 2) +
                    pol_par_dc[2][3][0 + arrIdxStart + arrIdxShift] * Power(time, 3) +
                    pol_par_dc[2][4][0 + arrIdxStart + arrIdxShift] * Power(time, 4));

            for (Int_t t_it = 0; t_it < 4; t_it++) 
                if (time_b[j] >= t_dc[t_it][1 + arrIdxStart + arrIdxShift] && time_b[j] < t_dc[t_it + 1][1 + arrIdxStart + arrIdxShift]) {
                    func_nr_b = t_it;
                    break;
                }

            time = time_b[j];

            if (func_nr_b == 1 || func_nr_b == 2) d_b = scale * (pol_par_dc[1][0][1 + arrIdxStart + arrIdxShift] + pol_par_dc[1][1][1 + arrIdxStart + arrIdxShift] * time +
                    pol_par_dc[1][2][1 + arrIdxStart + arrIdxShift] * Power(time, 2) +
                    pol_par_dc[1][3][1 + arrIdxStart + arrIdxShift] * Power(time, 3) +
                    pol_par_dc[1][4][1 + arrIdxStart + arrIdxShift] * Power(time, 4));

            else if (func_nr_b == 0) d_b = 0;

            else if (func_nr_b == 3) d_b = scale * (pol_par_dc[2][0][1 + arrIdxStart + arrIdxShift] + pol_par_dc[2][1][1 + arrIdxStart + arrIdxShift] * time +
                    pol_par_dc[2][2][1 + arrIdxStart + arrIdxShift] * Power(time, 2) +
                    pol_par_dc[2][3][1 + arrIdxStart + arrIdxShift] * Power(time, 3) +
                    pol_par_dc[2][4][1 + arrIdxStart + arrIdxShift] * Power(time, 4));

            a_pm[0] = wirenr_a[i] - 119 + d_a;
            a_pm[1] = wirenr_a[i] - 119 - d_a;
            b_pm[0] = wirenr_b[j] - 118.5 + d_b;
            b_pm[1] = wirenr_b[j] - 118.5 - d_b;

            Double_t dmin = LDBL_MAX;
            for (Int_t k = 0; k < 2; k++)
                for (Int_t m = 0; m < 2; m++)
                    if (Abs(a_pm[k] - b_pm[m]) < dmin) {
                        dmin = Abs(a_pm[k] - b_pm[m]);
                        _ab[0][pair] = a_pm[k];
                        _ab[1][pair] = b_pm[m];
                    }

            CompareDaDb(d_a, sigm_ab[0][pair]);
            CompareDaDb(d_b, sigm_ab[1][pair]);

            pair++;

            used_a[i] = kTRUE;
            used_b[j] = kTRUE;
        }
    return pair;
}

Int_t BmnDchTrackFinder::ReconstructionSingle(Int_t dchID, TString wire, TString lay, Int_t single, Int_t it,
        Double_t* wirenr, Double_t* time_, Bool_t* used,
        Float_t** _single, Float_t** sigm_single) {

    const Int_t arrIdxStart = (wire == "x") ? 0 : (wire == "y") ? 2 : (wire == "u") ? 4 : 6;

    const Int_t arrIdx1 = (lay == "a") ? 0 : 1;
    const Int_t arrIdx2 = (dchID == 2) ? 8 : 0;
    const Double_t coeff = (lay == "a") ? 119 : 118.5;

    for (Int_t i = 0; i < it; ++i) {
        if (used[i])
            continue;

        Int_t func_nr = -1;
        for (Int_t t_it = 0; t_it < 4; t_it++) {
            if (time_[i] >= t_dc[t_it][0 + arrIdxStart + arrIdx1 + arrIdx2] && time_[i] < t_dc[t_it + 1][0 + arrIdxStart + arrIdx1 + arrIdx2]) {
                func_nr = t_it;
                break;
            }
        }
        Double_t time = time_[i];
        Double_t d = 0;

        if (func_nr == 1 || func_nr == 2) d = scale * (pol_par_dc[1][0][0 + arrIdx1 + arrIdx2] + pol_par_dc[1][1][0 + arrIdx1 + arrIdx2] * time +
                pol_par_dc[1][2][0 + arrIdxStart + arrIdx1 + arrIdx2] * Power(time, 2) +
                pol_par_dc[1][3][0 + arrIdxStart + arrIdx1 + arrIdx2] * Power(time, 3) +
                pol_par_dc[1][4][0 + arrIdxStart + arrIdx1 + arrIdx2] * Power(time, 4));

        else if (func_nr == 0) d = 0;

        else if (func_nr == 3) d = scale * (pol_par_dc[2][0][0 + arrIdx1 + arrIdx2] + pol_par_dc[2][1][0 + arrIdx1 + arrIdx2] * time +
                pol_par_dc[2][2][0 + arrIdxStart + arrIdx1 + arrIdx2] * Power(time, 2) +
                pol_par_dc[2][3][0 + arrIdxStart + arrIdx1 + arrIdx2] * Power(time, 3) +
                pol_par_dc[2][4][0 + arrIdxStart + arrIdx1 + arrIdx2] * Power(time, 4));

        _single[0 + arrIdx1][single] = wirenr[i] - coeff + d;
        _single[0 + arrIdx1][single + 1] = wirenr[i] - coeff - d;

        CompareDaDb(d, sigm_single[0 + arrIdx1][single], sigm_single[0 + arrIdx1][single + 1]);

        single += 2;
    }
    return single;
}
ClassImp(BmnDchTrackFinder)

