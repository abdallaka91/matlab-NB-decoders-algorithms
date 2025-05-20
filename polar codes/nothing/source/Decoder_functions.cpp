
#include "Decoder_functions.h"
#include "struct.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <vector>
#include <iostream>
#include <queue>
#include <limits>
#include <cmath>

void PoAwN::decoding::Channel_LLR(const vector<vector<softdata_t>> &chan_observ,
                                  const vector<vector<uint16_t>> &bin_symb_seq,
                                  uint16_t q,
                                  softdata_t sigma,
                                  vector<vector<softdata_t>> &chan_LLR)
{
    const int N = chan_observ.size();
    const int q1 = chan_observ[0].size();

    constexpr const softdata_t two_pow16 = softdata_t(1 << 16);

    vector<softdata_t> hard_decison(q1, two_pow16 - 1);

    const softdata_t fct = softdata_t(2) / (sigma * sigma);
    for (int i = 0; i < N; i++)
    {
        softdata_t mn_llr = std::numeric_limits<softdata_t>::max() / 2;
        for (int j = 0; j < q; j++)
        {
            softdata_t temp = 0;

            for (int k = 0; k < q1; k++)
            {
                // true == 1 and false == 0
                hard_decison[k] = softdata_t(uint8_t(chan_observ[i][k] <= 0));
                temp += chan_observ[i][k] * ((softdata_t)bin_symb_seq[j][k] - hard_decison[k]);
            }
            temp *= fct;
            chan_LLR[i][j] = temp;
            if (temp < mn_llr)
                mn_llr = temp;
        }
        if (mn_llr > std::numeric_limits<softdata_t>::min())
            for (int j = 0; j < q; j++)
            {
                chan_LLR[i][j] -= mn_llr;
            }
    }
}

// void PoAwN::decoding::ECN_EMS_bubble(const decoder_t &theta_1,
//                                      const decoder_t &phi_1,
//                                      const decoder_parameters &dec_param,
//                                      uint16_t coef,
//                                      const vector<vector<uint16_t>> &ADDGF,
//                                      const vector<vector<uint16_t>> &DIVGF,
//                                      decoder_t &theta)
// {
//     theta.intrinsic_LLR.reserve(dec_param.nm);
//     theta.intrinsic_GF.reserve(dec_param.nm);
//     bool rel_theta = (theta_1.intrinsic_LLR[dec_param.Zc] > phi_1.intrinsic_LLR[dec_param.Zc]);
//     decoder_t phi_1_p = phi_1;
//     for (int i = 0; i < phi_1.intrinsic_GF.size(); i++)
//         phi_1_p.intrinsic_GF[i] = DIVGF[phi_1_p.intrinsic_GF[i]][coef];
//     vector<softdata_t> srtr(dec_param.nb);
//     vector<uint16_t> srtrg(dec_param.nb);
//     vector<vector<uint16_t>> sij(2, vector<uint16_t>(dec_param.nb));
//     vector<uint16_t> a_gf = rel_theta ? theta_1.intrinsic_GF : phi_1_p.intrinsic_GF;
//     vector<softdata_t> a = rel_theta ? theta_1.intrinsic_LLR : phi_1_p.intrinsic_LLR;
//     vector<uint16_t> b_gf = rel_theta ? phi_1_p.intrinsic_GF : theta_1.intrinsic_GF;
//     vector<softdata_t> b = rel_theta ? phi_1_p.intrinsic_LLR : theta_1.intrinsic_LLR;
//     uint16_t nH = a_gf.size() < dec_param.nH ? a_gf.size() : dec_param.nH;
//     uint16_t nL = b_gf.size() < dec_param.nL ? b_gf.size() : dec_param.nL;
//     vector<vector<bool>> Ti(nH, vector<bool>(nL, false));
//     for (int i = 0; i < dec_param.nb; ++i)
//     {
//         srtr[i] = a[i];
//         srtrg[i] = ADDGF[a_gf[i]][b_gf[0]];
//         sij[0][i] = i;
//         sij[1][i] = 0;
//         Ti[i][0] = true;
//     }
//     int cnt = 0, nop = 0;
//     while (nop < dec_param.nopM)
//     {
//         ++nop;
//         auto min_it = min_element(srtr.begin(), srtr.end());
//         int n = distance(srtr.begin(), min_it);
//         int i = sij[0][n];
//         int j = sij[1][n];
//         if (find(theta.intrinsic_GF.begin(), theta.intrinsic_GF.end(), srtrg[n]) == theta.intrinsic_GF.end())
//         {
//             theta.intrinsic_LLR.push_back(*min_it);
//             theta.intrinsic_GF.push_back(srtrg[n]);
//             ++cnt;
//         }
//         if (i == nH - 1 || j == nL - 1 || cnt == dec_param.nm)
//             break;
//         int H = (i == 0) ? 1 : (i == dec_param.nb - 1) ? 0
//                                                        : 1;
//         int Hb = 1 - H;
//         int i1 = (!Ti[i + Hb][j + H]) ? i + Hb : i + H;
//         int j1 = (!Ti[i + Hb][j + H]) ? j + H : j + Hb;
//         Ti[i1][j1] = true;
//         srtr[n] = a[i1] + b[j1];
//         srtrg[n] = ADDGF[a_gf[i1]][b_gf[j1]];
//         sij[0][n] = i1;
//         sij[1][n] = j1;
//     }
// }

void PoAwN::decoding::ECN_EMS(const decoder_t &theta_1,
                              const decoder_t &phi_1,
                              const vector<vector<uint16_t>> &ADDGF,
                              const vector<vector<uint16_t>> &DIVGF,
                              const decoder_parameters &dec_param,
                              const uint16_t coef,
                              decoder_t &theta)
{
    decoder_t phi_1_p = phi_1;
    for (int i = 0; i < phi_1.intrinsic_GF.size(); i++)
        phi_1_p.intrinsic_GF[i] = DIVGF[phi_1_p.intrinsic_GF[i]][coef];

    vector<uint16_t> a_gf, b_gf;
    vector<softdata_t> a, b;

    a_gf.assign(theta_1.intrinsic_GF.begin(), theta_1.intrinsic_GF.begin() + dec_param.nH);
    b_gf.assign(phi_1_p.intrinsic_GF.begin(), phi_1_p.intrinsic_GF.begin() + dec_param.nL);
    a.assign(theta_1.intrinsic_LLR.begin(), theta_1.intrinsic_LLR.begin() + dec_param.nH);
    b.assign(phi_1_p.intrinsic_LLR.begin(), phi_1_p.intrinsic_LLR.begin() + dec_param.nL);

    uint16_t nH = dec_param.nH;
    uint16_t nL = dec_param.nL;
    uint16_t N1 = nH * nL;
    decoder_t bubble;
    bubble.intrinsic_GF.resize(nH * nL);
    bubble.intrinsic_LLR.resize(nH * nL);
    vector<vector<uint16_t>> sort_bub_idxs(dec_param.nm, vector<uint16_t>(2));
    vector<vector<uint16_t>> temp_idxs(dec_param.q, vector<uint16_t>(2, 0));

    for (int i = 0; i < nH; i++)
        for (int j = 0; j < nL; ++j)
        {
            bubble.intrinsic_GF[i * nL + j] = ADDGF[a_gf[i]][b_gf[j]];
            bubble.intrinsic_LLR[i * nL + j] = a[i] + b[j];
        }

    vector<softdata_t> temp_llr(dec_param.q, std::numeric_limits<softdata_t>::max() / 2);
    vector<uint16_t> temp_GF(dec_param.q);
    iota(temp_GF.begin(), temp_GF.end(), 0);

    for (int i = 0; i < N1; i++)
        if (bubble.intrinsic_LLR[i] < temp_llr[bubble.intrinsic_GF[i]])
            temp_llr[bubble.intrinsic_GF[i]] = bubble.intrinsic_LLR[i];

    partial_sort(temp_GF.begin(), temp_GF.begin() + dec_param.nm, temp_GF.end(), [&](int i, int j)
                 { return temp_llr[i] < temp_llr[j]; });

    theta.intrinsic_LLR.assign(dec_param.q, temp_llr[temp_GF[dec_param.nm - 1]] + dec_param.offset);
    for (int i = 0; i < dec_param.nm; i++)
    {
        theta.intrinsic_GF[i] = temp_GF[i];
        theta.intrinsic_LLR[i] = temp_llr[temp_GF[i]];
    }

    for (int i = dec_param.nm; i < dec_param.q; i++)
        theta.intrinsic_GF[i] = temp_GF[i];
}

void PoAwN::decoding::ECN_EMS_L(const decoder_t &theta_1,
                                const decoder_t &phi_1,
                                const vector<vector<uint16_t>> &ADDGF,
                                const vector<vector<uint16_t>> &DIVGF,
                                const decoder_parameters &dec_param,
                                const uint16_t coef,
                                const uint16_t ucap_theta,
                                const uint16_t ucap_phi,
                                decoder_t &theta,
                                vector<uint16_t> &Bt1)
{
    decoder_t phi_1_p = phi_1;
    bool rel_theta = (theta_1.intrinsic_LLR[dec_param.Zc] > phi_1.intrinsic_LLR[dec_param.Zc]);
    if (dec_param.sig_mod == "BPSK")
        for (int i = 0; i < phi_1.intrinsic_GF.size(); i++)
            phi_1_p.intrinsic_GF[i] = DIVGF[phi_1_p.intrinsic_GF[i]][coef];

    vector<uint16_t> a_gf, b_gf;
    vector<softdata_t> a, b;
    if (rel_theta)
    {
        a_gf=theta_1.intrinsic_GF;;
        a=theta_1.intrinsic_LLR;
        b_gf=phi_1_p.intrinsic_GF;
        b=phi_1_p.intrinsic_LLR;
    }
    else
    {
        b_gf=theta_1.intrinsic_GF;
        b=theta_1.intrinsic_LLR;
        a_gf=phi_1_p.intrinsic_GF;
        a=phi_1_p.intrinsic_LLR;
    }


    uint16_t N1 = dec_param.nH * dec_param.nL;

    decoder_t bubble;
    bubble.intrinsic_GF.resize(dec_param.nH * dec_param.nL);
    bubble.intrinsic_LLR.resize(dec_param.nH * dec_param.nL);
    vector<vector<uint16_t>> idxs(dec_param.nH * dec_param.nL, vector<uint16_t>(2, 0));

    for (int i = 0; i < dec_param.nH; i++)
        for (int j = 0; j < dec_param.nL; ++j)
        {
            idxs[i * dec_param.nL + j][0] = i;
            idxs[i * dec_param.nL + j][1] = j;
            bubble.intrinsic_GF[i * dec_param.nL + j] = ADDGF[a_gf[i]][b_gf[j]];
            bubble.intrinsic_LLR[i * dec_param.nL + j] = a[i] + b[j];
        }

    vector<softdata_t> temp_llr(dec_param.q, std::numeric_limits<softdata_t>::max() / 2);
    vector<uint16_t> temp_GF(dec_param.q);
    iota(temp_GF.begin(), temp_GF.end(), 0);

    for (int i = 0; i < N1; i++)
        if (bubble.intrinsic_LLR[i] < temp_llr[bubble.intrinsic_GF[i]])
        {
            temp_llr[bubble.intrinsic_GF[i]] = bubble.intrinsic_LLR[i];
        }

    partial_sort(temp_GF.begin(), temp_GF.begin() + dec_param.nm, temp_GF.end(), [&](int i, int j)
                 { return temp_llr[i] < temp_llr[j]; });

    theta.intrinsic_LLR.assign(dec_param.nm, temp_llr[temp_GF[dec_param.nm - 1]] + dec_param.offset);
    for (int i = 0; i < dec_param.nm; i++)
    {
        theta.intrinsic_GF[i] = temp_GF[i];
        theta.intrinsic_LLR[i] = temp_llr[temp_GF[i]];
    }


    bool brk1 = false;
    for (int i = 0; i < dec_param.nH; i++)
    {
        for (int j = 0; j < dec_param.nL; j++)
        {
            if ((rel_theta && a_gf[i] == ucap_theta && b_gf[j] == ucap_phi) ||
                (!rel_theta && a_gf[i] == ucap_phi && b_gf[j] == ucap_theta))
            {
                Bt1[0] = i;
                Bt1[1] = j;
                brk1 = true;
                break;
            }
        }
        if (brk1)
            break;
    }
    int aa=-1;
}

void PoAwN::decoding::ECN_PA(const decoder_t &theta_1,
                             const decoder_t &phi_1,
                             const vector<vector<uint16_t>> &ADDGF,
                             const vector<vector<uint16_t>> &DIVGF,
                             const decoder_parameters &dec_param,
                             const uint16_t coef,
                             const uint16_t l,
                             const uint16_t s,
                             decoder_t &theta,
                             uint16_t &nm)

{
    vector<vector<uint16_t>> Bubb_Indicator = dec_param.Bubb_Indicator[l][s];
    uint16_t nH = dec_param.ns[l][s][0];
    uint16_t nL = dec_param.ns[l][s][1];
    nm = nL;
    // nm = std::max({dec_param.ns[l][s][1], dec_param.ns[l+1][2*s][1], dec_param.ns[l+1][2*s+1][1]});
    decoder_t phi_1_p = phi_1;
    bool rel_theta = (theta_1.intrinsic_LLR[dec_param.Zc] > phi_1.intrinsic_LLR[dec_param.Zc]);
    if (dec_param.sig_mod == "BPSK")
        for (int i = 0; i < phi_1.intrinsic_GF.size(); i++)
            phi_1_p.intrinsic_GF[i] = DIVGF[phi_1_p.intrinsic_GF[i]][coef];
    vector<uint16_t> a_gf, b_gf;
    vector<softdata_t> a, b;
    if (rel_theta)
    {
        a_gf=theta_1.intrinsic_GF;;
        a=theta_1.intrinsic_LLR;
        b_gf=phi_1_p.intrinsic_GF;
        b=phi_1_p.intrinsic_LLR;
    }
    else
    {
        b_gf=theta_1.intrinsic_GF;
        b=theta_1.intrinsic_LLR;
        a_gf=phi_1_p.intrinsic_GF;
        a=phi_1_p.intrinsic_LLR;
    }

    uint16_t sz_bub = Bubb_Indicator[0].size();
    decoder_t bubbles;
    bubbles.intrinsic_GF.resize(sz_bub);
    bubbles.intrinsic_LLR.resize(sz_bub);

    uint16_t i1, i2;
    for (int i = 0; i < sz_bub; i++)
    {
        i1 = Bubb_Indicator[0][i];
        i2 = Bubb_Indicator[1][i];
        bubbles.intrinsic_GF[i] = ADDGF[a_gf[i1]][b_gf[i2]];
        bubbles.intrinsic_LLR[i] = a[i1] + b[i2];
    }

    uint16_t q = dec_param.q;

    vector<softdata_t> temp_llr(q);
    vector<uint16_t> temp_GF(q);

    temp_llr.assign(q, std::numeric_limits<softdata_t>::max() / 2);
    temp_GF.resize(q);
    iota(temp_GF.begin(), temp_GF.end(), 0);

    for (int i = 0; i < sz_bub; i++)
        if (bubbles.intrinsic_LLR[i] < temp_llr[bubbles.intrinsic_GF[i]])
            temp_llr[bubbles.intrinsic_GF[i]] = bubbles.intrinsic_LLR[i];

    partial_sort(temp_GF.begin(), temp_GF.begin() + nm, temp_GF.end(), [&](int i, int j)
                 { return temp_llr[i] < temp_llr[j]; });

    theta.intrinsic_LLR.assign(nm, temp_llr[temp_GF[nm - 1]] + dec_param.offset);
    for (int i = 0; i < nm; i++)
    {
        theta.intrinsic_GF[i] = temp_GF[i];
        theta.intrinsic_LLR[i] = temp_llr[temp_GF[i]];
    }
}

void PoAwN::decoding::LLR_sort(const vector<vector<softdata_t>> &chan_LLR,
                               const uint16_t nm,
                               vector<decoder_t> &chan_LLR_sorted,
                               const softdata_t ofst)
{
    int N = chan_LLR.size();
    int q = chan_LLR[0].size();
    vector<uint16_t> temp_GF(q);
    vector<softdata_t> temp_llr_arr(q);
    vector<uint16_t> temp_gf_arr(q);
    for (int k = 0; k < N; ++k)
    {
        temp_llr_arr = chan_LLR[k];
        iota(temp_GF.begin(), temp_GF.end(), 0);
        partial_sort(temp_GF.begin(), temp_GF.begin() + nm, temp_GF.end(),
                     [&temp_llr_arr](int i1, int i2)
                     { return temp_llr_arr[i1] < temp_llr_arr[i2]; });
        for (int i = 0; i < nm; ++i)
        {
            chan_LLR_sorted[k].intrinsic_LLR[i] = temp_llr_arr[temp_GF[i]];
            chan_LLR_sorted[k].intrinsic_GF[i] = temp_GF[i];
        }

    }
}

// void PoAwN::decoding::ECN_AEMS_bubble(const decoder_t &theta_1,
//                                       const decoder_t &phi_1,
//                                       const decoder_parameters &dec_param,
//                                       uint16_t coef,
//                                       const vector<vector<uint16_t>> &ADDGF,
//                                       const vector<vector<uint16_t>> &DIVGF,
//                                       decoder_t &theta)
// {
//     theta.intrinsic_LLR.reserve(dec_param.nm);
//     theta.intrinsic_GF.reserve(dec_param.nm);
//     bool rel_theta = (theta_1.intrinsic_LLR[dec_param.Zc] > phi_1.intrinsic_LLR[dec_param.Zc]);
//     decoder_t phi_1_p = phi_1;
//     for (int i = 0; i < phi_1.intrinsic_GF.size(); i++)
//         phi_1_p.intrinsic_GF[i] = DIVGF[phi_1_p.intrinsic_GF[i]][coef];
//     vector<softdata_t> srtr(dec_param.nb);
//     vector<uint16_t> srtrg(dec_param.nb);
//     vector<vector<uint16_t>> sij(2, vector<uint16_t>(dec_param.nb));
//     vector<uint16_t> a_gf = rel_theta ? theta_1.intrinsic_GF : phi_1_p.intrinsic_GF;
//     vector<softdata_t> a = rel_theta ? theta_1.intrinsic_LLR : phi_1_p.intrinsic_LLR;
//     vector<uint16_t> b_gf = rel_theta ? phi_1_p.intrinsic_GF : theta_1.intrinsic_GF;
//     vector<softdata_t> b = rel_theta ? phi_1_p.intrinsic_LLR : theta_1.intrinsic_LLR;
//     uint16_t nH = a_gf.size() < dec_param.nH ? a_gf.size() : dec_param.nH;
//     uint16_t nL = b_gf.size() < dec_param.nL ? b_gf.size() : dec_param.nL;
//     vector<vector<bool>> Ti(nH, vector<bool>(nL, false));
//     for (int i = 0; i < dec_param.nb; ++i)
//     {
//         srtr[i] = a[i];
//         srtrg[i] = ADDGF[a_gf[i]][b_gf[0]];
//         sij[0][i] = i;
//         sij[1][i] = 0;
//         Ti[i][0] = true;
//     }
//     int cnt = 0, nop = 0;
//     while (nop < dec_param.nopM)
//     {
//         ++nop;
//         auto min_it = min_element(srtr.begin(), srtr.end());
//         int n = distance(srtr.begin(), min_it);
//         int i = sij[0][n];
//         int j = sij[1][n];
//         if (find(theta.intrinsic_GF.begin(), theta.intrinsic_GF.end(), srtrg[n]) == theta.intrinsic_GF.end())
//         {
//             theta.intrinsic_LLR.push_back(*min_it);
//             theta.intrinsic_GF.push_back(srtrg[n]);
//             ++cnt;
//         }
//         if (i == nH - 1 || j == nL - 1 || cnt == dec_param.nm)
//             break;
//         int H = (i == 0) ? 1 : (i == dec_param.nb - 1) ? 0
//                                                        : 1;
//         int Hb = 1 - H;
//         int i1 = (!Ti[i + Hb][j + H]) ? i + Hb : i + H;
//         int j1 = (!Ti[i + Hb][j + H]) ? j + H : j + Hb;
//         Ti[i1][j1] = true;
//         srtr[n] = a[i1] + b[j1];
//         srtrg[n] = ADDGF[a_gf[i1]][b_gf[j1]];
//         sij[0][n] = i1;
//         sij[1][n] = j1;
//     }
// }

void PoAwN::decoding::VN_update(const decoder_t &theta_1,
                                const decoder_t &phi_1,
                                const vector<vector<uint16_t>> &ADDGF,
                                const vector<vector<uint16_t>> &DIVGF,
                                const decoder_parameters &dec_param,
                                uint16_t coef,
                                uint16_t hard_decision,
                                decoder_t &phi)
{
    uint16_t nm = dec_param.nm, q = dec_param.q;
    vector<softdata_t> theta1_llr(q, theta_1.intrinsic_LLR[nm-1]+dec_param.offset);
    vector<softdata_t> phi1_llr(q, phi_1.intrinsic_LLR[nm-1]+dec_param.offset);
    vector<softdata_t> temp_llr(q);
    for (int i = 0; i < nm; i++)
    {
        phi1_llr[phi_1.intrinsic_GF[i]] = phi_1.intrinsic_LLR[i];
        theta1_llr[ADDGF[hard_decision][theta_1.intrinsic_GF[i]]] = theta_1.intrinsic_LLR[i];
    }

    softdata_t mn_llr = std::numeric_limits<softdata_t>::max();

    if (dec_param.sig_mod == "BPSK")
        for (int i = 0; i < q; i++)
        {
            temp_llr[DIVGF[i][coef]] = theta1_llr[i] + phi1_llr[i];
            if (temp_llr[DIVGF[i][coef]] < mn_llr)
                mn_llr = temp_llr[DIVGF[i][coef]];
        }
    else
        for (int i = 0; i < q; i++)
        {
            temp_llr[i] = theta1_llr[i] + phi1_llr[i];
            if (temp_llr[i] < mn_llr)
                mn_llr = temp_llr[i];
        }

    if (std::abs(mn_llr > 1e-5))
        for (int i = 0; i < q; i++)
            temp_llr[i] -= mn_llr;

    vector<uint16_t> temp_GF(q);
    iota(temp_GF.begin(), temp_GF.end(), 0);
    partial_sort(temp_GF.begin(), temp_GF.begin() + nm, temp_GF.end(),
                 [&temp_llr](int i1, int i2)
                 { return temp_llr[i1] < temp_llr[i2]; });

    for (int i = 0; i < nm; ++i)
    {
        phi.intrinsic_LLR[i] = temp_llr[temp_GF[i]];
        phi.intrinsic_GF[i] = temp_GF[i];
    }
}

void PoAwN::decoding::VN_update_PA(const decoder_t &theta_1,
                                   const decoder_t &phi_1,
                                   const vector<vector<uint16_t>> &ADDGF,
                                   const vector<vector<uint16_t>> &DIVGF,
                                   const decoder_parameters &dec_param,
                                   uint16_t coef,
                                   const uint16_t l,
                                   const uint16_t s,
                                   const uint16_t df_nm,
                                   uint16_t hard_decision,
                                   decoder_t &phi)
{

    uint16_t nm = 0;
    uint16_t q = dec_param.q;
    // if (l < dec_param.n - 1)
    //     nm = dec_param.ns[l + 1][2 * s + 1][1];
    // if (nm == 0)
    //     nm = dec_param.ns[l][s][1];
    // if (nm == 0)
    //     nm = df_nm;
    if (l == 0)
        nm = dec_param.ns[l][s][1];
    else
        nm = dec_param.ns[l - 1][(int)floor((double)s / 2)][1];
    vector<softdata_t> theta1_llr(q, theta_1.intrinsic_LLR[nm-1]);
    vector<softdata_t> phi1_llr(q, phi_1.intrinsic_LLR[nm-1]);
    vector<softdata_t> temp_llr(q);
    for (int i = 0; i < q; i++)
    {
        phi1_llr[phi_1.intrinsic_GF[i]] = phi_1.intrinsic_LLR[i];
        theta1_llr[ADDGF[hard_decision][theta_1.intrinsic_GF[i]]] = theta_1.intrinsic_LLR[i];
    }

    softdata_t mn_llr = std::numeric_limits<softdata_t>::max();

    if (dec_param.sig_mod == "BPSK")
        for (int i = 0; i < q; i++)
        {
            temp_llr[DIVGF[i][coef]] = theta1_llr[i] + phi1_llr[i];
            if (temp_llr[DIVGF[i][coef]] < mn_llr)
                mn_llr = temp_llr[DIVGF[i][coef]];
        }
    else
        for (int i = 0; i < q; i++)
        {
            temp_llr[i] = theta1_llr[i] + phi1_llr[i];
            if (temp_llr[i] < mn_llr)
                mn_llr = temp_llr[i];
        }

    if (std::abs(mn_llr > 1e-5))
        for (int i = 0; i < q; i++)
            temp_llr[i] -= mn_llr;

    vector<uint16_t> temp_GF(q);
    iota(temp_GF.begin(), temp_GF.end(), 0);
    partial_sort(temp_GF.begin(), temp_GF.begin() + nm, temp_GF.end(),
                 [&temp_llr](int i1, int i2)
                 { return temp_llr[i1] < temp_llr[i2]; });

    for (int i = 0; i < nm; ++i)
    {
        phi.intrinsic_LLR[i] = temp_llr[temp_GF[i]];
        phi.intrinsic_GF[i] = temp_GF[i];
    }
}

void PoAwN::decoding::decode_SC_PA(const decoder_parameters &dec_param,
                                   const vector<vector<uint16_t>> &ADDGF,
                                   const vector<vector<uint16_t>> &MULGF,
                                   const vector<vector<uint16_t>> &DIVGF,
                                   vector<vector<decoder_t>> &L,
                                   vector<uint16_t> &info_sec_rec)
{
    vector<vector<bool>> Roots = dec_param.Roots_V;
    uint16_t MxUS = dec_param.MxUS, n = dec_param.n, N = dec_param.N;
    int l = 0, s = 0;
    uint16_t hard_decsion, temp_coef, i1, i2, i3, SZc, SZc1, l1;
    vector<uint16_t> Root;

    // vector<vector<uint16_t>> V(n + 1, vector<uint16_t>(dec_param.N, dec_param.MxUS));
    // V[n] = (vector<uint16_t>(dec_param.N, dec_param.frozen_val));
    // for (uint16_t i = 0; i < dec_param.K; i++)
    //     V[n][dec_param.reliab_sequence[i]] = dec_param.MxUS;
    vector<vector<uint16_t>> V = dec_param.ufrozen;
    uint16_t temp_GF_p, theta_gf;
    decoder_t theta_t({0}, {0});
    uint16_t nm = dec_param.MxUS;
    while (l > -1)
    {
        if (Roots[l][s])
        {
            if (s % 2 == 1)
            {
                l -= 1;
                s = (s - 1) / 2;
                Root = dec_param.Roots_indices[l][s];
                SZc = Root.size();
                SZc1 = SZc >> 1;
                for (uint16_t t = 0; t < SZc1; t++)
                {
                    l1 = l + 1;
                    i1 = Root[t], i2 = Root[t + SZc1], i3 = dec_param.coefs_id[l][s][t];
                    temp_coef = dec_param.polar_coeff[n - l - 1][i3];
                    V[l][i1] = ADDGF[V[l1][i1]][V[l1][i2]];
                    V[l][i2] = MULGF[V[l1][i2]][temp_coef];
                }
                Roots[l][s] = true;
            }
            else
            {
                l -= 1;
                s /= 2;
            }
        }

        else if (Roots[l + 1][2 * s])
        {
            Root = dec_param.Roots_indices[l][s];
            SZc = Root.size();
            SZc1 = SZc >> 1;
            for (uint16_t t = 0; t < SZc1; t++)
            {
                i3 = dec_param.coefs_id[l][s][t];
                temp_coef = dec_param.polar_coeff[n - l - 1][i3];
                hard_decsion = V[l + 1][Root[t]];
                VN_update_PA(L[l][Root[t]], L[l][Root[t + SZc1]], ADDGF, DIVGF, dec_param,
                             temp_coef, l, s, nm, hard_decsion, L[l + 1][Root[t + SZc1]]);
            }
            l += 1;
            s = 2 * s + 1;
            if (l == n)
            {
                Roots[n][s] = true;
                if (V[n][s] == MxUS)
                    V[n][s] = L[n][s].intrinsic_GF[0];
                if (s == N - 1)
                    break;
            }
        }

        else
        {
            Root = dec_param.Roots_indices[l][s];
            SZc = Root.size();
            SZc1 = SZc >> 1;
            decoder_t phi_1, theta_1;
            int cnt0 = 0;

            for (uint16_t t = 0; t < SZc1; t++)
            {
                if (dec_param.ufrozen[l + 1][Root[t]] != dec_param.frozen_val)
                {
                    i3 = dec_param.coefs_id[l][s][t];
                    temp_coef = dec_param.polar_coeff[n - l - 1][i3];
                    if (l < n - 1)
                    {
                        ECN_PA(L[l][Root[t]], L[l][Root[t + SZc1]], ADDGF, DIVGF, dec_param,
                               temp_coef, l, s, L[l + 1][Root[t]], nm);
                    }
                    else
                    {
                        temp_GF_p = DIVGF[L[l][Root[t + SZc1]].intrinsic_GF[0]][temp_coef];
                        theta_gf = L[l][Root[t]].intrinsic_GF[0];
                        theta_t.intrinsic_GF[0] = ADDGF[theta_gf][temp_GF_p];
                        L[l + 1][Root[t]] = theta_t;
                    }
                }
                else
                {
                    V[l + 1][Root[t]] = dec_param.frozen_val;
                    cnt0++;
                }
            }

            l = l + 1;
            s = 2 * s;
            if (cnt0 == SZc1)
            {
                Roots[l][s] = true;
            }
            if (l == n)
            {
                Roots[n][s] = true;
                if (V[n][s] == MxUS)
                    V[n][s] = L[n][s].intrinsic_GF[0];
            }
        }
    }
    for (uint16_t i = 0; i < dec_param.K; i++)
        info_sec_rec[i] = V[n][dec_param.reliab_sequence[i]];
}

void PoAwN::decoding::decode_SC_bubble_gen(const decoder_parameters &dec_param,
                                           const vector<vector<uint16_t>> &ADDGF,
                                           const vector<vector<uint16_t>> &MULGF,
                                           const vector<vector<uint16_t>> &DIVGF,
                                           vector<vector<decoder_t>> &L,
                                           vector<uint16_t> &info_sec_rec,
                                           vector<vector<vector<vector<uint16_t>>>> &Bt)
{
    uint16_t MxUS = dec_param.MxUS, n = dec_param.n, N = dec_param.N;
    vector<vector<bool>> Roots(n + 1);
    for (int i = 0; i < n; i++)
        Roots[i] = dec_param.clst_frozen[i];

    Roots[n].assign(N, false);
    for (int i = dec_param.K; i < N; i++)
    {
        Roots[n][dec_param.reliab_sequence[i]] = true;
    }

    int l = 0, s = 0;
    uint16_t hard_decsion, temp_coef, i1, i2, i3, SZc, SZc1, l1;
    vector<uint16_t> Root;

    vector<vector<uint16_t>> V = dec_param.ufrozen;
    bool c1;
    while (l > -1)
    {
        if (Roots[l][s])
        {
            if (s % 2 == 1)
            {
                l -= 1;
                s = (s - 1) / 2;
                Root = dec_param.Roots_indices[l][s];
                SZc = Root.size();
                SZc1 = SZc >> 1;

                for (uint16_t t = 0; t < SZc1; t++)
                {
                    l1 = l + 1;
                    i1 = Root[t], i2 = Root[t + SZc1], i3 = dec_param.coefs_id[l][s][t];
                    temp_coef = dec_param.polar_coeff[n - l - 1][i3];
                    V[l][i1] = ADDGF[V[l1][i1]][V[l1][i2]];
                    V[l][i2] = MULGF[V[l1][i2]][temp_coef];
                }
                Roots[l][s] = true;
            }
            else
            {
                l -= 1;
                s /= 2;
            }
        }
        else if (Roots[l + 1][2 * s])
        {
            Root = dec_param.Roots_indices[l][s];
            SZc = Root.size();
            SZc1 = SZc >> 1;
            for (uint16_t t = 0; t < SZc1; t++)
            {
                i3 = dec_param.coefs_id[l][s][t];
                temp_coef = dec_param.polar_coeff[n - l - 1][i3];
                hard_decsion = V[l + 1][Root[t]];
                bool cnd1 = hard_decsion != dec_param.ucap[l + 1][Root[t]];
                VN_update(L[l][Root[t]], L[l][Root[t + SZc1]], ADDGF, DIVGF, dec_param,
                          temp_coef, hard_decsion, L[l + 1][Root[t + SZc1]]);
                          int aa=1;
                // for (int y = 0; y < dec_param.nm; y++)
                // {
                //     c1 = (L[l + 1][Root[t + SZc1]].intrinsic_GF[y] == dec_param.ucap[l + 1][Root[t + SZc1]]);
                //     if (c1)
                //     {
                //         dec_param.cnd1[l][Root[t + SZc1]] = y;
                //         break;
                //     }
                // }
            }
            l += 1;
            s = 2 * s + 1;
            if (l == n)
            {
                Roots[n][s] = true;
                if (V[n][s] == MxUS)
                    V[n][s] = L[n][s].intrinsic_GF[0];
                if (s == N - 1)
                    break;
            }
        }
        else
        {
            Root = dec_param.Roots_indices[l][s];
            SZc = Root.size();
            SZc1 = SZc >> 1;

            for (uint16_t t = 0; t < SZc1; t++)
            {

                i3 = dec_param.coefs_id[l][s][t];
                temp_coef = dec_param.polar_coeff[n - l - 1][i3];
                // if (l < n - 1)
                ECN_EMS_L(L[l][Root[t]], L[l][Root[t + SZc1]], ADDGF, DIVGF, dec_param, temp_coef,
                          dec_param.ucap[l][Root[t]], dec_param.ucap[l][Root[t + SZc1]], L[l + 1][Root[t]], Bt[l][s][t]);
                          int aa=-1;

            }
            l = l + 1;
            s = 2 * s;
            if (l == n)
            {
                Roots[n][s] = true;
                if (V[n][s] == MxUS)
                    V[n][s] = L[n][s].intrinsic_GF[0];
            }
        }
    }
    info_sec_rec.resize(dec_param.K, dec_param.MxUS);
    for (uint16_t i = 0; i < dec_param.K; i++)
    {
        info_sec_rec[i] = V[n][dec_param.reliab_sequence[i]];
    }
}
void PoAwN::decoding::frozen_lay_pos(const decoder_parameters &dec_param,
                                     vector<vector<uint16_t>> &ufrozen,
                                     vector<vector<bool>> &clst_frozen)
{
    clst_frozen.resize(dec_param.n);
    ufrozen.resize(dec_param.n + 1, vector<uint16_t>(dec_param.N, dec_param.MxUS));
    for (int i = dec_param.N - 1; i >= dec_param.K; i--)
        ufrozen[dec_param.n][dec_param.reliab_sequence[i]] = dec_param.frozen_val;
    uint16_t i1, i2;
    for (int l = dec_param.n; l > 0; l--)
    {
        clst_frozen[l - 1].resize(1 << (l - 1), false);
        int sz1 = dec_param.clusts_CNs[l - 1].size();
        for (int s = 0; s < sz1; s++)
        {
            int sz2 = dec_param.clusts_CNs[l - 1][s].size();
            int cnt2 = 0;
            for (int k = 0; k < sz2; k++)
            {
                i1 = dec_param.clusts_CNs[l - 1][s][k];
                i2 = dec_param.clusts_VNs[l - 1][s][k];
                if (ufrozen[l][i1] == dec_param.frozen_val && ufrozen[l][i2] == dec_param.frozen_val)
                {
                    ufrozen[l - 1][i1] = dec_param.frozen_val;
                    ufrozen[l - 1][i2] = dec_param.frozen_val;
                    cnt2++;
                }
            }
            if (cnt2 == sz2)
                clst_frozen[l - 1][s] = true;
        }
    }
}
