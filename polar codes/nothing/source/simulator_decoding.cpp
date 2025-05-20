#include <cmath>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "Decoder_functions.h"
#include "GF_tools.h"
#include "init.h"
#include "struct.h"
#include "tools.h"
#include "HelperFunc.h"
#include <fstream>
#include <iomanip>
#include <cstring>
#include <string>
#include <iomanip>
#include <algorithm>
#include "channel.h"
#include <filesystem>
#include <omp.h>
#include <atomic>

using namespace PoAwN::structures;
using namespace PoAwN::tools;
using namespace PoAwN::init;
using namespace PoAwN::decoding;
using namespace PoAwN::channel;
using std::array;
using std::cout;
using std::endl;
using std::stod;
using std::stoi;
using std::string;
using std::vector;

int main(int argc, char *argv[])
{

    if (argc != 10)
    {
        cout << "validate: NbMonteCarlo, SNR, sig_mod(BPSK, CCSK_bin, CCSK_NB), q, N, K, nH, nL, offset," << std::endl;
        return 1;
    }
    uint16_t q, N, K, n, nL, nH, nm, nb, Zc, nopM, p, frozen_val = 0;
    softdata_t offset;
    int NbMonteCarlo = stoi(argv[1]);
    float EbN0 = stod(argv[2]);
    string sig_mod = argv[3];
    std::transform(sig_mod.begin(), sig_mod.end(), sig_mod.begin(), ::toupper);
    q = stoi(argv[4]);
    p = log2(q);
    N = stoi(argv[5]);
    K = stoi(argv[6]);
    n = log2(N);
    nH = stoi(argv[7]);
    nL = stoi(argv[8]);
    nm = nL;
    Zc = 2;
    offset = stod(argv[9]);
    int FER_STOP=100;

    base_code_t code_param(N, K, n, q, p, frozen_val);
    code_param.sig_mod = sig_mod;

    int gf_rand_SEED = 0;
    float nse_rand_SEED = 1.2544;
    bool repeatable_randgen = 0;

    table_GF table;

    cout << "Loading code_param..." << endl;
    LoadCode(code_param, EbN0);

    cout << "OK!, " << "Loading tables..." << endl;
    // void LoadTables(base_code_t & code, table_GF & table,  const uint16_t *GF_polynom_primitive)
    LoadTables(code_param, table, GF_polynom_primitive.data());
    cout << "Done!" << endl;
    cout << "Simulation starts..." << endl;

    decoder_parameters dec_param(code_param, offset, nm, nL, nH, nb, Zc, nopM);
    LoadBubblesIndcatorlists(dec_param, EbN0);

    vector<vector<std::array<int, 2>>> ns0(n);
    dec_param.ns.resize(n);

    for (int i = 0; i < n; i++)
    {
        ns0[i].resize(1 << i);
        dec_param.ns[i].resize(1 << i);
        for (int j = 0; j < 1 << i; j++)
        {
            dec_param.ns[i][j].resize(2);
            for (int k = 0; k < dec_param.Bubb_Indicator[i][j][0].size(); k++)
            {
                if (dec_param.Bubb_Indicator[i][j][0][k] + 1 > ns0[i][j][0])
                    ns0[i][j][0] = dec_param.Bubb_Indicator[i][j][0][k] + 1;
                if (dec_param.Bubb_Indicator[i][j][1][k] + 1 > ns0[i][j][1])
                    ns0[i][j][1] = dec_param.Bubb_Indicator[i][j][1][k] + 1;
            }
            dec_param.ns[i][j][0] = ns0[i][j][0];
            dec_param.ns[i][j][1] = ns0[i][j][1];
        }
    }

    dec_param.Roots_V.resize(n + 1);
    dec_param.Roots_indices.resize(n);
    dec_param.clusts_CNs.resize(n);
    dec_param.clusts_VNs.resize(n);
    dec_param.coefs_id.resize(n);

    dec_param.Roots_V[n].resize(1U << n, false);
    for (uint16_t l = 0; l < n; l++)
    {
        dec_param.Roots_V[l].resize(1U << l, false);
        dec_param.Roots_indices[l].resize(pow(2, l));
        dec_param.clusts_CNs[l].resize(pow(2, l));
        dec_param.clusts_VNs[l].resize(pow(2, l));
        dec_param.coefs_id[l].resize(pow(2, l));
        for (uint16_t s = 0; s < dec_param.Roots_V[l].size(); s++)
        {
            uint16_t sz1 = N >> (l + 1U), sz2 = sz1 << 1U;
            dec_param.clusts_CNs[l][s].resize(sz1);
            dec_param.clusts_VNs[l][s].resize(sz1);
            dec_param.coefs_id[l][s].resize(sz1);
            for (uint16_t t = 0; t < sz1; ++t)
            {
                dec_param.clusts_CNs[l][s][t] = s * sz2 + (t << 1U) - (t % sz1);
                dec_param.clusts_VNs[l][s][t] = dec_param.clusts_CNs[l][s][t] + (N >> (l + 1U));
                dec_param.coefs_id[l][s][t] = dec_param.clusts_VNs[l][s][t] - (s + 1) * sz1;
            }
        }
        for (uint16_t s = 0; s < dec_param.Roots_indices[l].size(); s++)
        {

            uint16_t sz1 = N >> l;
            dec_param.Roots_indices[l][s].resize(sz1);
            for (uint16_t t = 0; t < sz1; ++t)
                dec_param.Roots_indices[l][s][t] = s * sz1 + t;
        }
    }
    frozen_lay_pos(dec_param, dec_param.ufrozen, dec_param.clst_frozen);

    CCSK_seq ccsk_seq;
    vector<vector<uint16_t>> CCSK_rotated_codes(q, vector<uint16_t>());
    if (code_param.sig_mod == "CCSK_BIN")
        create_ccsk_rotated_table(ccsk_seq.CCSK_bin_seq[code_param.p - 2], ccsk_seq.CCSK_bin_seq[code_param.p - 2].size(), CCSK_rotated_codes);
    else if (code_param.sig_mod == "CCSK_NB")
        create_ccsk_rotated_table(ccsk_seq.CCSK_GF_seq[code_param.p - 2], ccsk_seq.CCSK_GF_seq[code_param.p - 2].size(), CCSK_rotated_codes);

    decoder_t temp_dec;
    temp_dec.intrinsic_LLR.reserve(dec_param.nm);
    temp_dec.intrinsic_GF.reserve(dec_param.nm);

    vector<vector<decoder_t>> L(n + 1, vector<decoder_t>(N));
    vector<uint16_t> info_sec_rec(K, dec_param.MxUS);
    vector<uint16_t> KSYMB(K);

    dec_param.ucap.resize(n + 1, vector<uint16_t>(N, dec_param.MxUS));
    dec_param.ucap[n].assign(N, dec_param.frozen_val);

    vector<vector<vector<int16_t>>> hst1(n, vector<vector<int16_t>>(N, vector<int16_t>(dec_param.nm, 0)));

    q = code_param.q;
    p = code_param.p;
    vector<vector<softdata_t>> bin_mod_dict;
    if (code_param.sig_mod == "CCSK_BIN")
    {
        bin_mod_dict.resize(q, vector<softdata_t>(q, 0));

        for (int i = 0; i < q; i++)
            for (int j = 0; j < q; j++)
                bin_mod_dict[i][j] = (CCSK_rotated_codes[i][j] == 0) ? 1 : -1;
    }
    else if (code_param.sig_mod == "BPSK")
    {
        bin_mod_dict.resize(q, vector<softdata_t>(p, 0));
        for (int i = 0; i < q; i++)
            for (int j = 0; j < p; j++)
                bin_mod_dict[i][j] = (table.BINDEC[i][j] == 0) ? 1 : -1;
    }
    // unsigned int FER = 0, i0;
    // for (i0 = 1; i0 <= NbMonteCarlo; ++i0)
    // {
    //     for (int i = 0; i <= n; i++)
    //         for (int j = 0; j < N; j++)
    //             L[i][j] = decoder_t(vector<softdata_t>(q), vector<uint16_t>(q));
    //     if (code_param.sig_mod == "CCSK_BIN")
    //         EncodeChanBPSK_BinCCSK(dec_param, table, EbN0, CCSK_rotated_codes, L[0], KSYMB, bin_mod_dict);
    //     else if (code_param.sig_mod == "CCSK_NB")
    //         EncodeChanGF_CCSK(dec_param, table, EbN0, CCSK_rotated_codes, L[0], KSYMB);
    //     else
    //         EncodeChanBPSK_BinCCSK(dec_param, table, EbN0, table.BINDEC, L[0], KSYMB, bin_mod_dict);

    //     decode_SC_PA(dec_param, table.ADDGF, table.MULGF, table.DIVGF, L, info_sec_rec);

    //     for (uint16_t i = 0; i < dec_param.K; i++)
    //     {
    //         if (KSYMB[i] != info_sec_rec[i])
    //         {
    //             FER++;
    //             break;
    //         }
    //     }
    //     if (FER >= 100)
    //         break;
    //     if ((i0 % 100 == 0 && i0 > 0))
    //         cout << "\rSNR: " << EbN0 << " dB, FER = " << FER << "/" << (float)i0 << " = " << (float)FER / (float)i0 << std::flush;
    // }
    // i0--;
    // cout << "\rSNR: " << EbN0 << " dB, FER = " << FER << "/" << (float)i0 << " = " << (float)FER / (float)i0 << std::flush;

    // cout << endl;

    unsigned int FER_out = 0, gen_frames_out = 0;
    std::atomic<int> global_counter(0);
    std::atomic<int> FER(0);
    std::atomic<bool> stop(false);
    unsigned base_seed = std::chrono::system_clock::now().time_since_epoch().count();
    // const int base_seed = 42;

#pragma omp parallel
    {
        PoAwN::structures::decoder_parameters dec_param_local = dec_param;
        int thread_id = omp_get_thread_num();
        std::mt19937 gen(thread_id + base_seed);
        std::vector<uint16_t> KSYMB(K);
        std::vector<uint16_t> info_sec_rec(K, dec_param_local.MxUS);
        std::vector<std::vector<decoder_t>> L(n + 1, std::vector<decoder_t>(N));
        for (int i = 0; i <= n; i++)
        for (int j = 0; j < N; j++)
            L[i][j] = decoder_t(std::vector<softdata_t>(nm), std::vector<uint16_t>(nm));

        while (!stop.load())
        {
            bool succ_dec = true;
            if (code_param.sig_mod == "CCSK_BIN")
                EncodeChanBPSK_BinCCSK(gen, dec_param_local, table, EbN0, CCSK_rotated_codes, L[0], KSYMB, bin_mod_dict);
            else if (code_param.sig_mod == "CCSK_NB")
                EncodeChanGF_CCSK(gen, dec_param_local, table, EbN0, CCSK_rotated_codes, L[0], KSYMB);
            else
                EncodeChanBPSK_BinCCSK(gen, dec_param_local, table, EbN0, table.BINDEC, L[0], KSYMB, bin_mod_dict);

            decode_SC_PA(dec_param_local, table.ADDGF, table.MULGF, table.DIVGF, L, info_sec_rec);

            for (uint16_t i = 0; i < dec_param_local.K; i++)
            {
                if (KSYMB[i] != info_sec_rec[i])
                {
                    succ_dec = false;
                    break;
                }
            }
            if (global_counter.load() >= NbMonteCarlo || FER.load() >= FER_STOP)
            {
                stop.store(true);
            }
            if (stop.load())
                break;

            global_counter.fetch_add(1);
            if (!succ_dec)
                FER.fetch_add(1);
            if ((global_counter.load() % 100) == 0)
            {
#pragma omp critical
                {
                    FER_out = FER.load();
                    gen_frames_out = global_counter.load();
                    std::cout << "\rSNR: " << EbN0 << " dB, FER = " << FER_out
                              << "/" << gen_frames_out << " = "
                              << (float)FER_out / gen_frames_out << std::flush;
                }
            }
        }
    }

    FER_out = FER.load();
    gen_frames_out = global_counter.load();

    std::cout << "\rSNR: " << EbN0 << " dB, FER = " << FER_out
              << "/" << gen_frames_out << " = "
              << (float)FER_out / gen_frames_out << std::flush << std::endl;
}

