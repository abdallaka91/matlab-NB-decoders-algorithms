#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "Decoder_functions.h"
#include "GF_tools.h"
#include "init.h"
#include "struct.h"
#include "tools.h"
#include "channel.h"
#include "HelperFunc.h"
#include <random>

#define PI 3.14159265358979323846
#define RepRndGn false

using namespace PoAwN::structures;
using namespace PoAwN::tools;
using namespace PoAwN::init;
using namespace PoAwN::decoding;
using std::array;
using std::cout;
using std::endl;
using std::stod;
using std::stoi;
using std::string;
using std::vector;

void PoAwN::channel::EncodeChanBPSK_BinCCSK(std::mt19937 &gen,
                                            decoder_parameters &dec_param,
                                            const table_GF &table,
                                            const float SNR,
                                            const vector<vector<uint16_t>> &bin_table,
                                            vector<decoder_t> &chan_LLR_sorted,
                                            vector<uint16_t> &KSYMB,
                                            const vector<vector<softdata_t>> &bin_mod_dict)
{
    uint16_t N = dec_param.N, K = dec_param.K, q = dec_param.q;
    uint16_t nm = dec_param.nm;
    float sigma = sqrt(1.0 / (pow(10, SNR / 10.0))); // N0/2 or N0?
    vector<uint16_t> NSYMB(N);
    // RandomSymbGenerator(K, q, RepRndGn, 0, KSYMB);
    // std::mt19937 gen(0);
    std::uniform_int_distribution<int> unif_dist(0, q - 1);
    for (uint16_t k = 0; k < K; k++)
    {
        KSYMB[k] = (uint16_t)unif_dist(gen);
    }
    for (int i = 0; i < K; i++)
        dec_param.ucap[dec_param.n][dec_param.reliab_sequence[i]] = KSYMB[i];
    Encoder(table.ADDGF, table.MULGF, dec_param.polar_coeff, dec_param.ucap, NSYMB);

    vector<vector<softdata_t>> noisy_sig(N, vector<softdata_t>(bin_table[0].size(), (softdata_t)0.0));

    for (int i = 0; i < int(noisy_sig.size()); i++)
        noisy_sig[i] = bin_mod_dict[NSYMB[i]];
    // awgn_channel_noise(sigma, RepRndGn, 0, noisy_sig);
    {
        uint16_t q1 = noisy_sig[0].size();
        vector<vector<softdata_t>> noise_table(N, vector<softdata_t>(q1, 0));
        std::normal_distribution<double> norm_dist(0, sigma);
        {
            for (int i = 0; i < N; i++)
            {

                for (int j = 0; j < q1; j++)
                {
                    noise_table[i][j] = (softdata_t)norm_dist(gen);
                }
            }
        }

        for (int i = 0; i < N; i++)
            for (int j = 0; j < q1; j++)
                noisy_sig[i][j] += noise_table[i][j];
    }
    vector<vector<softdata_t>> chan_LLR(N, vector<softdata_t>(q, 0));
    Channel_LLR(noisy_sig, bin_table, q, sigma, chan_LLR);
    LLR_sort(chan_LLR, nm, chan_LLR_sorted, dec_param.offset);
}

void PoAwN::channel::EncodeChanGF_CCSK(std::mt19937 &gen,
                                       decoder_parameters &dec_param,
                                       const table_GF &table,
                                       const float SNR,
                                       const vector<vector<uint16_t>> &CCSK_rotated_codes,
                                       vector<decoder_t> &chan_LLR_sorted,
                                       vector<uint16_t> &KSYMB)
{
    uint16_t N = dec_param.N, K = dec_param.K, q = dec_param.q, csk_sz = CCSK_rotated_codes[0].size();
    uint16_t nm = dec_param.nm;
    float sigma = sqrt(1.0 / (pow(10, SNR / 10.0))); // N0/2 or N0?
    vector<uint16_t> NSYMB(N);
    vector<vector<uint16_t>> KBIN(K);

    // RandomBinaryGenerator(K, q, CCSK_rotated_codes, RepRndGn, 0, KBIN, KSYMB);
    // std::mt19937 gen(0);
    std::uniform_int_distribution<int> unif_dist(0, q - 1);
    for (uint16_t k = 0; k < K; k++)
    {
        KSYMB[k] = (uint16_t)unif_dist(gen);
    }
    for (int i = 0; i < K; i++)
        dec_param.ucap[dec_param.n][dec_param.reliab_sequence[i]] = KSYMB[i];
    Encoder(table.ADDGF, table.MULGF, dec_param.polar_coeff, dec_param.ucap, NSYMB);

    softdata_t modulation_I[q];
    softdata_t modulation_Q[q];
    for (int i = 0; i < q; i++)
    {
        modulation_I[i] = sin(i * 2 * PI / 64);
        modulation_Q[i] = cos(i * 2 * PI / 64);
    }
    vector<vector<vector<softdata_t>>> NBIN(N, vector<vector<softdata_t>>(csk_sz, vector<softdata_t>(2)));
    vector<vector<vector<softdata_t>>> rnd_noise(N, vector<vector<softdata_t>>(csk_sz, vector<softdata_t>(2)));
    vector<vector<softdata_t>> noisy_sym(csk_sz, vector<softdata_t>(2));
    vector<vector<softdata_t>> chan_LLR(N, vector<softdata_t>(q, 0));

    softdata_t min1;
    std::normal_distribution<double> norm_dist(0, sigma);

    for (int n = 0; n < N; n++)
    {
        min1 = std::numeric_limits<softdata_t>::max() / 2;
        // AWGN_gen(0, sigma, RepRndGn, 0, rnd_noise[n]);
        for (int i = 0; i < csk_sz; i++)
        {
            NBIN[n][i][0] = modulation_I[CCSK_rotated_codes[NSYMB[n]][i]];
            NBIN[n][i][1] = modulation_Q[CCSK_rotated_codes[NSYMB[n]][i]];
            // noisy_sym[i][0] = NBIN[n][i][0] + rnd_noise[n][i][0];
            // noisy_sym[i][1] = NBIN[n][i][1] + rnd_noise[n][i][1];
            noisy_sym[i][0] = NBIN[n][i][0] + (softdata_t)norm_dist(gen);
            noisy_sym[i][1] = NBIN[n][i][1] + (softdata_t)norm_dist(gen);
        }

        for (int k = 0; k < q; k++)
        {
            for (int i = 0; i < csk_sz; i++)
                chan_LLR[n][k] += pow(noisy_sym[i][0] - modulation_I[CCSK_rotated_codes[k][i]], 2) + pow(noisy_sym[i][1] - modulation_Q[CCSK_rotated_codes[k][i]], 2);
            chan_LLR[n][k] /= 2 * pow(sigma, 2);
        }
        for (int k = 0; k < q; k++)
        {
            chan_LLR[n][k] = chan_LLR[n][k] + 10 * log10(2 / pow(sigma, 2));
            if (chan_LLR[n][k] < min1)
                min1 = chan_LLR[n][k];
        }
        for (int k = 0; k < q; k++)
            chan_LLR[n][k] -= min1;
    }
    LLR_sort(chan_LLR, nm, chan_LLR_sorted, dec_param.offset);

    // save_intrinsic_LLR(chan_LLR, "LLR_to_test.txt");
    // appendSequenceToFile("SEQ_to_test.txt", KSYMB);
}
