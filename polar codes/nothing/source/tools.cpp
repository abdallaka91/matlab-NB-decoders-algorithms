#include "tools.h"
#include "struct.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
#include <random>
// #include "debugging_tools.h"

using std::vector;

void PoAwN::tools::circshift(vector<uint16_t> &CCSK_seq, int shift)
{
    const uint16_t q = CCSK_seq.size();
    const auto shift_masked = shift & (q - 1);
    rotate(CCSK_seq.begin(), CCSK_seq.begin() + (q - shift_masked), CCSK_seq.end());
}

void PoAwN::tools::create_ccsk_rotated_table(const vector<uint16_t> &CCSK_seq,
                                             const uint16_t q,
                                             vector<vector<uint16_t>> &ccsk_rotated_table)
{
    std::vector<uint16_t> eta0(CCSK_seq.begin(), CCSK_seq.begin() + q);
    ccsk_rotated_table[0] = eta0;

    for (int i = 1; i < q; ++i)
    {
        circshift(eta0, -1);
        ccsk_rotated_table[i] = eta0;
    }
}
void PoAwN::tools::RandomBinaryGenerator(const uint16_t K,
                                         const uint16_t q,
                                         const vector<vector<uint16_t>> &bin_table,
                                         const bool repeatable,
                                         const int SEED,
                                         vector<vector<uint16_t>> &KBIN,
                                         vector<uint16_t> &KSYMB)
{
    for (uint16_t k = 0; k < K; k++)
    {
        static std::mt19937 gen(repeatable ? SEED : std::random_device{}());
        static std::uniform_int_distribution<int> dist(0, q - 1);
        uint16_t randv = (uint16_t)dist(gen);
        KSYMB[k] = randv;
        KBIN[k] = bin_table[randv];
    }
}

void PoAwN::tools::RandomSymbGenerator(const uint16_t K,
                                       const uint16_t q,
                                       const bool repeatable,
                                       const int SEED,
                                       vector<uint16_t> &KSYMB)
{
    for (uint16_t k = 0; k < K; k++)
    {
        static std::mt19937 gen(repeatable ? SEED : std::random_device{}());
        static std::uniform_int_distribution<int> dist(0, q - 1);
        uint16_t randv = (uint16_t)dist(gen);
        KSYMB[k] = randv;
    }
}

void PoAwN::tools::AWGN_gen(double MEAN,
                            double STD,
                            bool repeatable,
                            double SEED,
                            vector<vector<softdata_t>> &noise_table)
{
    static std::mt19937 gen(repeatable ? std::mt19937(SEED) : std::mt19937(std::random_device{}()));
    int N = noise_table.size();
    for (int i = 0; i < N; i++)
    {
        int q1 = noise_table[i].size();
        for (int j = 0; j < q1; j++)
        {
            std::normal_distribution<double> dist(MEAN, STD);
            noise_table[i][j] = (softdata_t)dist(gen);
        }
    }
}

void PoAwN::tools::awgn_channel_noise(const double sigma,
                                      const bool repeatable,
                                      const double SEED,
                                      vector<vector<softdata_t>> &noisy_sig)
{
    uint16_t N = noisy_sig.size();
    uint16_t q1 = noisy_sig[0].size();
    vector<vector<softdata_t>> noise_table;
    noise_table.resize(N, vector<softdata_t>(q1, 0));
    AWGN_gen(SEED, sigma, repeatable, SEED, noise_table);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < q1; j++)
            noisy_sig[i][j] += noise_table[i][j];
}

void PoAwN::tools::Encoder(const vector<vector<uint16_t>> &ADDGF, const vector<vector<uint16_t>> &MULGF,
                           const vector<vector<uint16_t>> &polar_coeff,
                           vector<vector<uint16_t>> &ucap, vector<uint16_t> &NSYMB)
{
    uint16_t N = ucap[0].size();
    uint16_t n = ucap.size() - 1;

    uint16_t a, b;
    uint16_t tmp_add;
    uint16_t tmp_mul;
    for (int l = n - 1; l >= 0; l--)
    {
        for (uint16_t k = 0; k < N / 2; k++)
        {
            uint16_t pw1 = pow(2, n - l - 1);
            a = 2 * k - k % pw1;
            b = pw1 + 2 * k - k % pw1;
            tmp_add = ADDGF[ucap[l + 1][a]][ucap[l + 1][b]];
            ucap[l][a] = tmp_add;
            tmp_mul = MULGF[ucap[l + 1][b]][polar_coeff[n - l - 1][k]];
            ucap[l][b] = tmp_mul;
        }
    }
    for (uint16_t i = 0; i < N; ++i)
        NSYMB[i] = ucap[0][i];
}
float PoAwN::tools::My_drand48(int *initialise)
{
    static thread_local std::mt19937 generator(std::random_device{}());
    static thread_local std::uniform_real_distribution<float> distribution(0.0f, 1.0f);

    if (*initialise == -1)
    {
        generator.seed(std::random_device{}());
        *initialise = 0;
    }

    return distribution(generator);
}

// void PoAwN::tools::inv_Encoder(const vector<vector<uint16_t>> &ADDGF, const vector<vector<uint16_t>> &DIVGF,
//                                const vector<vector<uint16_t>> &polar_coeff,
//                                const vector<uint16_t> NSYMB, vector<uint16_t> &u_symb)
// {
//     uint16_t N = NSYMB.size();
//     uint16_t n = log2(N);
//     vector<vector<uint16_t>> temp_symb(N, vector<uint16_t>(n + 1, 0));
//     for (uint16_t i = 0; i < N; ++i)
//         temp_symb[i][n] = NSYMB[i];

//     uint16_t a, b;
//     uint16_t tmp_add;
//     uint16_t tmp_div;
//     for (int l = n - 1; l >= 0; l--)
//     {
//         for (uint16_t k = 0; k < N / 2; k++)
//         {
//             uint16_t pw1 = pow(2, l );
//             a = 2 * k - k % pw1;
//             b = pw1 + 2 * k - k % pw1;
//             tmp_div = DIVGF[temp_symb[b][l + 1]][polar_coeff[l][k]];
//             tmp_add = ADDGF[temp_symb[a][l + 1]][tmp_div];
//             temp_symb[a][l] = tmp_add;
//             temp_symb[b][l] = tmp_div;
//         }
//     }
//     for (uint16_t i = 0; i < N; ++i)
//         u_symb[i] = temp_symb[i][0];
// }
