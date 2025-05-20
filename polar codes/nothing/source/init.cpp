#include "GF_tools.h"
#include "struct.h"
#include "tools.h"
#include <fstream>
#include <iomanip>
#include <cstring>
#include <iostream>
#include "init.h"
#include <sstream>
#include <cmath>
void PoAwN::init::LoadCode(PoAwN::structures::base_code_t &code, float SNR)
{
    code.Rate = (float)code.K / (float)code.N;
    code.reliab_sequence.resize(code.q, 0);

    std::ostringstream fname;
    std::string mat_direct;
    if (code.sig_mod == "BPSK")
        mat_direct = "./matrices/bpsk/N";
    else if (code.sig_mod == "CCSK_BIN")
        mat_direct = "./matrices/ccsk_bin/N";
    else
        mat_direct = "./matrices/ccsk_nb/N";

    fname << mat_direct << code.N << "/mat_N" << code.N << "_GF"
          << code.q << "_SNR" << std::fixed << std::setprecision(3) << SNR
          << ".txt";

    std::cout << fname.str() << std::endl;
    std::ifstream opfile(fname.str());
    if (!opfile)
    {
        std::cerr << "Sequence Unavailable!!" << std::endl;
        exit(-1010);
    }

    int tmp;
    int temp_cnt = 0;
    code.reliab_sequence.resize(code.N);
    for (int k = 0; k < 2; k++) // if k<1 then read the first reliability sequence existed in the file, if k<2 then cosider the second.
    // I did this loop because my files contain 2  sequnces, the first generated taking Entropies as measurments, while the second is
    //  generated from probability of errors. what has been noticed is at very low SNR (FER<0.1) taking probability of error is better
    {
        temp_cnt = 0;
        for (int i = 0; i < code.N; i++)
        {
            if (opfile >> tmp) 
            {
                code.reliab_sequence[i] = static_cast<uint16_t>(tmp);
                temp_cnt++;
            }
            else
            {
                std::cerr << "Invalid or missing value at index " << i << " of Reliabilities sequence nb  " << k << std::endl;
                exit(EXIT_FAILURE);
            }

        }
        if (temp_cnt < code.N)
        {
            std::cerr << "Reliabilities sequence nb " << k << " is not complete" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    if (code.sig_mod == "BPSK")
    {
        code.polar_coeff.resize(code.n, std::vector<uint16_t>(code.N / 2));
        // Read polar coefficients
        for (int i = 0; i < code.n; i++)
        {
            for (int j = 0; j < code.N / 2; j++)
            {
                opfile >> tmp;
                code.polar_coeff[i][j] = tmp;
            }
        }
    }
    else
    {
        code.polar_coeff.resize(code.n, std::vector<uint16_t>(code.N / 2));
        for (int i = 0; i < code.n; i++)
        {
            for (int j = 0; j < code.N / 2; j++)
            {
                code.polar_coeff[i][j] = 1;
            }
        }
    }
}

void PoAwN::init::LoadBubblesIndcatorlists(PoAwN::structures::decoder_parameters &dec,
                                           const float SNR)
{
    uint16_t n = dec.n, nH = dec.nH, nL = dec.nL;
    std::ostringstream fname;
    std::string mat_direct;
    if (dec.sig_mod == "BPSK")
        mat_direct = "./BubblesPattern/bpsk/N";
    else if (dec.sig_mod == "CCSK_BIN")
        mat_direct = "./BubblesPattern/ccsk_bin/N";
    else
        mat_direct = "./BubblesPattern/ccsk_nb/N";

    fname << mat_direct << dec.N
          << "/BubblesIndicatorsLists"
          << "/bubbles_N" << dec.N
          << "_K" << dec.K
          << "_GF" << dec.q
          << "_SNR" << std::fixed << std::setprecision(3) << SNR << "_"
          << dec.nH << "x" << dec.nL
          << "_Bt_lsts.txt";
    std::string filename = fname.str();

    std::ifstream file(filename);
    std::vector<std::string> lines;
    std::string line;
    int linecount = 0;
    while (std::getline(file, line))
    {
        if (!line.empty())
        {
            lines.push_back(line);
            linecount++;
        }
    }
    if (linecount < (1u << n) - 1)
    {
        std::cerr << "File data is not enough, it should  contain " << (1u << n) - 1
                  << " lines, each contains the list of (i,j) coordinates of cluster s at layer l, the lines format should be: l s, i0 j0 i1 j1..." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    int line_cnt0, line_cnt = 0;
    dec.Bubb_Indicator.resize(n);
    for (int l = 0; l < n; l++)
    {
        dec.Bubb_Indicator[l].resize(1 << l);
        line_cnt0 = (1u << l) - 1;
        for (int s = 0; s < (1u << l); s++)
        {
            line_cnt = line_cnt0 + s;
            dec.Bubb_Indicator[l][s].resize(2);
            {
                std::istringstream iss(lines[line_cnt]);
                std::vector<int> numbers;
                std::string token;
                int skipped = 0;

                while (iss >> token)
                {

                    if (skipped < 2)
                    {
                        skipped++;
                        continue;
                    }

                    std::string clean;
                    for (char c : token)
                        if (std::isdigit(c))
                            clean += c;
                        else if (!clean.empty())
                        {
                            numbers.push_back(std::stoi(clean));
                            clean.clear();
                        }
                    if (!clean.empty())
                        numbers.push_back(std::stoi(clean));
                }

                for (size_t i = 0; i < numbers.size(); i += 2)
                {
                    dec.Bubb_Indicator[l][s][0].push_back(numbers[i]);
                    dec.Bubb_Indicator[l][s][1].push_back(numbers[i + 1]);
                }
            }
        }
    }
}

// void PoAwN::init::LoadBubblesIndcatorlists(PoAwN::structures::decoder_parameters &dec,
//                                            const float SNR,
//                                            const float Pt)
// {
//     char fname[70];
//     sprintf(fname, "./sizes/N%d/bubblesizemap_N%d_K%d_GF%d.txt", dec.N, dec.N, dec.K, dec.q);
//     FILE *opfile;
//     opfile = fopen(fname, "r");
//     int N = dec.N, tmp;

//     for (int i = 0; i < N - 1; i++)
//     {
//         fscanf(opfile, "%d", &tmp);
//         if (i == 0)
//             dec.nm = tmp;
//     }
//     dec.Bubb_Indicator.resize(dec.n);
//     dec.Bubb_Indicator[0].resize(1);
//     int p = 0, p1 = -1;
//     for (int d = 0; d < dec.N - 1; d++)
//     {
//         p1++;
//         dec.Bubb_Indicator[p][p1].resize(2);

//         for (int n = 0; n < dec.nm; n++)
//         {
//             for (int q = 0; q < dec.nm; q++)
//             {
//                 fscanf(opfile, "%d", &tmp);
//                 if (tmp == 1)
//                 {
//                     dec.Bubb_Indicator[p][p1][0].push_back(n);
//                     dec.Bubb_Indicator[p][p1][1].push_back(q);
//                 }
//             }
//         }
//         if (d == pow(2,p+1)-2)
//         {
//             p1 = -1;
//             p++;
//             if(p==dec.n)
//                 break;
//             dec.Bubb_Indicator[p].resize(1 << (p));
//         }
//     }
//     fclose(opfile);
// }
void PoAwN::init::LoadTables(PoAwN::structures::base_code_t &code,
                             PoAwN::structures::table_GF &table,
                             const uint16_t *GF_polynom_primitive)
{
    uint16_t prim_pol = GF_polynom_primitive[code.p - 2];
    PoAwN::GFtools::GF_bin_seq_gen(code.q, prim_pol, table.BINGF, table.BINDEC);
    PoAwN::GFtools::GF_bin2GF(table.BINGF, table.DECGF, table.GFDEC); // x^-inf=GFDEC[0]=0, x^0=GFDEC[1]=1, , x^1=GFDEC[2]=2,.., , x^(q-2)=GFDEC[q-1] and GFDEC[DECGF[i]]=i
    PoAwN::GFtools::GF_add_mat_gen(table.BINGF, table.DECGF, table.GFDEC, table.ADDGF, table.ADDDEC);
    PoAwN::GFtools::GF_mul_mat_gen(table.DECGF, table.MULGF, table.MULDEC);
    PoAwN::GFtools::GF_div_mat_gen(table.DECGF, table.DIVGF, table.DIVDEC);
}