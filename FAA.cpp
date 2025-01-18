
/*-------------------------------------------------------*/
/*
    Library: New Project, win32, static library, no headers. Include all src files, choose MFC, compile.
    Include compiled library to project, include header files.
*/
/*-------------------------------------------------------*/

/* packages */
#include <NTL/mat_GF2.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <cstdint>
#include <numeric>
#include <iostream>

/* namespaces */
using namespace ::NTL;
using namespace ::std;

/* new types */
typedef uint64_t u64;       // 64-bit
typedef unsigned int u32;   // 32-bit
typedef unsigned short u16; // 16-bit
typedef unsigned char u8;   //  8-bit

/* global variables */
/*----------------------- SETUP -------------------------*/

u32 n, d, e;
u32 N, D, E;

/*-------------------------------------------------------*/

/* declaration of functions */

u32 Algorithm1(u8 *TT, u8 **SolutionLst);
void hadamard(u8 *A, u8 *B);
u32 degree(u8 *A);
u32 hamming(u32 x);
void dispEq(u8 *CT, string &ANF);
void ArrayToString(u8 *A, string &str);

u32 rnd();
u32 tf0(u64 &x);

float factln(int nn);
float bico(int nn, int k);
float gammln(float xx);
void binary(u32 x, u32 s);

/*-------------------------------------------------------*/
/* main */

int main()
{
    u32 i, j;
    u32 dim, degf, degg, degh;
    u32 seed;
    time_t start, end;
    double diffms;
    char c;
    u8 error = 0;

    /*-------------------------------------------------------*/

    seed = time(NULL);
    srand(seed);

    /*-------------------------------------------------------*/

    printf("\nFAA Equation Finder\n");
    printf("-------------------\n\n");
    printf("Written by Simon Fischer, www.simonfischer.ch, 2008\n\n");

    printf("Input:  Boolean function f with n variables, degree e and degree d\n");
    printf("Output: Boolean functions g and h with f*g=h and deg g<=e and deg h<=d\n\n");

    printf("Read input.txt and write output.txt\n");

    printf("\nChoose n: ");
    cin >> n;
    printf("Choose e: ");
    cin >> e;
    printf("Choose d: ");
    cin >> d;

    string InputFileName = "input.txt";
    string OutputFileName = "outputi_18_5_9_linux.txt";
    // printf("Choose input file: ");  cin >> InputFileName;
    // printf("Choose output file: "); cin >> OutputFileName;

    /*-------------------------------------------------------*/

    // Compute N,E,D
    N = (u32)1 << n;
    for (E = 0, i = 0; i <= e; i++)
    {
        E += (u32)bico(n, i);
    }
    for (D = 0, i = 0; i <= d; i++)
    {
        D += (u32)bico(n, i);
    }

    // Prepare Array for Solutions CTf
    u8 **SolutionLst = new u8 *[E];
    for (i = 0; i < E; i++)
    {
        SolutionLst[i] = new u8[N];
    }

    u8 *TTf = new u8[N];
    u8 *CTf = new u8[N];
    u8 *TTg = new u8[N];
    u8 *CTg = new u8[N];
    u8 *TTh = new u8[N];
    u8 *CTh = new u8[N];
    string ANF = "";
    string table = "";
    string text = "";

    /*-------------------------------------------------------*/

    FILE *fInput = fopen(InputFileName.c_str(), "r");
    FILE *fOutput = fopen(OutputFileName.c_str(), "w");

    if (fInput == NULL)
    {
        printf("\nError: Could not open input file\n");
        error = 1;
    }
    else
    {
        i = 0;
        while ((c = getc(fInput)) != EOF)
        {
            if (c == '0')
            {
                if (i < N)
                {
                    TTf[i] = 0;
                }
                i++;
            }
            if (c == '1')
            {
                if (i < N)
                {
                    TTf[i] = 1;
                }
                i++;
            }
        }
        if (i != N)
        {
            error = 1;
            printf("\nError: Number of 0 and 1 in input file does not correspond to 2^n\n");
        }
    }

    if (error == 0)
    {

        /*-------------------------------------------------------*/

        // Compute the solutions CTg
        printf("\nSet up and solve matrix of size %u x %u", (N - D), E);

        start = clock();
        dim = Algorithm1(TTf, SolutionLst);
        end = clock();

        diffms = ((end - start) * 1000) / CLOCKS_PER_SEC;
        printf(" (%.0f ms)\n\n", diffms);

        printf("\nNumber of independent solutions = %u\n\n\n", dim);

        // Transform TTf to CTf and compute degree of f
        hadamard(TTf, CTf);
        degf = degree(CTf);

        // Display TT,CT,ANF of f
        table = "TT(f) = [";
        ArrayToString(TTf, table);
        table += "]";
        // printf("Input f (deg f=%u)\n",degf);
        // printf("-------\n\n");
        // printf(table.c_str()); printf("\n");
        fprintf(fOutput, "Input f (deg f=%u)\n", degf);
        fprintf(fOutput, "-------\n\n");
        fprintf(fOutput, table.c_str());
        fprintf(fOutput, "\n");

        table = "CT(f) = [";
        ArrayToString(CTf, table);
        table += "]";
        // printf(table.c_str()); printf("\n");
        fprintf(fOutput, table.c_str());
        fprintf(fOutput, "\n");

        ANF = "f = ";
        dispEq(CTf, ANF);
        // printf(ANF.c_str()); printf("\n\n");
        fprintf(fOutput, ANF.c_str());
        fprintf(fOutput, "\n\n");

        for (i = 0; i < dim; i++)
        {
            // Choose a solution of CTg
            for (j = 0; j < N; j++)
            {
                CTg[j] = SolutionLst[i][j];
            }

            // Compute all remaining CT and TT of g and h
            hadamard(CTg, TTg);
            degg = degree(CTg);
            for (j = 0; j < N; j++)
            {
                TTh[j] = TTf[j] * TTg[j];
            }
            hadamard(TTh, CTh);
            degh = degree(CTh);

            // printf("Solution %u (deg g=%u, deg h=%u)\n",i+1,degg,degh);
            // printf("--------\n\n");
            fprintf(fOutput, "Solution %u (deg g=%u, deg h=%u)\n", i + 1, degg, degh);
            fprintf(fOutput, "--------\n\n");

            // Display TT,CT,ANF of g
            table = "TT(g) = [";
            ArrayToString(TTg, table);
            table += "]";
            // printf(table.c_str()); printf("\n");
            fprintf(fOutput, table.c_str());
            fprintf(fOutput, "\n");

            table = "CT(g) = [";
            ArrayToString(CTg, table);
            table += "]";
            // printf(table.c_str()); printf("\n");
            fprintf(fOutput, table.c_str());
            fprintf(fOutput, "\n");

            ANF = "g = ";
            dispEq(CTg, ANF);
            // printf(ANF.c_str()); printf("\n\n");
            fprintf(fOutput, ANF.c_str());
            fprintf(fOutput, "\n\n");

            // Display TT,CT,ANF of h
            table = "TT(h) = [";
            ArrayToString(TTh, table);
            table += "]";
            // printf(table.c_str()); printf("\n");
            fprintf(fOutput, table.c_str());
            fprintf(fOutput, "\n");

            table = "CT(h) = [";
            ArrayToString(CTh, table);
            table += "]";
            // printf(table.c_str()); printf("\n");
            fprintf(fOutput, table.c_str());
            fprintf(fOutput, "\n");

            ANF = "h = ";
            dispEq(CTh, ANF);
            // printf(ANF.c_str()); printf("\n\n");
            fprintf(fOutput, ANF.c_str());
            fprintf(fOutput, "\n\n");
        }
    } // if InputFile

    delete (SolutionLst);
    delete (TTf);
    delete (CTf);
    delete (TTg);
    delete (CTg);
    delete (TTh);
    delete (CTh);

    // fclose(fInput);
    fclose(fOutput);

    printf("\n\n");
    system("pause");
    return 0;
}

/*-------------------------------------------------------*/

void ArrayToString(u8 *A, string &str)
{
    u32 i;
    char buffer[33];
    for (i = 0; i < N; i++)
    {
        str += std::to_string(A[i]);
    }
}

/*-------------------------------------------------------*/

// Compute Walsh-Hadamard transfrom from A to B (e.g. TT to CT, or CT to TT)
void hadamard(u8 *A, u8 *B)
{
    u32 i, j;
    for (i = 0; i < N; i++)
    {
        B[i] = 0;
        for (j = 0; j < N; j++)
        {
            B[i] ^= ((i | j) == i) & A[j];
        }
    }
}

/*-------------------------------------------------------*/

// Compute degree of CT A
u32 degree(u8 *A)
{
    u32 i, deg, tmp;
    deg = 0;
    for (i = 0; i < N; i++)
    {
        if (A[i] == 1)
        {
            tmp = hamming(i);
            if (tmp > deg)
            {
                deg = tmp;
            }
        }
    }
    return (deg);
}

/*-------------------------------------------------------*/

u32 hamming(u32 x)
{
    u32 weight = 0, i;
    for (i = 0; i < 32; i++, x >>= 1)
    {
        weight += (u32)(x & 1);
    }
    return weight;
}

/*-------------------------------------------------------*/

// Gibt einen String aus, der die ersten "n" Bit der Zahl "x" enthält
void binary(u32 x, u32 s)
{
    u32 y = x;
    for (u32 i = 0; i < s; i++)
    {
        y = x >> (s - i - 1);
        if (y & 1)
            printf("1");
        else
            printf("0");
    }
}

/*-------------------------------------------------------*/

u32 tf0(u64 &x)
{
    // Update
    x = (x + ((x * x) | 5)) & 0xffffffffffffffff;

    // Outupt
    return x >> 32;
}

/*-------------------------------------------------------*/

u32 Algorithm1(u8 *TT, u8 **SolutionLst)
{
    u32 i, j, k;
    u32 input;
    u32 Gamma, Beta;
    u32 rows, cols;
    rows = N - D;
    cols = E;

    u32 *BetaLst = new u32[N];  // N-D
    u32 *GammaLst = new u32[N]; // E
    for (i = 0; i < N; i++)
    {
        BetaLst[i] = 0;
        GammaLst[i] = 0;
    }

    for (j = 0, i = 0; i < N; i++)
    {
        if (hamming(i) > d)
        {
            GammaLst[j] = i;
            j++;
        }
    }

    for (j = 0, i = 0; i < N; i++)
    {
        if (hamming(i) <= e)
        {
            BetaLst[j] = i;
            j++;
        }
    }

    mat_GF2 M, Mt, Mk;
    M.SetDims((long)rows, (long)cols);
    Mt.SetDims((long)cols, (long)rows);
    clear(M);
    clear(Mt);
    clear(Mk);

    for (i = 0; i < rows; i++)
    {
        Gamma = GammaLst[i];
        for (j = 0; j < cols; j++)
        {
            Beta = BetaLst[j]; // eventuell nur Beta<=Gamma (bitweise)
            for (input = 0, k = Beta; k <= Gamma; k++)
            {
                input ^= ((((Beta | k) == k) & ((k | Gamma) == Gamma)) & TT[k]) & 0x1;
            }
            input &= 0x1;
            M[i][j] = input;
            // printf("%u ",input); if(j==6) printf("\n");
        }
    }

    /* compute kernel with NTL */
    u32 dim;
    transpose(Mt, M);
    kernel(Mk, Mt);
    dim = Mk.NumRows();

    u32 index;

    for (i = 0; i < E; i++)
    {
        for (j = 0; j < N; j++)
            SolutionLst[i][j] = 0;
    }

    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < cols; j++)
        {
            index = BetaLst[j];
            if (Mk[i][j] == 1)
            {
                SolutionLst[i][index] = 1;
            }
            else
            {
                SolutionLst[i][index] = 0;
            }
        }

    } // for dim

    delete (BetaLst);
    delete (GammaLst);

    return dim;
}

/*-------------------------------------------------------*/

void dispEq(u8 *CT, string &ANF)
{

    u32 i, j;
    u32 monomial;
    char buffer[33];
    u32 flag = 0;
    u32 zero = 0;
    for (i = 0; i < N; i++)
    {
        if (CT[i] & 1)
        {
            zero = 1;
            monomial = i;
            if (flag == 1)
                ANF += " + ";
            if (monomial == 0)
                ANF += "1";
            for (j = 1; j <= n; j++)
            {
                if (monomial & 1)
                {
                    ANF += "x";
                    ANF += to_string(j);
                    ANF += " ";
                } //%d ",j);}
                monomial >>= 1;
            }
            flag = 1;
        }
    }
    if (zero == 0)
        ANF += "0";
    // printf("\n");
}

/*-------------------------------------------------------*/

/* from Numerical Recipes in C: binomial coef. computation */
float factln(int nn)
{
    /* Returns ln(nn!). */
    float gammln(float xx);

    static float a[101]; // A static array is automatically initialized to zero.
    if (nn <= 1)
        return 0.0;
    if (nn <= 100)
        return a[nn] ? a[nn] : (a[nn] = gammln(nn + 1.0)); // In range of table.
    else
        return gammln(nn + 1.0); // Out of range of table.
}
float bico(int nn, int k)
/* Returns the binomial coefficient (nn k) as a floating-point number. */
{
    float factln(int nn);
    return floor(0.5 + exp(factln(nn) - factln(k) - factln(nn - k)));
    /* The floor function cleans up roundoff error for smaller values of n and k. */
}
float gammln(float xx)
/* Returns the value ln[Γ(xx)] for xx > 0. */
{

    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677,
                            24.01409824083091, -1.231739572450155,
                            0.1208650973866179e-2, -0.5395239384953e-5};
    int j;
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++)
        ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

/*-------------------------------------------------------*/
