// seqsim_windows_fast.cpp - same output as seqsim_windows.cpp (same SEQSIM Score), computed exactly
// but much faster: for every position i of X, SEQSIM's Weight needs the longest match starting at X[i]
// anywhere in Y. That number is the "matching statistic" of X against Y, which a suffix automaton of
// Y gives in one linear pass (run on the reversed strings so matches START at i). Each pair then costs
// O(N) instead of O(N^2). Verified identical (max difference 0) against the original seqsim.cpp.
//
// Build: cl /O2 /openmp /EHsc seqsim_windows_fast.cpp     or   g++ -O3 -fopenmp seqsim_windows_fast.cpp -o seqsim_windows_fast
// Run:   seqsim_windows_fast window.fa out_prefix [binsize=2000]
// Writes <out>_bins.tsv, <out>_ff.csv (direct), <out>_fr.csv (bin i vs reverse complement of bin j)
#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>

struct SAM {
    std::vector<int> len, link;
    std::vector<int> nx;  // 4 per state, -1 = none
    int last;
    static int code(char c) { return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3; }
    explicit SAM(const std::string& s) {
        len.reserve(2 * s.size() + 2); link.reserve(2 * s.size() + 2); nx.reserve(8 * s.size() + 8);
        newState(0, -1); last = 0;
        for (char ch : s) extend(code(ch));
    }
    int newState(int l, int lk) { len.push_back(l); link.push_back(lk); for (int k = 0; k < 4; k++) nx.push_back(-1); return (int)len.size() - 1; }
    void extend(int c) {
        int cur = newState(len[last] + 1, 0);
        int p = last;
        while (p != -1 && nx[4 * p + c] == -1) { nx[4 * p + c] = cur; p = link[p]; }
        if (p == -1) link[cur] = 0;
        else {
            int q = nx[4 * p + c];
            if (len[p] + 1 == len[q]) link[cur] = q;
            else {
                int cl = newState(len[p] + 1, link[q]);
                for (int k = 0; k < 4; k++) nx[4 * cl + k] = nx[4 * q + k];
                while (p != -1 && nx[4 * p + c] == q) { nx[4 * p + c] = cl; p = link[p]; }
                link[q] = link[cur] = cl;
            }
        }
        last = cur;
    }
};

// Weight(X, Y) exactly as in seqsim.cpp: sqrt(sum_i (longest_i - 1)^2 over longest_i > 0)
// revX = reverse(X); samRevY = suffix automaton of reverse(Y)
static double WeightFast(const std::string& revX, const SAM& samRevY)
{
    double result = 0;
    int v = 0, l = 0;
    for (char ch : revX) {
        int c = SAM::code(ch);
        while (v != 0 && samRevY.nx[4 * v + c] == -1) { v = samRevY.link[v]; l = samRevY.len[v]; }
        if (samRevY.nx[4 * v + c] != -1) { v = samRevY.nx[4 * v + c]; l++; }
        else { v = 0; l = 0; }
        if (l > 0) result += (double)((l - 1) * (l - 1));
    }
    return sqrt(result);
}

static std::string revcomp(const std::string& s)
{
    std::string r(s.rbegin(), s.rend());
    for (char& ch : r) ch = ch == 'A' ? 'T' : ch == 'T' ? 'A' : ch == 'C' ? 'G' : 'C';
    return r;
}

int main(int argc, char** argv)
{
    if (argc < 3) { printf("Usage: %s window.fa out_prefix [binsize=2000]\n", argv[0]); return 1; }
    int bin = argc > 3 ? atoi(argv[3]) : 2000;
    FILE* f = fopen(argv[1], "r");
    if (!f) { printf("cannot open %s\n", argv[1]); return 1; }
    std::string seq; char buf[4096]; long start1 = 1;
    while (fgets(buf, sizeof buf, f)) {
        if (buf[0] == '>') { const char* c = strchr(buf, ':'); if (c) start1 = atol(c + 1); continue; }
        for (char* p = buf; *p; ++p) if (*p != '\n' && *p != '\r') seq.push_back((char)toupper((unsigned char)*p));
    }
    fclose(f);

    std::vector<std::string> rF, rR;   // reversed forward bins, reversed revcomp bins
    std::vector<int> idx;
    int nb = (int)(seq.size() / bin);
    for (int b = 0; b < nb; b++) {
        std::string s = seq.substr((size_t)b * bin, bin);
        if (s.find_first_not_of("ACGT") != std::string::npos) continue;
        std::string rc = revcomp(s);
        rF.push_back(std::string(s.rbegin(), s.rend()));
        rR.push_back(std::string(rc.rbegin(), rc.rend()));
        idx.push_back(b);
    }
    int n = (int)rF.size();
    printf("%s: %d bins of %d bp, %d usable\n", argv[1], nb, bin, n);
    std::vector<SAM> sF, sR;
    sF.reserve(n); sR.reserve(n);
    for (int k = 0; k < n; k++) { sF.emplace_back(rF[k]); sR.emplace_back(rR[k]); }
    double N = (double)bin;
    double denom = sqrt(N * (N - 1) * (2 * N - 1) / 6);

    std::vector<double> ff((size_t)n * n, 0.0), fr((size_t)n * n, 0.0);
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j >= i) {
                double a = (i == j) ? 100.0 : (WeightFast(rF[i], sF[j]) + WeightFast(rF[j], sF[i])) / (2.0 * denom) * 100.0;
                ff[(size_t)i * n + j] = ff[(size_t)j * n + i] = a;
            }
            // Score(F_i, R_j) = (Weight(F_i, R_j) + Weight(R_j, F_i)) / (2 denom) * 100
            fr[(size_t)i * n + j] = (WeightFast(rF[i], sR[j]) + WeightFast(rR[j], sF[i])) / (2.0 * denom) * 100.0;
        }
    }

    std::string o(argv[2]);
    FILE* b = fopen((o + "_bins.tsv").c_str(), "w");
    fprintf(b, "kept_index\tbin_index\tstart_1based\n");
    for (int k = 0; k < n; k++) fprintf(b, "%d\t%d\t%ld\n", k, idx[k], start1 + (long)idx[k] * bin);
    fclose(b);
    FILE* m1 = fopen((o + "_ff.csv").c_str(), "w");
    FILE* m2 = fopen((o + "_fr.csv").c_str(), "w");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(m1, j ? ",%.5f" : "%.5f", ff[(size_t)i * n + j]);
            fprintf(m2, j ? ",%.5f" : "%.5f", fr[(size_t)i * n + j]);
        }
        fprintf(m1, "\n"); fprintf(m2, "\n");
    }
    fclose(m1); fclose(m2);
    printf("wrote %s_{bins.tsv,ff.csv,fr.csv}\n", argv[2]);
    return 0;
}
