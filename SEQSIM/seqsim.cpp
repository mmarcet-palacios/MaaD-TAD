#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

#define MAX_LEN 30000

// Simplified weight calculation
// Performance could be optimized with suffix tree
double Weight(std::string const& X, std::string const& Y)
{
    int N = X.size();

    double result = 0;
    for (int i = 0; i < N; i++)
    {
        int longest = 0;
        for (int j = 0; j < N; j++)
        {
            int k = 0;
            while (i + k < N && j + k < N && X[i + k] == Y[j + k]) k++;
            if (k > longest) longest = k;
        }
        if (longest > 0) result += (longest - 1) * (longest - 1);
    }
    return sqrt(result);
}

double Score(std::string const& X, std::string const& Y)
{
    double N = (double)X.size();
    double xy = Weight(X, Y);
    double yx = Weight(Y, X);
    double denom = sqrt(N * (N - 1) * (2 * N - 1) / 6);
    return (xy + yx) / (2.0 * denom) * 100.0;
}

int main(int argc, char** argv)
{
    // Display usage instructions
    if (argc != 4)
    {
        printf("Usage: %s inputFile firstIdx lastIdx\n", argv[0]);
        printf("  inputFile: Text file with whitespace-separated pairs of gene names and gene strings\n");
        printf("             Gene strings consisting of letters outside of 'A','C','G','T' will be ignored\n");
        printf("  firstIdx : 0-based index of first gene to compute\n");
        printf("  lastIdx  : 0-based index of last gene to compute\n");
        exit(1);
    }

    // Validate command-line arguments
    std::string pramef = argv[1];
    int firstIdx = atoi(argv[2]);
    int lastIdx = atoi(argv[3]);
    if (firstIdx < 0)
    {
        printf("[ERROR] firstIdx must be non-negative\n");
        exit(1);
    }
    if (lastIdx < firstIdx)
    {
        printf("[ERROR] lastIdx must be greater or equal to firstIdx\n");
        exit(1);
    }

    // Read in input file
    FILE* input = fopen(pramef.c_str(), "r");
    if (!input)
    {
        printf("Unable to open inputFile %s\n", pramef.c_str());
        exit(1);
    }

    // Read in genes
    std::vector<std::pair<std::string, std::string>> genes;
    char allowed[] = "ACGT";
    char name[MAX_LEN + 1];
    char gene[MAX_LEN + 1];
    size_t len = 0;
    while (fscanf(input, " %3000s %3000s", name, gene) == 2)
    {
        size_t geneLen = strlen(gene);
        if (len == 0)
            len = geneLen;
        else if (geneLen != len)
        {
            printf("[ERROR] Gene strings must all be of the same length (Expected length %lu, not %lu)\n", len, geneLen);
            exit(1);
        }
        if (strspn(gene, allowed) != geneLen)
        {
            printf("[WARN] Skipping %s due to unrecognized letters outside of [%s]\n", name, allowed);
            continue;
        }
        genes.push_back(std::make_pair(name, gene));
    }
    int numGenes = genes.size();
    printf("Found %d genes of length %lu\n", numGenes, len);
    if (lastIdx >= numGenes)
    {
        printf("[ERROR] lastIdx must be less than %d\n", numGenes);
        exit(1);
    }

    // Open output file
    char outputFilename[100];
    sprintf(outputFilename, "%05d-%05d.csv", firstIdx, lastIdx);
    printf("Writing results to %s\n", outputFilename);
    FILE* output = fopen(outputFilename, "w");
    if (!output)
    {
        printf("[ERROR] Unable to open output file %s\n", outputFilename);
        exit(1);
    }

    // Print off header row
    fprintf(output, "Name");
    for (int colIdx = firstIdx; colIdx <= lastIdx; ++colIdx)
    {
        fprintf(output, ",%s", genes[colIdx].first.c_str());
    }
    fprintf(output, "\n");

    // Comptue results in parallel
    std::vector<std::vector<double>> results(numGenes);
    for (int rowIdx = 0; rowIdx < numGenes; rowIdx++)
        results[rowIdx].resize(lastIdx - firstIdx + 1);

#pragma omp parallel for collapse(2)
    for (int rowIdx = 0; rowIdx < numGenes; ++rowIdx)
        for (int colIdx = firstIdx; colIdx <= lastIdx; ++colIdx)
            results[rowIdx][colIdx - firstIdx] = Score(genes[rowIdx].second, genes[colIdx].second);

    for (int rowIdx = 0; rowIdx < numGenes; ++rowIdx)
    {
        fprintf(output, "%s", genes[rowIdx].first.c_str());
        for (int colIdx = firstIdx; colIdx <= lastIdx; ++colIdx)
            fprintf(output, ",%.5f", results[rowIdx][colIdx - firstIdx]);
        fprintf(output, "\n");
    }
    fclose(output);

    return 0;
}
