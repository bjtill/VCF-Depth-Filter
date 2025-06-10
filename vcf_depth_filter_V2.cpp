//BT June 4, 2025
//g++ -o vcf_dpfilter vcf_depth_filter.cpp -lz
//works with either .gz or not

//# Read compressed VCF, write compressed output
//./vcf_dpfilter --input variants.vcf.gz --output filtered.vcf.gz --min-depth 10 --max-depth 100

//# Read compressed VCF, write uncompressed output
//./vcf_dpfilter --input variants.vcf.gz --output filtered.vcf --min-depth 10 --max-depth 100

//# Read uncompressed VCF, write compressed output
//./vcf_dpfilter --input variants.vcf --output filtered.vcf.gz --min-depth 10 --max-depth 100


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <getopt.h>
#include <cstdlib>
#include <climits>
#include <zlib.h>

class VCFDepthFilter {
private:
    int minDepth;
    int maxDepth;
    std::string inputFile;
    std::string outputFile;

public:
    VCFDepthFilter() : minDepth(0), maxDepth(INT_MAX) {}

    void parseArguments(int argc, char* argv[]) {
        int opt;
        static struct option long_options[] = {
            {"min-depth", required_argument, 0, 'd'},
            {"max-depth", required_argument, 0, 'D'},
            {"input", required_argument, 0, 'i'},
            {"output", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        while ((opt = getopt_long(argc, argv, "d:D:i:o:h", long_options, NULL)) != -1) {
            switch (opt) {
                case 'd':
                    minDepth = std::atoi(optarg);
                    break;
                case 'D':
                    maxDepth = std::atoi(optarg);
                    break;
                case 'i':
                    inputFile = optarg;
                    break;
                case 'o':
                    outputFile = optarg;
                    break;
                case 'h':
                    printUsage();
                    exit(0);
                default:
                    printUsage();
                    exit(1);
            }
        }

        if (inputFile.empty()) {
            std::cerr << "Error: Input file is required\n";
            printUsage();
            exit(1);
        }

        if (outputFile.empty()) {
            outputFile = "filtered_" + inputFile;
        }
    }

    void printUsage() {
        std::cout << "Usage: vcf_filter [OPTIONS]\n"
                  << "Options:\n"
                  << "  --min-depth, -d <int>    Minimum depth threshold (default: 0)\n"
                  << "  --max-depth, -D <int>    Maximum depth threshold (default: unlimited)\n"
                  << "  --input, -i <file>       Input VCF file (required)\n"
                  << "  --output, -o <file>      Output VCF file (default: filtered_<input>)\n"
                  << "  --help, -h               Show this help message\n";
    }

    std::vector<std::string> split(const std::string& str, char delimiter) {
        std::vector<std::string> tokens;
        std::stringstream ss(str);
        std::string token;
        while (std::getline(ss, token, delimiter)) {
            tokens.push_back(token);
        }
        return tokens;
    }

    bool isMissingGenotype(const std::string& genotype) {
        return genotype == "./." || genotype == ".|." || genotype.empty();
    }

    int extractDepth(const std::string& format, const std::string& sampleData) {
        std::vector<std::string> formatFields = split(format, ':');
        std::vector<std::string> sampleFields = split(sampleData, ':');

        // Find DP field index
        int dpIndex = -1;
        for (int i = 0; i < formatFields.size(); i++) {
            if (formatFields[i] == "DP") {
                dpIndex = i;
                break;
            }
        }

        if (dpIndex == -1 || dpIndex >= sampleFields.size()) {
            return -1; // DP field not found or missing
        }

        try {
            return std::stoi(sampleFields[dpIndex]);
        } catch (const std::exception& e) {
            return -1; // Invalid depth value
        }
    }

    bool isGzipped(const std::string& filename) {
        return filename.size() >= 3 && filename.substr(filename.size() - 3) == ".gz";
    }

    bool readLine(gzFile& gzInput, std::ifstream& input, std::string& line, bool useGzip) {
        if (useGzip) {
            char buffer[65536];
            if (gzgets(gzInput, buffer, sizeof(buffer)) != NULL) {
                line = buffer;
                // Remove trailing newline if present
                if (!line.empty() && line.back() == '\n') {
                    line.pop_back();
                }
                return true;
            }
            return false;
        } else {
            return static_cast<bool>(std::getline(input, line));
        }
    }

    bool passesDepthFilter(const std::string& line) {
        std::vector<std::string> fields = split(line, '\t');
        
        if (fields.size() < 10) {
            return false; // Not enough fields for samples
        }

        std::string format = fields[8]; // FORMAT field
        
        // Check each sample (starting from field 9)
        for (int i = 9; i < fields.size(); i++) {
            std::string sampleData = fields[i];
            
            // Extract genotype (first field before first colon)
            std::string genotype = split(sampleData, ':')[0];
            
            // Skip missing genotypes
            if (isMissingGenotype(genotype)) {
                continue;
            }

            // Extract depth
            int depth = extractDepth(format, sampleData);
            
            // If we can't extract depth, treat as failed filter
            if (depth == -1) {
                return false;
            }

            // Check depth bounds
            if (depth < minDepth || depth > maxDepth) {
                return false;
            }
        }

        return true;
    }

    void filterVCF() {
        bool inputIsGzipped = isGzipped(inputFile);
        bool outputIsGzipped = isGzipped(outputFile);
        
        // Input handling
        std::ifstream input;
        gzFile gzInput = nullptr;
        
        if (inputIsGzipped) {
            gzInput = gzopen(inputFile.c_str(), "rb");
            if (!gzInput) {
                std::cerr << "Error: Cannot open compressed input file " << inputFile << std::endl;
                return;
            }
        } else {
            input.open(inputFile);
            if (!input.is_open()) {
                std::cerr << "Error: Cannot open input file " << inputFile << std::endl;
                return;
            }
        }

        // Output handling
        std::ofstream output;
        gzFile gzOutput = nullptr;
        
        if (outputIsGzipped) {
            gzOutput = gzopen(outputFile.c_str(), "wb");
            if (!gzOutput) {
                std::cerr << "Error: Cannot create compressed output file " << outputFile << std::endl;
                if (inputIsGzipped) gzclose(gzInput);
                else input.close();
                return;
            }
        } else {
            output.open(outputFile);
            if (!output.is_open()) {
                std::cerr << "Error: Cannot create output file " << outputFile << std::endl;
                if (inputIsGzipped) gzclose(gzInput);
                else input.close();
                return;
            }
        }

        std::string line;
        int totalVariants = 0;
        int passedVariants = 0;

        while (readLine(gzInput, input, line, inputIsGzipped)) {
            // Write header lines as-is
            if (line.empty() || line[0] == '#') {
                if (outputIsGzipped) {
                    gzprintf(gzOutput, "%s\n", line.c_str());
                } else {
                    output << line << std::endl;
                }
                continue;
            }

            totalVariants++;

            // Filter variant lines
            if (passesDepthFilter(line)) {
                if (outputIsGzipped) {
                    gzprintf(gzOutput, "%s\n", line.c_str());
                } else {
                    output << line << std::endl;
                }
                passedVariants++;
            }
        }

        // Close files
        if (inputIsGzipped) {
            gzclose(gzInput);
        } else {
            input.close();
        }
        
        if (outputIsGzipped) {
            gzclose(gzOutput);
        } else {
            output.close();
        }

        std::cout << "Filtering complete:\n";
        std::cout << "Total variants: " << totalVariants << std::endl;
        std::cout << "Passed variants: " << passedVariants << std::endl;
        std::cout << "Filtered variants: " << (totalVariants - passedVariants) << std::endl;
        std::cout << "Output written to: " << outputFile << std::endl;
    }

    void run(int argc, char* argv[]) {
        parseArguments(argc, argv);
        
        std::cout << "VCF Depth Filter\n";
        std::cout << "Input file: " << inputFile << std::endl;
        std::cout << "Output file: " << outputFile << std::endl;
        std::cout << "Min depth: " << minDepth << std::endl;
        std::cout << "Max depth: " << (maxDepth == INT_MAX ? "unlimited" : std::to_string(maxDepth)) << std::endl;
        std::cout << "\nProcessing...\n";

        filterVCF();
    }
};

int main(int argc, char* argv[]) {
    VCFDepthFilter filter;
    filter.run(argc, argv);
    return 0;
}
