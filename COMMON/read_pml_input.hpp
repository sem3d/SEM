#ifndef _READ_PML_INPUT_HPP_
#define _READ_PML_INPUT_HPP_

#include <string>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdlib>

struct CommonPmlSpec {
    int    n[6];             // number of element layers per side (0 = off). Order: x-, x+, y-, y+, z-, z+
    double step[6];          // total PML thickness per side (<=0 -> auto)
    double ratio[6];         // progression ratio or parameter per side
    std::string law[6];      // stretching law type: "linear", "geom", "power"
    int    npow;             // damping profile power (default 2)
    double Rc;               // target reflection coefficient (default 1e-3)
    double omegac;           // CPML parameter
    double kc;               // CPML parameter
    
    CommonPmlSpec() : npow(2), Rc(1e-3), omegac(0.0), kc(0.0) {
        for(int k=0;k<6;++k) { n[k]=0; step[k]=0.; ratio[k]=1.; law[k]="geom"; }
    }
    bool any() const { for(int k=0;k<6;++k) if(n[k]>0) return true; return false; }
};

// Custom line reader to replace dependency on mesher-specific getData_line
inline void read_line_helper(char** lineptr, size_t* n, FILE* stream) {
    if (!lineptr || !n || !stream) return;
    if (*lineptr == NULL || *n == 0) {
        *n = 512;
        *lineptr = (char*)malloc(*n);
    }
    size_t pos = 0;
    int c;
    while ((c = fgetc(stream)) != EOF) {
        if (pos + 1 >= *n) {
            *n *= 2;
            *lineptr = (char*)realloc(*lineptr, *n);
        }
        if (c == '\n') {
            break;
        }
        (*lineptr)[pos++] = c;
    }
    (*lineptr)[pos] = '\0';
}

inline bool read_common_pml_input(const char* fname, CommonPmlSpec& spec, bool is_2d) {
    FILE* f = fopen(fname, "r");
    if (!f) return false;
    
    char* buffer = NULL;
    size_t linesize = 0;
    while (true) {
        read_line_helper(&buffer, &linesize, f);
        if (feof(f) && (!buffer || buffer[0] == '\0')) break;
        if (!buffer || buffer[0] == '\0') continue;
        
        // Strip comment
        char* hash = strchr(buffer, '#');
        if (hash) *hash = '\0';
        
        char tok[64] = {0};
        if (sscanf(buffer, "%63s", tok) != 1) continue;
        
        if (!strcmp(tok, "pmlparams")) {
            double Rc = 1e-3, omegac = 0., kc = 0.;
            int npow = 2;
            // Support both 2D (4 parameters) and 3D (2 parameters) formats
            int c = sscanf(buffer, "%*s %d %lf %lf %lf", &npow, &Rc, &omegac, &kc);
            spec.npow = npow;
            spec.Rc = Rc;
            if (c >= 3) spec.omegac = omegac;
            if (c >= 4) spec.kc = kc;
            continue;
        }
        
        int s = -1;
        if (!strcmp(tok, "x-")) s = 0;
        else if (!strcmp(tok, "x+")) s = 1;
        else if (!strcmp(tok, "y-")) {
            if (is_2d) {
                printf("ERR pml.input: side '%s' is 3D-only; 2D has no y axis\n", tok);
                exit(1);
            }
            s = 2;
        }
        else if (!strcmp(tok, "y+")) {
            if (is_2d) {
                printf("ERR pml.input: side '%s' is 3D-only; 2D has no y axis\n", tok);
                exit(1);
            }
            s = 3;
        }
        else if (!strcmp(tok, "z-")) s = 4;
        else if (!strcmp(tok, "z+")) s = 5;
        else {
            printf("ERR pml.input: unknown side '%s'\n", tok);
            exit(1);
        }
        
        if (s >= 0) {
            char word1[64] = {0}, word2[64] = {0}, word3[64] = {0}, word4[64] = {0};
            int c = sscanf(buffer, "%*s %63s %63s %63s %63s", word1, word2, word3, word4);
            if (c < 1) continue;
            
            int nn = 1;
            double thick = 0., ratio = 1.;
            std::string law = "geom";
            
            // Check if word1 has a decimal point or exponent 'e'/'E' -> it is a float
            bool word1_is_float = (strchr(word1, '.') != NULL || strchr(word1, 'e') != NULL || strchr(word1, 'E') != NULL);
            
            if (word1_is_float) {
                // New format: <side> [thickness] [n_elements] [ratio] [law] (double then int)
                thick = atof(word1);
                if (c >= 2) nn = atoi(word2);
                if (c >= 3) ratio = atof(word3);
                if (c >= 4) law = word4;
            } else {
                // Old 2D format: <side> [n_elements] [thickness] [ratio] [law] (int then double)
                nn = atoi(word1);
                if (c >= 2) thick = atof(word2);
                if (c >= 3) ratio = atof(word3);
                if (c >= 4) law = word4;
            }
            
            if (thick < 0.) { printf("ERR pml.input: PML thickness for '%s' must be >= 0 (got %g)\n", tok, thick); exit(1); }
            if (nn < 1) { printf("ERR pml.input: element count for '%s' must be >= 1 (got %d)\n", tok, nn); exit(1); }
            if (ratio <= 0.) { printf("ERR pml.input: grading ratio for '%s' must be > 0 (got %g)\n", tok, ratio); exit(1); }
            
            spec.n[s] = nn;
            spec.step[s] = thick;
            spec.ratio[s] = ratio;
            spec.law[s] = law;
        }
    }
    if (buffer) free(buffer);
    fclose(f);
    return true;
}

#endif
