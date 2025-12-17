#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <direct.h>
#include <io.h>
#include <float.h>

// Maximum file definitions
#define MAX_FILES 64
#define MAX_SELECTED_FILES 20
#define PI 3.14159265358979323846

// Material database
#define NUM_MATERIALS 13  // Increased to 13 materials
#define MAX_CLUSTERS 10
#define MAX_ITERATIONS 100

typedef struct {
    char name[50];
    double pe_center;
    double zeff_center;
    int r, g, b; // Color values
    double alpha; // Li content (0-1)
    double beta;  // MnO2 content (0-1)
} Material;

// Cluster center
typedef struct {
    double pe_center;
    double zeff_center;
    int count;
    int assigned_material;
} ClusterCenter;

// Function to check if point is inside circle
int func(int x, int y, double ox, double oy, double r) {
    double d = sqrt(pow(x - ox, 2) + pow(y - oy, 2));
    return d <= r ? 1 : 0;
}

int main(void) {
    int i, j, k;
    FILE* fp = NULL;
    FILE* cp = NULL;
    errno_t error;
    double en[4];
    double** z = NULL;  // Dynamic 2D array
    double* a = NULL;   // Image data array
    int file_count = 0;  // Processed file count

    // Create results directory
    _mkdir("result");

    char filename[256];
    char line[256];
    char all_file_list[MAX_FILES][256] = { {0} };  // All files list
    double all_F[MAX_FILES] = { 0 };
    double all_G[MAX_FILES] = { 0 };
    double all_energy[MAX_FILES] = { 0 };  // Store energy for each file
    int total_files = 0;  // Total files

    // Read energy.txt file
    double base_energy[20] = {0};
    error = fopen_s(&fp, "energy.txt", "r");
    if (error != 0) {
        printf("Cannot open energy.txt\n");
        return 1;
    }
    for (i = 0; i < 20; i++) {
        if (fscanf_s(fp, "%lf", &base_energy[i]) != 1) {
            printf("Failed to read energy.txt\n");
            fclose(fp);
            return 1;
        }
        base_energy[i] *= 1000.0; // Convert to keV
    }
    fclose(fp);
    printf("Loaded 20 energy values from energy.txt\n");

    // Open setting1.txt
    error = fopen_s(&cp, "setting1.txt", "r");
    if (error != 0) {
        printf("Cannot open setting1.txt\n");
        return 1;
    }

    // Skip comment lines and read parameters
    for (k = 0; k < 4; k++) {
        if (fgets(line, sizeof(line), cp) == NULL) {
            printf("Settings file format error\n");
            fclose(cp);
            return 1;
        }
        // Skip comment lines and empty lines
        while (line[0] == '#' || line[0] == '\n') {
            if (fgets(line, sizeof(line), cp) == NULL) {
                printf("Settings file format error\n");
                fclose(cp);
                return 1;
            }
        }
        sscanf_s(line, "%lf", &en[k]);
    }

    // Set parameters
    double ox = en[1];
    double oy = en[2];
    double r = en[3];
    printf("Circle center: (%.2f, %.2f), radius: %.2f\n", ox, oy, r);
    
    // Read all file list, F, G values, use energy values from energy.txt
    int energy_index = 0;
    while (fscanf_s(cp, "%255s", filename, (unsigned)_countof(filename)) != EOF) {
        if (fscanf_s(cp, "%lf", &all_F[total_files]) != 1) break;
        if (fscanf_s(cp, "%lf", &all_G[total_files]) != 1) break;
        
        // Skip original energy value in file
        double dummy;
        if (fscanf_s(cp, "%lf", &dummy) != 1) break;
        
        // Use energy from energy.txt
        if (energy_index < 20) {
            all_energy[total_files] = base_energy[energy_index];
            energy_index++;
        } else {
            all_energy[total_files] = 0; // If more than 20 files, set to 0
        }

        strcpy_s(all_file_list[total_files], sizeof(all_file_list[total_files]), filename);
        printf("File %d: %s, F=%.2E, G=%.2E, Energy=%.2f keV\n", 
               total_files + 1, all_file_list[total_files], 
               all_F[total_files], all_G[total_files], all_energy[total_files]);
        total_files++;

        if (total_files >= MAX_FILES) break;
    }
    fclose(cp);
    printf("Total files loaded: %d\n", total_files);

    // Let user select number of files and specific files
    int num_files = 0;
    printf("Select number of files to process (0=all, 1-20): ");
    scanf_s("%d", &num_files);

    // Validate input
    if (num_files == 0) {
        // Select all files
        num_files = total_files;
        if (num_files > MAX_SELECTED_FILES) {
            printf("Note: Total files (%d) exceeds maximum selectable (%d), will select only first %d files.\n", total_files, MAX_SELECTED_FILES, MAX_SELECTED_FILES);
            num_files = MAX_SELECTED_FILES;
        }
        printf("Selected all %d files\n", num_files);
    }
    else if (num_files < 0 || num_files > MAX_SELECTED_FILES) {
        printf("Error: Invalid number of files. Please enter 0 (all) or an integer between 1 and %d.\n", MAX_SELECTED_FILES);
        return 1;
    }

    int selected_files[MAX_SELECTED_FILES] = { 0 };
    if (num_files == total_files || num_files == MAX_SELECTED_FILES) {
        // Automatically select first num_files files
        for (i = 0; i < num_files; i++) {
            selected_files[i] = i; // Index from 0
        }
    }
    else {
        printf("Select %d files (enter numbers, separated by spaces): ", num_files);
        for (i = 0; i < num_files; i++) {
            scanf_s("%d", &selected_files[i]);
            // Validate input
            if (selected_files[i] < 1 || selected_files[i] > total_files) {
                printf("Error: Invalid file selection. Please enter numbers between 1 and %d.\n", total_files);
                return 1;
            }
            selected_files[i]--; // Convert to 0-based index
        }
    }

    // Set final selected files
    file_count = num_files;
    char file_list[MAX_SELECTED_FILES][256] = { 0 };
    double F[MAX_SELECTED_FILES] = { 0 };
    double G[MAX_SELECTED_FILES] = { 0 };
    double f[MAX_SELECTED_FILES] = { 0 };
    double g[MAX_SELECTED_FILES] = { 0 };
    double selected_energy[MAX_SELECTED_FILES] = {0}; // Energy for selected files

    for (i = 0; i < file_count; i++) {
        int idx = selected_files[i];
        strcpy_s(file_list[i], sizeof(file_list[i]), all_file_list[idx]);
        F[i] = all_F[idx];
        G[i] = all_G[idx];
        selected_energy[i] = all_energy[idx]; // Store energy value
        f[i] = F[i] * 1e20;
        g[i] = G[i] * 1e20;
        printf("Selected file %d: %s, F=%.2E, G=%.2E, Energy=%.2f keV\n", 
               i + 1, file_list[i], F[i], G[i], selected_energy[i]);
    }

    // Ask user for current processing state
    char state[20];
    printf("Select processing state (enter 'fresh' or 'discharged'): ");
    scanf_s("%19s", state, (unsigned)_countof(state));

    // Validate input
    if (strcmp(state, "fresh") != 0 && strcmp(state, "discharged") != 0) {
        printf("Error: Invalid state selection. Please enter 'fresh' or 'discharged'.\n");
        return 1;
    }

    // Determine image dimensions based on file size
    const int width = 105;
    const int height = 105;
    const int total_pixels = width * height;
    const long expected_size = total_pixels * sizeof(double);
    printf("Image dimensions: %d x %d = %d pixels\n", width, height, total_pixels);
    printf("Expected file size: %ld bytes\n", expected_size);

    // Allocate memory
    a = (double*)malloc(total_pixels * sizeof(double));
    if (!a) {
        printf("Memory allocation failed: a\n");
        return 1;
    }

    // Allocate memory for z array
    z = (double**)malloc(file_count * sizeof(double*));
    if (!z) {
        printf("Memory allocation failed: z\n");
        free(a);
        return 1;
    }

    for (k = 0; k < file_count; k++) {
        z[k] = (double*)malloc(total_pixels * sizeof(double));
        if (!z[k]) {
            printf("Memory allocation failed: z[%d]\n", k);
            for (int i = 0; i < k; i++) free(z[i]);
            free(z);
            free(a);
            return 1;
        }
    }

    // Initialize other arrays
    double* sum_m = (double*)calloc(total_pixels, sizeof(double));
    double* Pe = NULL;
    double* Zeff = NULL;
    double* Pe_re = NULL;

    // Only allocate Pe, Zeff, Pe_re if more than 1 file
    if (file_count > 1) {
        Pe = (double*)malloc(total_pixels * sizeof(double));
        Zeff = (double*)malloc(total_pixels * sizeof(double));
        Pe_re = (double*)malloc(total_pixels * sizeof(double));

        if (!Pe || !Zeff || !Pe_re) {
            printf("Memory allocation failed\n");
            free(a);
            free(sum_m);
            for (k = 0; k < file_count; k++) free(z[k]);
            free(z);
            return 1;
        }
    }

    // Define region boundaries
    int region_min_x, region_max_x, region_min_y, region_max_y;
    if (strcmp(state, "fresh") == 0) {
        // Fresh state region definition
        region_min_x = 44;
        region_max_x = 65;
        region_min_y = 22;
        region_max_y = 85;
        printf("Fresh state region: x[%d-%d], y[%d-%d]\n", 
               region_min_x, region_max_x, region_min_y, region_max_y);
    }
    else { // discharged
        // Discharged state region definition
        region_min_x = 45;
        region_max_x = 67;
        region_min_y = 30;
        region_max_y = 83;
        printf("Discharged state region: x[%d-%d], y[%d-%d]\n", 
               region_min_x, region_max_x, region_min_y, region_max_y);
    }

    // Define cathode region - fixed range (based on paper Figure 2)
    int cathode_min_x = 47; 
    int cathode_max_x = 61;
    int cathode_min_y = 22;
    int cathode_max_y = 85;
    
    // Adjust discharged state cathode region (based on paper Figure 2b)
    if (strcmp(state, "discharged") == 0) {
        cathode_min_x = 49;
        cathode_max_x = 60;
        cathode_min_y = 31;
        cathode_max_y = 82;
    }
    
    printf("Cathode region: x[%d-%d], y[%d-%d]\n", 
           cathode_min_x, cathode_max_x, cathode_min_y, cathode_max_y);

    // Store region average linear attenuation coefficient for each file
    double region_avg[MAX_SELECTED_FILES] = {0};

    // Process files
    for (k = 0; k < file_count; k++) {
        error = fopen_s(&fp, file_list[k], "rb");
        if (error != 0) {
            printf("Cannot open file %s\n", file_list[k]);
            continue;
        }

        // Check file size
        _fseeki64(fp, 0, SEEK_END);
        __int64 file_size = _ftelli64(fp);
        _fseeki64(fp, 0, SEEK_SET);

        if (file_size != expected_size) {
            printf("File size mismatch: %s (actual: %lld, expected: %ld)\n",
                file_list[k], file_size, expected_size);
            fclose(fp);
            continue;
        }

        // Read data
        size_t read = fread(a, sizeof(double), total_pixels, fp);
        fclose(fp);

        if (read != total_pixels) {
            printf("File read incomplete: %s (read: %zu, expected: %d)\n",
                file_list[k], read, total_pixels);
            continue;
        }

        // Copy data to z array
        memcpy(z[k], a, total_pixels * sizeof(double));

        // Process pixels inside circle
        double Nin = 0;
        double m = 0;
        for (j = 0; j < height; j++) {
            for (i = 0; i < width; i++) {
                int in = func(i, j, ox, oy, r);
                if (in == 1) {
                    Nin += 1;
                    m += a[j * width + i];
                }
            }
        }

        if (Nin > 0) {
            m = m / Nin;
            printf("File %d circle average: %.4f\n", k + 1, m);
        }

        // Process specified region pixels
        double region_sum = 0;
        int region_count = 0;
        for (j = region_min_y; j <= region_max_y; j++) {
            for (i = region_min_x; i <= region_max_x; i++) {
                region_sum += a[j * width + i];
                region_count++;
            }
        }

        if (region_count > 0) {
            region_avg[k] = region_sum / region_count;
            printf("File %d region average: %.4f\n", k + 1, region_avg[k]);
        }
        else {
            region_avg[k] = 0;
            printf("File %d region has no valid pixels\n", k + 1);
        }

        for (i = 0; i < total_pixels; i++) {
            sum_m[i] += a[i];
        }
    }

    // Create energy-linear attenuation coefficient CSV file
    char mu_energy_filename[256];
    sprintf_s(mu_energy_filename, sizeof(mu_energy_filename), "result\\mu_vs_energy_%s.csv", state);
    error = fopen_s(&fp, mu_energy_filename, "w");
    if (error == 0) {
        fprintf(fp, "energy,mu_avg\n");
        for (k = 0; k < file_count; k++) {
            fprintf(fp, "%.6f,%.6f\n", selected_energy[k], region_avg[k]);
        }
        fclose(fp);
        printf("Created energy-linear attenuation file: %s\n", mu_energy_filename);
    }
    else {
        printf("Failed to create energy-linear attenuation file\n");
    }

    // Calculate average
    double* z_final_average = (double*)calloc(total_pixels, sizeof(double));
    for (i = 0; i < total_pixels; i++) {
        double sum_z = 0;
        for (k = 0; k < file_count; k++) {
            sum_z += z[k][i];
        }
        z_final_average[i] = sum_z / file_count;
    }
    
    // Calculate high-energy region average linear attenuation coefficient image
    double* high_energy_avg = (double*)calloc(total_pixels, sizeof(double));
    int high_energy_file_count = 0;

    // High-energy threshold (above 90keV)
    const double high_energy_threshold = 90.0;

    for (k = 0; k < file_count; k++) {
        if (selected_energy[k] >= high_energy_threshold) {
            for (i = 0; i < total_pixels; i++) {
                high_energy_avg[i] += z[k][i];
            }
            high_energy_file_count++;
        }
    }

    if (high_energy_file_count > 0) {
        for (i = 0; i < total_pixels; i++) {
            high_energy_avg[i] /= high_energy_file_count;
        }
        printf("High-energy average image calculated (using %d files)\n", high_energy_file_count);
    } else {
        printf("Warning: No high-energy files (>=%.1f keV) for average image calculation\n", high_energy_threshold);
        // If no high-energy files, fill with 0
        memset(high_energy_avg, 0, total_pixels * sizeof(double));
    }
    
    // Calculate high-energy region average linear attenuation coefficient in cathode region
    double cathode_avg = 0.0;
    int cathode_count = 0;
    for (j = cathode_min_y; j <= cathode_max_y; j++) {
        for (i = cathode_min_x; i <= cathode_max_x; i++) {
            int idx = j * width + i;
            cathode_avg += high_energy_avg[idx];
            cathode_count++;
        }
    }
    
    if (cathode_count > 0) {
        cathode_avg /= cathode_count;
        printf("Cathode region high-energy average linear attenuation: %.6f cm^-1\n", cathode_avg);
    } else {
        cathode_avg = 0.0;
        printf("Warning: Cathode region has no valid pixels\n");
    }
    
    // Save cathode region average as calibration value
    char cal_filename[256];
    sprintf_s(cal_filename, sizeof(cal_filename), "result\\calibration_%s.txt", state);
    error = fopen_s(&fp, cal_filename, "w");
    if (error == 0) {
        fprintf(fp, "%.6f\n", cathode_avg);
        fclose(fp);
        printf("Saved calibration value to: %s\n", cal_filename);
    } else {
        printf("Failed to save calibration value\n");
    }
    
    // Calculate SOC image (using actual range of current state)
    double* soc_image = (double*)malloc(total_pixels * sizeof(double));
    
    // Calculate actual linear attenuation coefficient range in cathode region for current state
    double min_mu = 1e10;
    double max_mu = -1e10;
    
    for (j = cathode_min_y; j <= cathode_max_y; j++) {
        for (i = cathode_min_x; i <= cathode_max_x; i++) {
            int idx = j * width + i;
            double mu = high_energy_avg[idx];
            if (mu < min_mu) min_mu = mu;
            if (mu > max_mu) max_mu = mu;
        }
    }
    
    // Ensure valid range
    if (max_mu <= min_mu) {
        max_mu = min_mu + 0.001; // Avoid division by zero
        printf("Warning: mu_max <= mu_min, range adjusted\n");
    }
    
    printf("Using current state range for SOC: mu_min=%.6f, mu_max=%.6f\n", min_mu, max_mu);
    
    // Calculate SOC
    for (i = 0; i < total_pixels; i++) {
        double mu = high_energy_avg[i];
        if (mu <= min_mu) {
            soc_image[i] = 0.0;
        } else if (mu >= max_mu) {
            soc_image[i] = 100.0;
        } else {
            soc_image[i] = (mu - min_mu) / (max_mu - min_mu) * 100.0;
        }
    }

    // Save SOC image
    char soc_filename[256];
    sprintf_s(soc_filename, sizeof(soc_filename), "result\\SOC_image_%s.csv", state);
    error = fopen_s(&fp, soc_filename, "w");
    if (error == 0) {
        fprintf(fp, "DataFormat,2\nmemo\nX (pixel), Y (pixel), Z\n");
        for (j = 0; j < height; j++) {
            for (i = 0; i < width; i++) {
                int idx = j * width + i;
                fprintf(fp, "%d,%d,%.2f\n", i, j, soc_image[idx]);
            }
        }
        fclose(fp);
        printf("Created SOC image file: %s\n", soc_filename);
    } else {
        printf("Failed to create SOC image file\n");
    }
    
    // Save cathode region SOC image
    char cathode_soc_filename[256];
    sprintf_s(cathode_soc_filename, sizeof(cathode_soc_filename), "result\\Cathode_SOC_%s.csv", state);
    error = fopen_s(&fp, cathode_soc_filename, "w");
    if (error == 0) {
        fprintf(fp, "DataFormat,2\nmemo\nX (pixel), Y (pixel), Z\n");
        for (j = cathode_min_y; j <= cathode_max_y; j++) {
            for (i = cathode_min_x; i <= cathode_max_x; i++) {
                int idx = j * width + i;
                fprintf(fp, "%d,%d,%.2f\n", i, j, soc_image[idx]);
            }
        }
        fclose(fp);
        printf("Created cathode region SOC file: %s\n", cathode_soc_filename);
    } else {
        printf("Failed to create cathode region SOC file\n");
    }
    
    // Save high-energy average image
    char high_energy_filename[256];
    sprintf_s(high_energy_filename, sizeof(high_energy_filename), "result\\HighEnergyMu_%s.raw", state);
    error = fopen_s(&fp, high_energy_filename, "wb");
    if (error == 0) {
        fwrite(high_energy_avg, sizeof(double), total_pixels, fp);
        fclose(fp);
        printf("Created high-energy average image file: %s\n", high_energy_filename);
    } else {
        printf("Failed to create high-energy average image file\n");
    }

    // Calculate Pe, Zeff, Pe_re (only when more than 1 file)
    if (file_count > 1) {
        printf("Calculating Pe, Zeff, and Pe_re...\n");

        // Pre-calculate sum_f and sum_g
        double sum_f = 0.0, sum_g = 0.0;
        for (k = 0; k < file_count; k++) {
            sum_f += f[k];
            sum_g += g[k];
        }

        for (i = 0; i < total_pixels; i++) {
            // 2 files use special formula
            if (file_count == 2) {
                double z1 = z[0][i];
                double z2 = z[1][i];
                double f1 = f[0];
                double f2 = f[1];
                double g1 = g[0];
                double g2 = g[1];

                // Prevent division by zero
                if (fabs(f1 * g2 - f2 * g1) < 1e-20) {
                    Pe[i] = 0.0;
                    Zeff[i] = 0.0;
                    Pe_re[i] = 0.0;
                    continue;
                }

                // Calculate a0 and a1
                double a0 = (z1 * g2 - z2 * g1) / (f1 * g2 - f2 * g1);
                double a1 = (f1 * z2 - f2 * z1) / (f1 * g2 - f2 * g1);

                // Calculate Pe and Zeff
                Pe[i] = a0 * 1e20;
                if (fabs(a0) > 1e-20) {
                    Zeff[i] = pow(fabs(a1 / a0), 0.25);
                }
                else {
                    Zeff[i] = 0.0;
                }

                // Calculate Pe_re
                double denom = (a1 * sum_f / a0) + sum_g;
                if (fabs(denom) > 1e-20) {
                    Pe_re[i] = sum_m[i] / denom * 1e20;
                }
                else {
                    Pe_re[i] = 0.0;
                }
            }
            // 3 or more files use least squares method
            else {
                double A00 = 0, A01 = 0, A02 = 0, A11 = 0, A12 = 0;
                for (int j = 0; j < file_count; j++) {
                    double fg_ratio = f[j] / g[j];
                    double z_g = z[j][i] / g[j];

                    A00 += 1.0;
                    A01 += fg_ratio;
                    A02 += z_g;
                    A11 += fg_ratio * fg_ratio;
                    A12 += fg_ratio * z_g;
                }

                double denom = A00 * A11 - A01 * A01;
                if (fabs(denom) > 1e-10) {
                    double a0 = (A02 * A11 - A01 * A12) / denom;
                    double a1 = (A00 * A12 - A01 * A02) / denom;

                    Pe[i] = a0 * 1e20;
                    if (fabs(a0) > 1e-20) {
                        Zeff[i] = pow(fabs(a1 / a0), 0.25);
                    }
                    else {
                        Zeff[i] = 0.0;
                    }

                    // Calculate Pe_re
                    double denom2 = (a1 * sum_f / a0) + sum_g;
                    if (fabs(denom2) > 1e-20) {
                        Pe_re[i] = sum_m[i] / denom2 * 1e20;
                    }
                    else {
                        Pe_re[i] = 0.0;
                    }
                }
                else {
                    Pe[i] = 0.0;
                    Zeff[i] = 0.0;
                    Pe_re[i] = 0.0;
                }
            }
        }

        // Output Pe and Zeff statistics for debugging
        double pe_min = DBL_MAX, pe_max = -DBL_MAX, pe_sum = 0;
        double zeff_min = DBL_MAX, zeff_max = -DBL_MAX, zeff_sum = 0;
        int valid_count = 0;
        
        for (i = 0; i < total_pixels; i++) {
            if (Pe[i] > 0 && Pe[i] < 1e30 && Zeff[i] > 0 && Zeff[i] < 100) {
                if (Pe[i] < pe_min) pe_min = Pe[i];
                if (Pe[i] > pe_max) pe_max = Pe[i];
                if (Zeff[i] < zeff_min) zeff_min = Zeff[i];
                if (Zeff[i] > zeff_max) zeff_max = Zeff[i];
                pe_sum += Pe[i];
                zeff_sum += Zeff[i];
                valid_count++;
            }
        }
        
        if (valid_count > 0) {
            printf("Pe statistics: min=%.2e, max=%.2e, avg=%.2e\n", pe_min, pe_max, pe_sum/valid_count);
            printf("Zeff statistics: min=%.2f, max=%.2f, avg=%.2f\n", zeff_min, zeff_max, zeff_sum/valid_count);
        }
    }

    // Material identification section - using complete lithiation material database
    if (file_count > 1) {
        printf("Performing material identification...\n");
        
        // Complete material database - includes LixMnO2 with different lithiation degrees
        Material materials[NUM_MATERIALS] = {
            // Main active materials
            {"LiMetal",      0.65e23,   3.0,  255, 0,   0,    1.0, 0.0},   // Red - Lithium metal (anode)
            {"MnO2",         1.20e23,   8.5,  0,   0,   255,  0.0, 1.0},   // Blue - Manganese dioxide (cathode fresh)
            
            // LixMnO2 compounds with different lithiation degrees
            {"Li0.1MnO2",    1.15e23,   8.2,  0,   50,  255,  0.1, 0.9},   // Light blue - low lithiation
            {"Li0.2MnO2",    1.10e23,   7.9,  0,   100, 255,  0.2, 0.8},   // Blue-green
            {"Li0.3MnO2",    1.05e23,   7.6,  0,   150, 255,  0.3, 0.7},   // Cyan-blue
            {"Li0.4MnO2",    1.00e23,   7.3,  0,   200, 255,  0.4, 0.6},   // Cyan-green
            {"Li0.5MnO2",    0.95e23,   7.0,  0,   255, 255,  0.5, 0.5},   // Cyan
            {"Li0.6MnO2",    0.90e23,   6.7,  0,   255, 200,  0.6, 0.4},   // Cyan-yellow
            {"Li0.7MnO2",    0.85e23,   6.4,  0,   255, 150,  0.7, 0.3},   // Yellow-green
            {"Li0.8MnO2",    0.80e23,   6.1,  0,   255, 100,  0.8, 0.2},   // Orange-yellow
            {"Li0.9MnO2",    0.75e23,   5.8,  0,   255, 50,   0.9, 0.1},   // Orange
            {"LiMnO2",       0.70e23,   5.5,  0,   255, 0,    1.0, 0.0},   // Red - fully lithiated
        };

        // Define material identification region
        int material_region_min_x, material_region_max_x, material_region_min_y, material_region_max_y;
        
        if (strcmp(state, "fresh") == 0) {
            material_region_min_x = 44;
            material_region_max_x = 65;
            material_region_min_y = 22;
            material_region_max_y = 85;
        } else {
            material_region_min_x = 45;
            material_region_max_x = 67;
            material_region_min_y = 30;
            material_region_max_y = 83;
        }
        
        printf("Material identification region: x[%d-%d], y[%d-%d]\n", 
               material_region_min_x, material_region_max_x, 
               material_region_min_y, material_region_max_y);
        
        // First analyze actual data range to adjust material database
        double min_pe = DBL_MAX, max_pe = -DBL_MAX;
        double min_zeff = DBL_MAX, max_zeff = -DBL_MAX;
        
        for (j = material_region_min_y; j <= material_region_max_y; j++) {
            for (i = material_region_min_x; i <= material_region_max_x; i++) {
                int idx = j * width + i;
                if (Pe[idx] > 0 && Pe[idx] < 1e30 && Zeff[idx] > 0 && Zeff[idx] < 100) {
                    if (Pe[idx] < min_pe) min_pe = Pe[idx];
                    if (Pe[idx] > max_pe) max_pe = Pe[idx];
                    if (Zeff[idx] < min_zeff) min_zeff = Zeff[idx];
                    if (Zeff[idx] > max_zeff) max_zeff = Zeff[idx];
                }
            }
        }
        
        printf("Actual data range: Pe=[%.2e, %.2e], Zeff=[%.2f, %.2f]\n", min_pe, max_pe, min_zeff, max_zeff);
        
        // Adjust material database parameters based on actual data range
        if (strcmp(state, "fresh") == 0) {
            // Fresh state: adjust main material parameters
            printf("Adjusting fresh state material parameters...\n");
            
            // Calculate actual Pe and Zeff in cathode region
            double cathode_pe_sum = 0.0, cathode_zeff_sum = 0.0;
            int cathode_count = 0;
            for (j = cathode_min_y; j <= cathode_max_y; j++) {
                for (i = cathode_min_x; i <= cathode_max_x; i++) {
                    int idx = j * width + i;
                    if (Pe[idx] > 0 && Pe[idx] < 1e30 && Zeff[idx] > 0 && Zeff[idx] < 100) {
                        cathode_pe_sum += Pe[idx];
                        cathode_zeff_sum += Zeff[idx];
                        cathode_count++;
                    }
                }
            }
            
            if (cathode_count > 0) {
                materials[1].pe_center = cathode_pe_sum / cathode_count;  // MnO2
                materials[1].zeff_center = cathode_zeff_sum / cathode_count;
                printf("Adjusted MnO2 parameters: Pe=%.2e, Zeff=%.2f\n", materials[1].pe_center, materials[1].zeff_center);
            }
            
            // Calculate actual Pe and Zeff in anode region (fresh state only)
            double anode_pe_sum = 0.0, anode_zeff_sum = 0.0;
            int anode_count = 0;
            int anode_region_min_x = 63, anode_region_max_x = 65;
            int anode_region_min_y = 22, anode_region_max_y = 85;
            
            for (j = anode_region_min_y; j <= anode_region_max_y; j++) {
                for (i = anode_region_min_x; i <= anode_region_max_x; i++) {
                    int idx = j * width + i;
                    if (Pe[idx] > 0 && Pe[idx] < 1e30 && Zeff[idx] > 0 && Zeff[idx] < 100) {
                        anode_pe_sum += Pe[idx];
                        anode_zeff_sum += Zeff[idx];
                        anode_count++;
                    }
                }
            }
            
            if (anode_count > 0) {
                materials[0].pe_center = anode_pe_sum / anode_count;  // LiMetal
                materials[0].zeff_center = anode_zeff_sum / anode_count;
                printf("Adjusted LiMetal parameters: Pe=%.2e, Zeff=%.2f\n", materials[0].pe_center, materials[0].zeff_center);
            }
        } else {
            // Discharged state: adjust LixMnO2 series parameters based on actual data
            printf("Adjusting discharged state material parameters...\n");
            
            // Calculate actual Pe and Zeff range in cathode region
            double cathode_pe_min = DBL_MAX, cathode_pe_max = -DBL_MAX;
            double cathode_zeff_min = DBL_MAX, cathode_zeff_max = -DBL_MAX;
            int cathode_count = 0;
            
            for (j = cathode_min_y; j <= cathode_max_y; j++) {
                for (i = cathode_min_x; i <= cathode_max_x; i++) {
                    int idx = j * width + i;
                    if (Pe[idx] > 0 && Pe[idx] < 1e30 && Zeff[idx] > 0 && Zeff[idx] < 100) {
                        if (Pe[idx] < cathode_pe_min) cathode_pe_min = Pe[idx];
                        if (Pe[idx] > cathode_pe_max) cathode_pe_max = Pe[idx];
                        if (Zeff[idx] < cathode_zeff_min) cathode_zeff_min = Zeff[idx];
                        if (Zeff[idx] > cathode_zeff_max) cathode_zeff_max = Zeff[idx];
                        cathode_count++;
                    }
                }
            }
            
            if (cathode_count > 0) {
                // Recalculate Pe and Zeff for LixMnO2 series based on actual range
                double pe_range = cathode_pe_max - cathode_pe_min;
                double zeff_range = cathode_zeff_max - cathode_zeff_min;
                
                for (int m = 2; m < NUM_MATERIALS; m++) {
                    double alpha = materials[m].alpha;
                    // Linear interpolation for Pe and Zeff
                    materials[m].pe_center = cathode_pe_min + alpha * pe_range;
                    materials[m].zeff_center = cathode_zeff_min + alpha * zeff_range;
                    printf("Adjusted %s parameters: Pe=%.2e, Zeff=%.2f, alpha=%.1f\n", 
                           materials[m].name, materials[m].pe_center, materials[m].zeff_center, alpha);
                }
            }
        }
        
        // Create material distribution image
        int* material_map = (int*)malloc(total_pixels * sizeof(int));
        int* material_count = (int*)calloc(NUM_MATERIALS, sizeof(int));
        
        // Initialize material_map as unknown
        for (i = 0; i < total_pixels; i++) {
            material_map[i] = -1; // Use -1 for unknown
        }
        
        // Enhanced matching algorithm: detect all materials
        int material_pixel_count = 0;
        for (j = material_region_min_y; j <= material_region_max_y; j++) {
            for (i = material_region_min_x; i <= material_region_max_x; i++) {
                int idx = j * width + i;
                double current_pe = Pe[idx];
                double current_zeff = Zeff[idx];
                double current_soc = soc_image[idx];
                
                // Check value validity
                if (current_pe < 1e10 || current_pe > 1e30 || current_zeff < 0.1 || current_zeff > 50) {
                    material_map[idx] = -1;
                    continue;
                }
                
                // Region priority identification
                int is_in_cathode = (i >= cathode_min_x && i <= cathode_max_x && j >= cathode_min_y && j <= cathode_max_y);
                int is_in_anode = (i >= 63 && i <= 65 && j >= 22 && j <= 85); // Fresh state anode region
                
                // Directly assign material type based on region and state
                int assigned_material = -1;
                
                if (strcmp(state, "fresh") == 0) {
                    if (is_in_cathode) {
                        assigned_material = 1; // MnO2
                    } else if (is_in_anode) {
                        assigned_material = 0; // LiMetal
                    }
                } else {
                    if (is_in_cathode) {
                        // Discharged state: use distance matching for LixMnO2 series
                        double min_distance = DBL_MAX;
                        for (int m = 2; m < NUM_MATERIALS; m++) {
                            // Calculate normalized Euclidean distance
                            double pe_distance = (current_pe - materials[m].pe_center) / 1e23;
                            double zeff_distance = (current_zeff - materials[m].zeff_center) / 5.0;
                            double distance = sqrt(pe_distance * pe_distance + zeff_distance * zeff_distance);
                            
                            if (distance < min_distance) {
                                min_distance = distance;
                                assigned_material = m;
                            }
                        }
                        
                        // Use loose matching threshold
                        if (min_distance > 2.0) {
                            assigned_material = -1; // unknown
                        }
                    }
                }
                
                // If not assigned by region, use distance matching for all materials
                if (assigned_material == -1) {
                    double min_distance = DBL_MAX;
                    
                    for (int m = 0; m < NUM_MATERIALS; m++) {
                        // Calculate normalized Euclidean distance
                        double pe_distance = (current_pe - materials[m].pe_center) / 1e23;
                        double zeff_distance = (current_zeff - materials[m].zeff_center) / 5.0;
                        double distance = sqrt(pe_distance * pe_distance + zeff_distance * zeff_distance);
                        
                        if (distance < min_distance) {
                            min_distance = distance;
                            assigned_material = m;
                        }
                    }
                    
                    // Use loose matching threshold
                    if (min_distance > 3.0) {
                        assigned_material = -1; // unknown
                    }
                }
                
                material_map[idx] = assigned_material;
                if (assigned_material != -1) {
                    material_count[assigned_material]++;
                }
                material_pixel_count++;
            }
        }
        
        // Output material distribution results
        char material_csv_filename[256];
        sprintf_s(material_csv_filename, sizeof(material_csv_filename), "result\\MaterialDistribution_%s.csv", state);
        error = fopen_s(&fp, material_csv_filename, "w");
        if (error == 0) {
            fprintf(fp, "X,Y,Material,Pe,Zeff,SOC,Alpha,Beta,Pe_Ref,Zeff_Ref\n");
            for (j = material_region_min_y; j <= material_region_max_y; j++) {
                for (i = material_region_min_x; i <= material_region_max_x; i++) {
                    int idx = j * width + i;
                    int material_idx = material_map[idx];
                    if (material_idx >= 0 && material_idx < NUM_MATERIALS) {
                        fprintf(fp, "%d,%d,%s,%.6e,%.6f,%.2f,%.2f,%.2f,%.6e,%.6f\n", 
                            i, j, materials[material_idx].name, Pe[idx], Zeff[idx], soc_image[idx],
                            materials[material_idx].alpha, materials[material_idx].beta,
                            materials[material_idx].pe_center, materials[material_idx].zeff_center);
                    } else {
                        fprintf(fp, "%d,%d,Unknown,%.6e,%.6f,%.2f,0.0,0.0,0.0,0.0\n", 
                            i, j, Pe[idx], Zeff[idx], soc_image[idx]);
                    }
                }
            }
            fclose(fp);
            printf("Created material distribution file: %s\n", material_csv_filename);
        }
        
        // Output material statistics
        char material_stats_filename[256];
        sprintf_s(material_stats_filename, sizeof(material_stats_filename), "result\\MaterialStatistics_%s.csv", state);
        error = fopen_s(&fp, material_stats_filename, "w");
        if (error == 0) {
            fprintf(fp, "Material,Pe_Ref,Zeff_Ref,Alpha,Beta,Count,Percentage\n");
            int total_matched = 0;
            for (int m = 0; m < NUM_MATERIALS; m++) {
                total_matched += material_count[m];
            }
            
            for (int m = 0; m < NUM_MATERIALS; m++) {
                if (material_count[m] > 0) {
                    double percentage = (double)material_count[m] / total_matched * 100.0;
                    fprintf(fp, "%s,%.6e,%.6f,%.2f,%.2f,%d,%.2f\n", 
                        materials[m].name, materials[m].pe_center, materials[m].zeff_center,
                        materials[m].alpha, materials[m].beta, material_count[m], percentage);
                }
            }
            fclose(fp);
            printf("Created material statistics file: %s\n", material_stats_filename);
            
            // Console output for all material distributions
            printf("\n=== Material Identification Results ===\n");
            printf("Identification results:\n");
            for (int m = 0; m < NUM_MATERIALS; m++) {
                if (material_count[m] > 0) {
                    double percentage = (double)material_count[m] / total_matched * 100.0;
                    printf("  %s (α=%.1f, β=%.1f): %.1f%% (%d pixels)\n", 
                           materials[m].name, materials[m].alpha, materials[m].beta, 
                           percentage, material_count[m]);
                }
            }
        }
        
        // Output material distribution image (PPM format)
        char material_ppm_filename[256];
        sprintf_s(material_ppm_filename, sizeof(material_ppm_filename), "result\\MaterialMap_%s.ppm", state);
        error = fopen_s(&fp, material_ppm_filename, "w");
        if (error == 0) {
            fprintf(fp, "P3\n%d %d\n255\n", width, height);
            for (j = 0; j < height; j++) {
                for (i = 0; i < width; i++) {
                    int idx = j * width + i;
                    int material_idx = material_map[idx];
                    if (material_idx >= 0 && material_idx < NUM_MATERIALS) {
                        fprintf(fp, "%d %d %d ", 
                            materials[material_idx].r, 
                            materials[material_idx].g, 
                            materials[material_idx].b);
                    } else {
                        // Unknown material shown as gray
                        fprintf(fp, "200 200 200 ");
                    }
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("Created material distribution image: %s\n", material_ppm_filename);
        }
        
        // Discharged state special processing: analyze lithiation distribution in cathode region
        if (strcmp(state, "discharged") == 0) {
            printf("\n=== Discharged State Cathode Lithiation Analysis ===\n");
            
            // Create lithiation analysis file
            char lithiation_filename[256];
            sprintf_s(lithiation_filename, sizeof(lithiation_filename), "result\\LithiationAnalysis_%s.csv", state);
            error = fopen_s(&fp, lithiation_filename, "w");
            if (error == 0) {
                fprintf(fp, "X,Y,Pe,Zeff,SOC,Alpha,Beta,Material\n");
                
                int cathode_pixel_count = 0;
                double alpha_sum = 0.0;
                double beta_sum = 0.0;
                
                // Count pixels with different lithiation degrees
                int lithiation_count[12] = {0}; // Corresponding to 11 lithiation materials
                double lithiation_alpha[12] = {0};
                double lithiation_beta[12] = {0};
                
                for (j = cathode_min_y; j <= cathode_max_y; j++) {
                    for (i = cathode_min_x; i <= cathode_max_x; i++) {
                        int idx = j * width + i;
                        double current_pe = Pe[idx];
                        double current_zeff = Zeff[idx];
                        double current_soc = soc_image[idx];
                        int material_idx = material_map[idx];
                        
                        // Only analyze pixels identified as LixMnO2 series
                        if (material_idx >= 2 && material_idx < NUM_MATERIALS) {
                            double alpha = materials[material_idx].alpha;
                            double beta = materials[material_idx].beta;
                            
                            fprintf(fp, "%d,%d,%.6e,%.6f,%.2f,%.2f,%.2f,%s\n", 
                                i, j, current_pe, current_zeff, current_soc, 
                                alpha, beta, materials[material_idx].name);
                            
                            alpha_sum += alpha;
                            beta_sum += beta;
                            cathode_pixel_count++;
                            
                            // Count lithiation distribution
                            int lith_index = material_idx - 2;
                            lithiation_count[lith_index]++;
                            lithiation_alpha[lith_index] = alpha;
                            lithiation_beta[lith_index] = beta;
                        }
                    }
                }
                
                fclose(fp);
                printf("Created lithiation analysis file: %s\n", lithiation_filename);
                
                if (cathode_pixel_count > 0) {
                    double avg_alpha = alpha_sum / cathode_pixel_count;
                    double avg_beta = beta_sum / cathode_pixel_count;
                    
                    printf("Cathode region lithiation statistics:\n");
                    printf("  Average α (Li content): %.4f\n", avg_alpha);
                    printf("  Average β (MnO2 content): %.4f\n", avg_beta);
                    printf("  Analyzed pixels: %d\n", cathode_pixel_count);
                    printf("  Approximate formula: Li_{%.2f}MnO_{%.2f}\n", avg_alpha, 2.0 * avg_beta);
                    
                    // Output lithiation distribution
                    printf("\nLithiation distribution:\n");
                    for (int m = 0; m < 11; m++) {
                        if (lithiation_count[m] > 0) {
                            double percentage = (double)lithiation_count[m] / cathode_pixel_count * 100.0;
                            printf("  Li%.1fMnO2 (α=%.1f, β=%.1f): %.1f%% (%d pixels)\n", 
                                   lithiation_alpha[m], lithiation_alpha[m], lithiation_beta[m], 
                                   percentage, lithiation_count[m]);
                        }
                    }
                }
            }
            
            // Create lithiation scatter plot
            char scatter_filename[256];
            sprintf_s(scatter_filename, sizeof(scatter_filename), "result\\LithiationScatter_%s.csv", state);
            error = fopen_s(&fp, scatter_filename, "w");
            if (error == 0) {
                fprintf(fp, "Alpha,Beta,SOC,Material\n");
                
                for (j = cathode_min_y; j <= cathode_max_y; j++) {
                    for (i = cathode_min_x; i <= cathode_max_x; i++) {
                        int idx = j * width + i;
                        double current_pe = Pe[idx];
                        double current_zeff = Zeff[idx];
                        double current_soc = soc_image[idx];
                        int material_idx = material_map[idx];
                        
                        if (material_idx >= 2 && material_idx < NUM_MATERIALS) {
                            double alpha = materials[material_idx].alpha;
                            double beta = materials[material_idx].beta;
                            
                            fprintf(fp, "%.2f,%.2f,%.2f,%s\n", alpha, beta, current_soc, materials[material_idx].name);
                        }
                    }
                }
                
                fclose(fp);
                printf("Created lithiation scatter plot file: %s\n", scatter_filename);
                printf("  Scatter plot format: Alpha(Li content), Beta(MnO2 content), SOC, Material\n");
            }
        }
        
        free(material_map);
        free(material_count);
    }

    // Output other result files
    // For single file, output raw data
    if (file_count == 1) {
        // SingleEnergy.csv
        error = fopen_s(&fp, ".\\result\\SingleEnergy.csv", "w");
        if (error == 0) {
            fprintf(fp, "DataFormat,2\nmemo\nX (pixel), Y (pixel), Z\n");
            for (j = 0; j < height; j++) {
                for (i = 0; i < width; i++) {
                    fprintf(fp, "%d,%d,%.6f\n", i, j, z[0][j * width + i]);
                }
            }
            fclose(fp);
            printf("Created: result\\SingleEnergy.csv\n");
        }
        else {
            printf("Failed to create SingleEnergy.csv\n");
        }

        // SingleEnergy.raw
        error = fopen_s(&fp, ".\\result\\SingleEnergy.raw", "wb");
        if (error == 0) {
            fwrite(z[0], sizeof(double), total_pixels, fp);
            fclose(fp);
            printf("Created: result\\SingleEnergy.raw\n");
        }
    }
    // For multiple files, output Pe, Zeff, etc.
    else {
        // Pe.csv
        error = fopen_s(&fp, ".\\result\\Pe.csv", "w");
        if (error == 0) {
            fprintf(fp, "DataFormat,2\nmemo\nX (pixel), Y (pixel), Z\n");
            for (j = 0; j < height; j++) {
                for (i = 0; i < width; i++) {
                    fprintf(fp, "%d,%d,%.6f\n", i, j, Pe[j * width + i]);
                }
            }
            fclose(fp);
            printf("Created: result\\Pe.csv\n");
        }
        else {
            printf("Failed to create Pe.csv\n");
        }

        // Zeff.csv
        error = fopen_s(&fp, ".\\result\\Zeff.csv", "w");
        if (error == 0) {
            fprintf(fp, "DataFormat,2\nmemo\nX (pixel), Y (pixel), Z\n");
            for (j = 0; j < height; j++) {
                for (i = 0; i < width; i++) {
                    fprintf(fp, "%d,%d,%.6f\n", i, j, Zeff[j * width + i]);
                }
            }
            fclose(fp);
            printf("Created: result\\Zeff.csv\n");
        }
    }

    // Statistical calculations
    printf("Calculating statistical results...\n");
    double sumValue = 0, sumZeff = 0, sumPe = 0, sumPere = 0;
    double* sumImg = (double*)calloc(file_count, sizeof(double));
    double aveValue = 0, aveZeff = 0, avePe = 0, avePere = 0;
    double s2_Value = 0, s2_Zeff = 0, s2_Pe = 0, s2_Pere = 0;
    double s_Value = 0, s_Zeff = 0, s_Pe = 0, s_Pere = 0;
    double Nin = 0;

    for (j = 0; j < height; j++) {
        for (i = 0; i < width; i++) {
            int in = func(i, j, ox, oy, r);
            if (in == 1) {
                Nin += 1;

                if (file_count == 1) {
                    sumValue += z[0][j * width + i];
                }
                else {
                    sumZeff += Zeff[j * width + i];
                    sumPe += Pe[j * width + i];
                    if (Pe_re != NULL) {
                        sumPere += Pe_re[j * width + i];
         
