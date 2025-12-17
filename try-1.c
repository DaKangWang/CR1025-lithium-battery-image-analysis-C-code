#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <direct.h>
#include <io.h>

// Maximum number of files definitions
#define MAX_FILES 64
#define MAX_SELECTED_FILES 20

// Circle inside check function
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
    double all_energy[MAX_FILES] = { 0 };  // Store energy value for each file
    int total_files = 0;  // Total number of files

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
    printf("Circle center: (%.2f, %.2f), Radius: %.2f\n", ox, oy, r);

    // Read all files list, F, G values and energy values
    while (fscanf_s(cp, "%255s", filename, (unsigned)_countof(filename)) != EOF) {
        if (fscanf_s(cp, "%lf", &all_F[total_files]) != 1) break;
        if (fscanf_s(cp, "%lf", &all_G[total_files]) != 1) break;
        if (fscanf_s(cp, "%lf", &all_energy[total_files]) != 1) break;  // Read energy value

        strcpy_s(all_file_list[total_files], sizeof(all_file_list[total_files]), filename);
        printf("File %d: %s, F=%.2E, G=%.2E, Energy=%.2f keV\n", 
               total_files + 1, all_file_list[total_files], 
               all_F[total_files], all_G[total_files], all_energy[total_files]);
        total_files++;

        if (total_files >= MAX_FILES) break;
    }
    fclose(cp);
    printf("Total %d files loaded\n", total_files);

    // Let user select number of files and specific files
    int num_files = 0;
    printf("Select number of files to process (0=All, 1-20): ");
    scanf_s("%d", &num_files);

    // Validate input
    if (num_files == 0) {
        // Select all files
        num_files = total_files;
        if (num_files > MAX_SELECTED_FILES) {
            printf("Note: Total files(%d) exceed maximum selectable(%d), will only select first %d files.\n", total_files, MAX_SELECTED_FILES, MAX_SELECTED_FILES);
            num_files = MAX_SELECTED_FILES;
        }
        printf("Selected all %d files\n", num_files);
    }
    else if (num_files < 0 || num_files > MAX_SELECTED_FILES) {
        printf("Error: Invalid number of files. Please enter 0 (all) or integer between 1 and %d.\n", MAX_SELECTED_FILES);
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

    // Determine image size based on file size
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
        region_max_x = 68;
        region_min_y = 22;
        region_max_y = 85;
        printf("Fresh state region: x[%d-%d], y[%d-%d]\n", 
               region_min_x, region_max_x, region_min_y, region_max_y);
    }
    else { // discharged
        // Discharged state region definition
        region_min_x = 45;
        region_max_x = 69;
        region_min_y = 31;
        region_max_y = 82;
        printf("Discharged state region: x[%d-%d], y[%d-%d]\n", 
               region_min_x, region_max_x, region_min_y, region_max_y);
    }

    // Store regional average linear attenuation coefficient for each file
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
            printf("File %d specified region average: %.4f\n", k + 1, region_avg[k]);
        }
        else {
            region_avg[k] = 0;
            printf("File %d specified region has no valid pixels\n", k + 1);
        }

        for (i = 0; i < total_pixels; i++) {
            sum_m[i] += a[i];
        }
    }

    // Create energy vs linear attenuation coefficient CSV file
    char mu_energy_filename[256];
    sprintf_s(mu_energy_filename, sizeof(mu_energy_filename), "result\\mu_vs_energy_%s.csv", state);
    error = fopen_s(&fp, mu_energy_filename, "w");
    if (error == 0) {
        fprintf(fp, "energy,mu_avg\n");
        for (k = 0; k < file_count; k++) {
            fprintf(fp, "%.6f,%.6f\n", selected_energy[k], region_avg[k]);
        }
        fclose(fp);
        printf("Created energy vs linear attenuation coefficient file: %s\n", mu_energy_filename);
    }
    else {
        printf("Failed to create energy vs linear attenuation coefficient file\n");
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

    // Calculate Pe, Zeff, Pe_re (only when more than 1 file)
    if (file_count > 1) {
        printf("Calculating Pe, Zeff and Pe_re...\n");

        // Pre-calculate sum_f and sum_g
        double sum_f = 0.0, sum_g = 0.0;
        for (k = 0; k < file_count; k++) {
            sum_f += f[k];
            sum_g += g[k];
        }

        for (i = 0; i < total_pixels; i++) {
            // Special formula for 2 files
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
            // Use least squares method for 3 or more files
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
    }

    // Output results (CSV format)
    printf("Outputting result files...\n");

    // Define regions based on battery state
    FILE* anode_file = NULL, * cathode_file = NULL, * separator_file = NULL;
    FILE* casing_left_file = NULL, * casing_right_file = NULL;

    // Define region boundaries
    int anode_min_x, anode_max_x, anode_min_y, anode_max_y;
    int cathode_min_x, cathode_max_x, cathode_min_y, cathode_max_y;
    int separator_min_x, separator_max_x, separator_min_y, separator_max_y;
    int casing_left_min_x, casing_left_max_x, casing_left_min_y, casing_left_max_y;
    int casing_right_min_x, casing_right_max_x, casing_right_min_y, casing_right_max_y;

    // Region counters
    int anode_count = 0, cathode_count = 0, separator_count = 0;
    int casing_left_count = 0, casing_right_count = 0;

    if (strcmp(state, "fresh") == 0) {
        // Fresh state region definition
        // Anode
        anode_min_x = 63; anode_max_x = 65;
        anode_min_y = 22; anode_max_y = 85;

        // Cathode
        cathode_min_x = 47; cathode_max_x = 61;
        cathode_min_y = 22; cathode_max_y = 85;

        // Separator - modified to horizontal axis 61-62
        separator_min_x = 61; separator_max_x = 62;
        separator_min_y = 22; separator_max_y = 85;

        // Left casing
        casing_left_min_x = 44; casing_left_max_x = 46;
        casing_left_min_y = 4; casing_left_max_y = 103;

        // Right casing
        casing_right_min_x = 66; casing_right_max_x = 68;
        casing_right_min_y = 20; casing_right_max_y = 91;

        printf("Fresh state region definitions:\n");
        printf("  Anode: x[%d-%d], y[%d-%d]\n", anode_min_x, anode_max_x, anode_min_y, anode_max_y);
        printf("  Cathode: x[%d-%d], y[%d-%d]\n", cathode_min_x, cathode_max_x, cathode_min_y, cathode_max_y);
        printf("  Separator: x[%d-%d], y[%d-%d]\n", separator_min_x, separator_max_x, separator_min_y, separator_max_y);
        printf("  Left casing: x[%d-%d], y[%d-%d]\n", casing_left_min_x, casing_left_max_x, casing_left_min_y, casing_left_max_y);
        printf("  Right casing: x[%d-%d], y[%d-%d]\n", casing_right_min_x, casing_right_max_x, casing_right_min_y, casing_right_max_y);

        // Create anode file
        error = fopen_s(&anode_file, "result/Anode_Li_metal.csv", "w");
        if (error == 0) {
            // Decide header based on file count
            if (file_count == 1) {
                fprintf(anode_file, "x,y,Value,Region\n");
            }
            else {
                fprintf(anode_file, "x,y,Pe,Zeff,Region\n");
            }
        }
        else {
            printf("Failed to create Anode_Li_metal.csv\n");
        }
    }
    else { // discharged
        // Discharged state region definition
        // Anode disappeared
        anode_min_x = -1; anode_max_x = -1;
        anode_min_y = -1; anode_max_y = -1;

        // Cathode
        cathode_min_x = 49; cathode_max_x = 60;
        cathode_min_y = 31; cathode_max_y = 82;

        // Separator
        separator_min_x = 61; separator_max_x = 66;
        separator_min_y = 30; separator_max_y = 83;

        // Left casing
        casing_left_min_x = 45; casing_left_max_x = 47;
        casing_left_min_y = 9; casing_left_max_y = 101;

        // Right casing
        casing_right_min_x = 67; casing_right_max_x = 69;
        casing_right_min_y = 28; casing_right_max_y = 86;

        printf("Discharged state region definitions:\n");
        printf("  Anode: disappeared\n");
        printf("  Cathode: x[%d-%d], y[%d-%d]\n", cathode_min_x, cathode_max_x, cathode_min_y, cathode_max_y);
        printf("  Separator: x[%d-%d], y[%d-%d]\n", separator_min_x, separator_max_x, separator_min_y, separator_max_y);
        printf("  Left casing: x[%d-%d], y[%d-%d]\n", casing_left_min_x, casing_left_max_x, casing_left_min_y, casing_left_max_y);
        printf("  Right casing: x[%d-%d], y[%d-%d]\n", casing_right_min_x, casing_right_max_x, casing_right_min_y, casing_right_max_y);
    }

    // Create common region files
    error = fopen_s(&cathode_file, "result/Cathode_MnO2.csv", "w");
    if (error == 0) {
        // Decide header based on file count
        if (file_count == 1) {
            fprintf(cathode_file, "x,y,Value,Region\n");
        }
        else {
            fprintf(cathode_file, "x,y,Pe,Zeff,Region\n");
        }
    }
    else {
        printf("Failed to create Cathode_MnO2.csv\n");
    }

    error = fopen_s(&separator_file, "result/Separator.csv", "w");
    if (error == 0) {
        // Decide header based on file count
        if (file_count == 1) {
            fprintf(separator_file, "x,y,Value,Region\n");
        }
        else {
            fprintf(separator_file, "x,y,Pe,Zeff,Region\n");
        }
    }
    else {
        printf("Failed to create Separator.csv\n");
    }

    error = fopen_s(&casing_left_file, "result/SUS_casing_left.csv", "w");
    if (error == 0) {
        // Decide header based on file count
        if (file_count == 1) {
            fprintf(casing_left_file, "x,y,Value,Region\n");
        }
        else {
            fprintf(casing_left_file, "x,y,Pe,Zeff,Region\n");
        }
    }
    else {
        printf("Failed to create SUS_casing_left.csv\n");
    }

    error = fopen_s(&casing_right_file, "result/SUS_casing_right.csv", "w");
    if (error == 0) {
        // Decide header based on file count
        if (file_count == 1) {
            fprintf(casing_right_file, "x,y,Value,Region\n");
        }
        else {
            fprintf(casing_right_file, "x,y,Pe,Zeff,Region\n");
        }
    }
    else {
        printf("Failed to create SUS_casing_right.csv\n");
    }

    for (j = 0; j < height; j++) {
        for (i = 0; i < width; i++) {
            int idx = j * width + i;

            // Anode region (fresh state only)
            if (strcmp(state, "fresh") == 0 &&
                i >= anode_min_x && i <= anode_max_x &&
                j >= anode_min_y && j <= anode_max_y) {
                if (anode_file) {
                    if (file_count == 1) {
                        fprintf(anode_file, "%d,%d,%.6f,Anode\n", i, j, z[0][idx]);
                    }
                    else {
                        fprintf(anode_file, "%d,%d,%.6f,%.6f,Anode\n", i, j, Pe[idx], Zeff[idx]);
                    }
                    anode_count++;
                }
            }
            // Cathode region
            if (i >= cathode_min_x && i <= cathode_max_x &&
                j >= cathode_min_y && j <= cathode_max_y) {
                if (cathode_file) {
                    if (file_count == 1) {
                        fprintf(cathode_file, "%d,%d,%.6f,Cathode\n", i, j, z[0][idx]);
                    }
                    else {
                        fprintf(cathode_file, "%d,%d,%.6f,%.6f,Cathode\n", i, j, Pe[idx], Zeff[idx]);
                    }
                    cathode_count++;
                }
            }
            // Separator region
            if (i >= separator_min_x && i <= separator_max_x &&
                j >= separator_min_y && j <= separator_max_y) {
                if (separator_file) {
                    if (file_count == 1) {
                        fprintf(separator_file, "%d,%d,%.6f,Separator\n", i, j, z[0][idx]);
                    }
                    else {
                        fprintf(separator_file, "%d,%d,%.6f,%.6f,Separator\n", i, j, Pe[idx], Zeff[idx]);
                    }
                    separator_count++;
                }
            }
            // Left casing region
            if (i >= casing_left_min_x && i <= casing_left_max_x &&
                j >= casing_left_min_y && j <= casing_left_max_y) {
                if (casing_left_file) {
                    if (file_count == 1) {
                        fprintf(casing_left_file, "%d,%d,%.6f,Casing_Left\n", i, j, z[0][idx]);
                    }
                    else {
                        fprintf(casing_left_file, "%d,%d,%.6f,%.6f,Casing_Left\n", i, j, Pe[idx], Zeff[idx]);
                    }
                    casing_left_count++;
                }
            }
            // Right casing region
            if (i >= casing_right_min_x && i <= casing_right_max_x &&
                j >= casing_right_min_y && j <= casing_right_max_y) {
                if (casing_right_file) {
                    if (file_count == 1) {
                        fprintf(casing_right_file, "%d,%d,%.6f,Casing_Right\n", i, j, z[0][idx]);
                    }
                    else {
                        fprintf(casing_right_file, "%d,%d,%.6f,%.6f,Casing_Right\n", i, j, Pe[idx], Zeff[idx]);
                    }
                    casing_right_count++;
                }
            }
        }
    }

    // Close region files
    if (anode_file) fclose(anode_file);
    if (cathode_file) fclose(cathode_file);
    if (separator_file) fclose(separator_file);
    if (casing_left_file) fclose(casing_left_file);
    if (casing_right_file) fclose(casing_right_file);

    // Output region statistics
    printf("\nRegion data statistics:\n");
    printf("  Anode region points: %d\n", anode_count);
    printf("  Cathode region points: %d\n", cathode_count);
    printf("  Separator region points: %d\n", separator_count);
    printf("  Left casing region points: %d\n", casing_left_count);
    printf("  Right casing region points: %d\n", casing_right_count);
    printf("Region data files created\n");

    // Output other result files
    // Output raw data when single file
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
    // Output Pe, Zeff etc. when multiple files
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

        // Pe_re.csv
        error = fopen_s(&fp, ".\\result\\Pe_re.csv", "w");
        if (error == 0) {
            fprintf(fp, "DataFormat,2\nmemo\nX (pixel), Y (pixel), Z\n");
            for (j = 0; j < height; j++) {
                for (i = 0; i < width; i++) {
                    fprintf(fp, "%d,%d,%.6f\n", i, j, Pe_re[j * width + i]);
                }
            }
            fclose(fp);
            printf("Created: result\\Pe_re.csv\n");
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

        // Pe.raw
        error = fopen_s(&fp, ".\\result\\Pe.raw", "wb");
        if (error == 0) {
            fwrite(Pe, sizeof(double), total_pixels, fp);
            fclose(fp);
            printf("Created: result\\Pe.raw\n");
        }

        // Zeff.raw
        error = fopen_s(&fp, ".\\result\\Zeff.raw", "wb");
        if (error == 0) {
            fwrite(Zeff, sizeof(double), total_pixels, fp);
            fclose(fp);
            printf("Created: result\\Zeff.raw\n");
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
                    sumPere += Pe_re[j * width + i];
                }

                for (k = 0; k < file_count; k++) {
                    sumImg[k] += z[k][j * width + i];
                }
            }
        }
    }

    // Calculate averages
    if (Nin > 0) {
        if (file_count == 1) {
            aveValue = sumValue / Nin;
        }
        else {
            aveZeff = sumZeff / Nin;
            avePe = sumPe / Nin;
            avePere = sumPere / Nin;
        }

        // Calculate standard deviation
        for (j = 0; j < height; j++) {
            for (i = 0; i < width; i++) {
                int in = func(i, j, ox, oy, r);
                if (in == 1) {
                    double diff;

                    if (file_count == 1) {
                        diff = z[0][j * width + i] - aveValue;
                        s2_Value += diff * diff;
                    }
                    else {
                        diff = Zeff[j * width + i] - aveZeff;
                        s2_Zeff += diff * diff;

                        diff = Pe[j * width + i] - avePe;
                        s2_Pe += diff * diff;

                        diff = Pe_re[j * width + i] - avePere;
                        s2_Pere += diff * diff;
                    }
                }
            }
        }

        if (file_count == 1) {
            s_Value = sqrt(s2_Value / Nin);
        }
        else {
            s_Zeff = sqrt(s2_Zeff / Nin);
            s_Pe = sqrt(s2_Pe / Nin);
            s_Pere = sqrt(s2_Pere / Nin);
        }
    }

    // Output statistical results
    error = fopen_s(&fp, ".\\result\\result.csv", "w");
    if (error == 0) {
        fprintf(fp, "Number of files: %d\n", file_count);
        fprintf(fp, "Processing state: %s\n", state);

        if (file_count == 1) {
            fprintf(fp, "Parameter,Average,StdDev,SNR\n");
            fprintf(fp, "RawData,%.6f,%.6f,%.6f\n", aveValue, s_Value, aveValue / s_Value);
        }
        else {
            fprintf(fp, "Parameter,Average,StdDev,SNR\n");
            fprintf(fp, "Pe,%.6f,%.6f,%.6f\n", avePe, s_Pe, avePe / s_Pe);
            fprintf(fp, "Zeff,%.6f,%.6f,%.6f\n", aveZeff, s_Zeff, aveZeff / s_Zeff);
            fprintf(fp, "Pe_re,%.6f,%.6f,%.6f\n", avePere, s_Pere, avePere / s_Pere);
        }

        // Add region information
        fprintf(fp, "\nRegion statistics:\n");
        fprintf(fp, "Region,Points\n");
        fprintf(fp, "Anode,%d\n", anode_count);
        fprintf(fp, "Cathode,%d\n", cathode_count);
        fprintf(fp, "Separator,%d\n", separator_count);
        fprintf(fp, "LeftCasing,%d\n", casing_left_count);
        fprintf(fp, "RightCasing,%d\n", casing_right_count);

        fclose(fp);
        printf("Created: result\\result.csv\n");
    }

    // Free memory
    free(a);
    free(sum_m);
    free(z_final_average);
    free(sumImg);

    if (file_count > 1) {
        free(Pe);
        free(Zeff);
        free(Pe_re);
    }

    for (k = 0; k < file_count; k++) {
        free(z[k]);
    }
    free(z);

    printf("Processing completed!\n");
    return 0;
}
