// Released under MIT License.
// Copyright (c) 2022-2023 Ladislav Bartos

#include <unistd.h>
#include <groan.h>

const char VERSION[] = "v2023/01/11";

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;

// axis in which the cross-section shall be calculated
typedef enum axis {none = 0, xaxis = 1, yaxis = 2, zaxis = 3} axis_t;

/*
 * Get user-defined grid dimension.
 */
int get_range(const char *optarg, float array[2])
{
    if (sscanf(optarg, "%f-%f", &array[0], &array[1]) != 2 && 
        sscanf(optarg, "%f - %f", &array[0], &array[1]) != 2 &&
        sscanf(optarg, "%f %f", &array[0], &array[1]) != 2) {
        fprintf(stderr, "Could not understand grid dimension specifier.\n");
        return 1;
    }

    return 0;
}

/*
 * Parses command line arguments.
 * Returns zero, if parsing has been successful. Else returns non-zero.
 */
int get_arguments(
        int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **output_file,
        char **atoms,
        char **radiifile,
        axis_t *axis,
        float *array_dimx,
        float *array_dimy,
        float *array_dimz,
        int *grid_density)
{
    int gro_specified = 0, atoms_specified = 0;

    int opt = 0;
    while((opt = getopt(argc, argv, "c:f:n:o:s:r:i:x:y:z:d:h")) != -1) {
        switch (opt) {
        // help
        case 'h':
            return 1;
        // gro file to read
        case 'c':
            *gro_file = optarg;
            gro_specified = 1;
            break;
        // xtc file to read
        case 'f':
            *xtc_file = optarg;
            break;
        // ndx file to read
        case 'n':
            *ndx_file = optarg;
            break;
        // output file
        case 'o':
            *output_file = optarg;
            break;
        // atoms to work with
        case 's':
            *atoms = optarg;
            atoms_specified = 1;
            break;
        case 'r':
            *radiifile = optarg;
            break;
        // axis along which the cross sectional area should be calculated
        case 'i':
            if (!strcmp(optarg, "x")) *axis = xaxis;
            else if (!strcmp(optarg, "y")) *axis = yaxis;
            else if (!strcmp(optarg, "z")) *axis = zaxis;
            else {
                fprintf(stderr, "Unknown axis %s\n", optarg);
                *axis = none;
                return 1;
            }
            break;
        // grid dimensions
        case 'x':
            if (get_range(optarg, array_dimx) != 0) return 1;
            break;
        case 'y':
            if (get_range(optarg, array_dimy) != 0) return 1;
            break;
        case 'z':
            if (get_range(optarg, array_dimz) != 0) return 1;
            break;
        // grid density
        case 'd':
            sscanf(optarg, "%d", grid_density);
            if (*grid_density <= 0) {
                fprintf(stderr, "Grid density must be > 0.\n");
                return 1;
            }
            break;
        default:
            //fprintf(stderr, "Unknown command line option: %c.\n", opt);
            return 1;
        }
    }

    if (!gro_specified || !atoms_specified) {
        fprintf(stderr, "Gro file and atoms specification must always be supplied.\n");
        return 1;
    }
    return 0;
}

void print_usage(const char *program_name)
{
    printf("Usage: %s -c GRO_FILE -s SELECTION [OPTION]...\n", program_name);
    printf("\nOPTIONS\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read (optional)\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-o STRING        output file (default: area.xvg)\n");
    printf("-s STRING        selection of atoms which cross-section shall be calculated\n");
    printf("-r STRING        file containing van der Waals radii of relevant atoms (default: vdwradii.dat)\n");
    printf("-i CHAR          axis in which the cross-section shall be calculated (default: z)\n");
    printf("-x FLOAT-FLOAT   grid dimension in x axis (default: box size from gro file)\n");
    printf("-y FLOAT-FLOAT   grid dimension in y axis (default: box size from gro file)\n");
    printf("-z FLOAT-FLOAT   grid dimension in z axis (default: box size from gro file)\n");
    printf("-d INTEGER       density of the grid used for calculation in points per nm (default: 10)\n");
    printf("\n");
}

/*
 * Prints parameters that the program will use for the calculation.
 */
void print_arguments(
        FILE *stream,
        const char *gro_file, 
        const char *xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *atoms,
        const char *radii_file,
        const axis_t axis,
        const float array_dimx[2],
        const float array_dimy[2],
        const float array_dimz[2],
        const int grid_density)
{
    fprintf(stream, "Parameters for Cross-Sectional Area calculation:\n");
    fprintf(stream, ">>> gro file:         %s\n", gro_file);
    if (xtc_file == NULL) fprintf(stream, ">>> xtc file:         ----\n");
    else fprintf(stream, ">>> xtc file:         %s\n", xtc_file);
    fprintf(stream, ">>> ndx file:         %s\n", ndx_file);
    fprintf(stream, ">>> output file:      %s\n", output_file);
    fprintf(stream, ">>> atoms:            %s\n", atoms);
    fprintf(stream, ">>> vdwradii file:    %s\n", radii_file);
    fprintf(stream, ">>> grid dimensions:  x: %.1f - %.1f nm, y: %.1f - %.1f nm, z: %.1f - %.1f nm\n", 
            array_dimx[0], array_dimx[1], 
            array_dimy[0], array_dimy[1], 
            array_dimz[0], array_dimz[1]);
    if (axis == xaxis) fprintf(stream, ">>> cross-axis:       %c\n", 'x');
    else if (axis == yaxis) fprintf(stream, ">>> cross-axis:       %c\n", 'y');
    else if (axis == zaxis) fprintf(stream, ">>> cross-axis:       %c\n", 'z');
    fprintf(stream, ">>> grid density:     %d points per nm\n\n", grid_density);
}

/*
 * Reads file containing van der Waals radii of atoms.
 */
dict_t *read_radii(const char *radiifile)
{
    FILE *file = fopen(radiifile, "r");
    if (file == NULL) return NULL;

    char line[1024] = "";
    dict_t *vdwradii = dict_create();

    while (fgets(line, 1024, file) != NULL) {
        char resname[512] = "";
        char atomname[512] = "";
        float radius = 0.0;

        if (sscanf(line, "%s %s %f\n", resname, atomname, &radius) != 3) {
            continue;
        }
        
        char identifier[1024] = "";
        sprintf(identifier, "%s-%s", resname, atomname);

        dict_set(vdwradii, identifier, &radius, sizeof(float));
    }

    fclose(file);
    return vdwradii;
}

/* 
 * Converts index of an array to coordinate.
 */
static inline float index2coor(int x, float minx, int grid_density)
{
    return (float) x / grid_density + minx;
}

/*
 * Converts coordinate to an index array.
 */
static inline size_t coor2index(float x, float minx, int grid_density)
{
    return (size_t) roundf((x - minx) * grid_density);
}

static inline void wrap_grid_size(long *x, long min_index, long max_index)
{
    if (*x < min_index) *x = min_index;
    if (*x > max_index) *x = max_index;
}

/*
 * Sets grid points to 1 if they are inside the vdw radius of any atom.
 * Returns zero if successful, non-zero in case of an error.
 */
int fill_grid(
        const system_t *system,
        const atom_selection_t *selection, 
        const dict_t *vdwradii,
        short *grid, 
        const size_t xsize,
        const size_t ysize,
        const size_t zsize,
        const float array_dimx[2],
        const float array_dimy[2],
        const float array_dimz[2],
        const int grid_density)
{
    // loop through all atoms of the selection
    for (size_t i = 0; i < selection->n_atoms; ++i) {

        atom_t *atom = selection->atoms[i];

        char id[1024] = "";
        sprintf(id, "%s-%s", atom->residue_name, atom->atom_name);
        float *radius = (float *) dict_get(vdwradii, id);
        if (radius == NULL) {
            fprintf(stderr, "Could not find vdw radius for atom %s of residue %s.\n", atom->atom_name, atom->residue_name);
            return 1;
        }

        // loop through nearby grid points
        long min_grid_x = coor2index(atom->position[0] - *radius, array_dimx[0], grid_density);
        long max_grid_x = coor2index(atom->position[0] + *radius, array_dimx[0], grid_density);
        long min_grid_y = coor2index(atom->position[1] - *radius, array_dimy[0], grid_density);
        long max_grid_y = coor2index(atom->position[1] + *radius, array_dimy[0], grid_density);
        long min_grid_z = coor2index(atom->position[2] - *radius, array_dimz[0], grid_density);
        long max_grid_z = coor2index(atom->position[2] + *radius, array_dimz[0], grid_density);

        wrap_grid_size(&min_grid_x, 0, (long) xsize - 1);
        wrap_grid_size(&max_grid_x, 0, (long) xsize - 1);
        wrap_grid_size(&min_grid_y, 0, (long) ysize - 1);
        wrap_grid_size(&max_grid_y, 0, (long) ysize - 1);
        wrap_grid_size(&min_grid_z, 0, (long) zsize - 1);
        wrap_grid_size(&max_grid_z, 0, (long) zsize - 1);

        //printf("%ld %ld %ld %ld %ld %ld\n", min_grid_x, max_grid_x, min_grid_y, max_grid_y, min_grid_z, max_grid_z);

        for (size_t x = (size_t) min_grid_x; x <= (size_t) max_grid_x; ++x) {
            for (size_t y = (size_t) min_grid_y; y <= (size_t) max_grid_y; ++y) {
                for (size_t z = (size_t) min_grid_z; z <= (size_t) max_grid_z; ++z) {
                    if (x >= xsize || y >= ysize || z >= zsize) continue;

                    vec_t real_coor = {index2coor(x, array_dimx[0], grid_density),
                                       index2coor(y, array_dimy[0], grid_density),
                                       index2coor(z, array_dimz[0], grid_density)};
                    if (distance3D(real_coor, atom->position, system->box) < *radius) {
                        grid[z * xsize * ysize + y * xsize + x] = 1;
                    }
                }
            }
        }
    } 
    return 0;
}

/*
 * Calculates area of a slice.
 */
float cross_section(
        const short *grid, 
        const size_t refindex,
        const axis_t axis,
        const size_t xsize, 
        const size_t ysize,
        const size_t zsize,
        const float array_dimx[2],
        const float array_dimy[2],
        const float array_dimz[2],
        const int grid_density)
{
    // calculate number of points inside the specified atoms in a slice
    size_t points_inside = 0;
    size_t total_points = 0;
    float total_area = 0.0;
    float grid_density_r = (float) 1 / grid_density;
    if (axis == xaxis) {
        for (size_t z = 0; z < zsize; ++z) {
            for (size_t y = 0; y < ysize; ++y) {
                points_inside += grid[z * xsize * ysize + y * xsize + refindex];
            }
        }
        total_points = zsize * ysize;
        total_area = (array_dimz[1] - array_dimz[0] + grid_density_r) * (array_dimy[1] - array_dimy[0] + grid_density_r);

    } else if (axis == yaxis) {
        for (size_t z = 0; z < zsize; ++z) {
            for (size_t x = 0; x < xsize; ++x) {
                points_inside += grid[z * xsize * ysize + refindex * xsize + x];
            }
        }
        total_points = zsize * xsize;
        total_area = (array_dimz[1] - array_dimz[0] + grid_density_r) * (array_dimx[1] - array_dimx[0] + grid_density_r);


    } else if (axis == zaxis) {
        for (size_t y = 0; y < ysize; ++y) {
            for (size_t x = 0; x < xsize; ++x) {
                points_inside += grid[refindex * xsize * ysize + y * xsize + x];
            }
        }
        total_points = xsize * ysize;
        total_area = (array_dimy[1] - array_dimy[0] + grid_density_r) * (array_dimx[1] - array_dimx[0] + grid_density_r);
    }

    // calculate the cross-section
    return ((float) (points_inside) / total_points) * total_area;
}

int main(int argc, char **argv)
{
    // get arguments
    char *gro_file = NULL;
    char *xtc_file = NULL;
    char *ndx_file = "index.ndx";
    char *output_file = "area.xvg";
    char *atoms = NULL;
    char *radii_file = "vdwradii.dat";
    axis_t axis = zaxis;
    float array_dimx[2] = {0.};
    float array_dimy[2] = {0.};
    float array_dimz[2] = {0.};
    int grid_density = 10;

    if (get_arguments(argc, argv, &gro_file, &xtc_file, &ndx_file, &output_file, &atoms, &radii_file, &axis, array_dimx, array_dimy, array_dimz, &grid_density) != 0) {
        print_usage(argv[0]);
        return 1;
    }

    printf("\n");

    // read van der waals file
    dict_t *vdwradii = read_radii(radii_file);
    if (vdwradii == NULL) {
        fprintf(stderr, "Could not read %s as radii file.\n", radii_file);
        return 1;
    }

    // read gro file
    system_t *system = load_gro(gro_file);
    if (system == NULL) {
        dict_destroy(vdwradii);
        return 1;
    }

    // if array dimensions were not set, get them from gro file
    if (array_dimx[0] == 0 && array_dimx[1] == 0) {
        array_dimx[1] = system->box[0];  
    }
    if (array_dimy[0] == 0 && array_dimy[1] == 0) {
        array_dimy[1] = system->box[1];
    }
    if (array_dimz[0] == 0 && array_dimz[1] == 0) {
        array_dimz[1] = system->box[2];
    }

    // check that the array dimensions don't have nonsensical values
    if (array_dimx[0] >= array_dimx[1] || array_dimy[0] >= array_dimy[1] || array_dimz[0] >= array_dimz[1]) {
        fprintf(stderr, "Nonsensical array dimensions.\n");
        dict_destroy(vdwradii);
        free(system);
        return 1;
    }

    print_arguments(stdout, gro_file, xtc_file, ndx_file, output_file, atoms, radii_file, axis, array_dimx, array_dimy, array_dimz, grid_density);

    // try opening output file
    FILE *output = fopen(output_file, "w");
    if (output == NULL) {
        fprintf(stderr, "Output file could not be opened.\n");
        dict_destroy(vdwradii);
        free(system);
        return 1;
    }

    // try reading ndx file (ignore if this fails)
    dict_t *ndx_groups = read_ndx(ndx_file, system);

    // select all atoms
    atom_selection_t *all = select_system(system);

    // select specified atoms
    atom_selection_t *selected = smart_select(all, atoms, ndx_groups);
    free(all);
    all = NULL;

    // check that the selection has been successful
    if (selected == NULL || selected->n_atoms == 0) {
        fprintf(stderr, "No atoms ('%s') found.\n", atoms);

        dict_destroy(ndx_groups);
        dict_destroy(vdwradii);
        free(system);
        free(selected);
        fclose(output);
        return 1;
    }

    // prepare 3D grid
    size_t xsize = (size_t) roundf( (array_dimx[1] - array_dimx[0]) * grid_density ) + 1;
    size_t ysize = (size_t) roundf( (array_dimy[1] - array_dimy[0]) * grid_density ) + 1;
    size_t zsize = (size_t) roundf( (array_dimz[1] - array_dimz[0]) * grid_density ) + 1;
    size_t n_points = xsize * ysize * zsize;
    short *grid = calloc(n_points, sizeof(short));

    if (grid == NULL) {
        fprintf(stderr, "Could not allocate memory (grid too large?)\n");
        dict_destroy(ndx_groups);
        dict_destroy(vdwradii);
        free(system);
        free(selected);
        fclose(output);
        return 1;
    }

    // write header for the output file
    fprintf(output, "# Generated with crosection (C Cross-Sectional Area Calculator) %s\n", VERSION);
    fprintf(output, "# Command line: ");
    for (int i = 0; i < argc; ++i) {
        fprintf(output, "%s ", argv[i]);
    }
    fprintf(output, "\n");

    fprintf(output, "@  title \"Cross-sectional area\"\n");
    if (axis == xaxis) fprintf(output, "@  xaxis label \"x coordinate [nm]\"\n");
    else if (axis == yaxis) fprintf(output, "@  xaxis label \"y coordinate [nm]\"\n");
    else if (axis == zaxis) fprintf(output, "@  xaxis label \"z coordinate [nm]\"\n");
    
    fprintf(output, "@  yaxis label \"cross-sectional area [nm^2]\"\n");
    fprintf(output, "@  s1 legend \"%s\"\n", atoms);

    // variables for the cross-sectional area calculation
    size_t dim = 0;
    float *array = NULL;
    if (axis == xaxis) {
        dim = xsize;
        array = array_dimx;
    }
    else if (axis == yaxis) {
        dim = ysize;
        array = array_dimy;
    }
    else if (axis == zaxis) {
        dim = zsize;
        array = array_dimz;
    }

    // if there is no xtc file supplied, analyze only the gro file
    if (xtc_file == NULL) {
        if (fill_grid(system, selected, vdwradii, grid, xsize, ysize, zsize, array_dimx, array_dimy, array_dimz, grid_density) != 0) {
            dict_destroy(ndx_groups);
            dict_destroy(vdwradii);
            free(system);
            free(selected);
            fclose(output);
            free(grid);
            return 1;
        }
        
        for (size_t i = 0; i < dim; ++i) {
            fprintf(output, "%f %f\n", index2coor(i, array[0], grid_density), cross_section(grid, i, axis, xsize, ysize, zsize, array_dimx, array_dimy, array_dimz, grid_density));
        }
    } else {
        XDRFILE *xtc = xdrfile_open(xtc_file, "r");
        if (xtc == NULL) {
            fprintf(stderr, "File %s could not be read as an xtc file.\n", xtc_file);
            dict_destroy(ndx_groups);
            dict_destroy(vdwradii);
            free(system);
            free(selected);
            fclose(output);
            free(grid);
            return 1;
        }

        // check that the gro file and the xtc file match each other
        if (!validate_xtc(xtc_file, (int) system->n_atoms)) {
            fprintf(stderr, "Number of atoms in %s does not match %s.\n", xtc_file, gro_file);
            xdrfile_close(xtc);
            dict_destroy(ndx_groups);
            dict_destroy(vdwradii);
            free(system);
            free(selected);
            fclose(output);
            free(grid);
            return 1;
        }

        float *area = calloc(dim, sizeof(float));
        size_t n_frames = 0;

        while (read_xtc_step(xtc, system) == 0) {
            // print info about the progress of reading
            if ((int) system->time % PROGRESS_FREQ == 0) {
                printf("Step: %d. Time: %.0f ps\r", system->step, system->time);
                fflush(stdout);
            }

            // fill grid
            if (fill_grid(system, selected, vdwradii, grid, xsize, ysize, zsize, array_dimx, array_dimy, array_dimz, grid_density) != 0) {
                dict_destroy(ndx_groups);
                dict_destroy(vdwradii);
                free(system);
                free(selected);
                xdrfile_close(xtc);
                fclose(output);
                free(grid);
                free(area);
                return 1;
            }


            // calculate cross-section
            for (size_t i = 0; i < dim; ++i) {
                area[i] += cross_section(grid, i, axis, xsize, ysize, zsize, array_dimx, array_dimy, array_dimz, grid_density);
            }

            n_frames++;   
            memset(grid, 0, xsize * ysize * zsize * sizeof(short));
        }

        xdrfile_close(xtc);
        for (size_t i = 0; i < dim; ++i) {
            fprintf(output, "%f %f\n", index2coor(i, array[0], grid_density), area[i] / n_frames);
        }

        free(area);
        printf("\n");

    }

    dict_destroy(vdwradii);
    dict_destroy(ndx_groups);
    free(selected);
    free(system);
    free(grid);
    fclose(output);

    return 0;
}