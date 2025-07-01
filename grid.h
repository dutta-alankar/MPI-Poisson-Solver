#ifndef GRID_H
#define GRID_H

#include "vector.h"
#include <mpi.h>

// Represents a single cell in the grid
typedef struct {
    Vector3D pos; // Position of the cell center
    double density; // Density at the cell center
    double potential; // Potential at the cell center
} Cell;

// Represents the 3D grid, distributed across MPI processes
typedef struct {
    int global_size; // Global number of cells in each dimension
    double box_size; // Physical size of the simulation box
    Cell *cells; // Array of local cells (including ghost cells)

    int local_size_x, local_size_y, local_size_z; // Dimensions of the local grid (without ghost cells)
    int offset_x, offset_y, offset_z; // Offset of the local grid in the global grid

    int nghosts; // Number of ghost cells on each side

    int rank, size; // MPI rank and size
    MPI_Comm comm; // MPI communicator
} Grid;

// Function prototypes
Grid* create_distributed_grid(int global_size, double box_size, int nghosts, MPI_Comm comm);
void free_grid(Grid *grid);
void initialize_density(Grid *grid);
void exchange_ghost_cells(Grid *grid);

// Helper to get 1D index from 3D local coordinates (including ghosts)
static inline int get_idx(const Grid *grid, int i, int j, int k) {
    int ly = grid->local_size_y;
    int lz = grid->local_size_z;
    return i * ly * lz + j * lz + k;
}

#endif // GRID_H
