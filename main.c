#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stddef.h> // For offsetof
#include "grid.h"
#include "tree.h"
#include "gravity.h"

#define GRID_SIZE 64
#define BOX_SIZE 2.0
#define OPENING_ANGLE 0.5
#define NGHOSTS 2

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Create a custom MPI datatype for the Cell struct
    MPI_Datatype MPI_Cell;
    int blocklengths[3] = {3, 1, 1};
    MPI_Aint displacements[3];
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    displacements[0] = offsetof(Cell, pos);
    displacements[1] = offsetof(Cell, density);
    displacements[2] = offsetof(Cell, potential);
    MPI_Type_create_struct(3, blocklengths, displacements, types, &MPI_Cell);
    MPI_Type_commit(&MPI_Cell);

    Grid *grid = create_distributed_grid(GRID_SIZE, BOX_SIZE, NGHOSTS, comm);
    initialize_density(grid);

    exchange_ghost_cells(grid);

    TreeNode *tree = build_tree(grid);

    int start_index = 0; // Start at the beginning of the local grid
    int end_index = grid->local_size_x * grid->local_size_y * grid->local_size_z; // End at the end of the local grid
    calculate_potential(grid, tree, OPENING_ANGLE, start_index, end_index);

    // Gather results
    if (rank == 0) {
        int *displs = (int *)malloc(size * sizeof(int));
        int *recvcounts = (int *)malloc(size * sizeof(int));
        int total_cells = GRID_SIZE * GRID_SIZE * GRID_SIZE;
        Cell *global_cells = (Cell *)malloc(total_cells * sizeof(Cell));

        for (int i = 0; i < size; ++i) {
            int local_size_x = GRID_SIZE / size + (i < GRID_SIZE % size ? 1 : 0);
            int offset_x = i * (GRID_SIZE / size) + (i < GRID_SIZE % size ? i : GRID_SIZE % size);
            recvcounts[i] = local_size_x * GRID_SIZE * GRID_SIZE;
            displs[i] = offset_x * GRID_SIZE * GRID_SIZE;
        }

        MPI_Gatherv(grid->cells + NGHOSTS * grid->local_size_y * grid->local_size_z, 
                    grid->local_size_x * grid->local_size_y * grid->local_size_z, 
                    MPI_Cell, 
                    global_cells, recvcounts, displs, MPI_Cell, 0, comm);

        // Output results from the global grid
        FILE *potential_file = fopen("potential.dat", "w");
        FILE *density_file = fopen("density.dat", "w");

        for (int i = 0; i < GRID_SIZE; ++i) {
            for (int j = 0; j < GRID_SIZE; ++j) {
                int index = i * GRID_SIZE * GRID_SIZE + j * GRID_SIZE + GRID_SIZE / 2;
                fprintf(potential_file, "%f %f %f\n", global_cells[index].pos.x, global_cells[index].pos.y, global_cells[index].potential);
                fprintf(density_file, "%f %f %f\n", global_cells[index].pos.x, global_cells[index].pos.y, global_cells[index].density);
            }
            fprintf(potential_file, "\n");
            fprintf(density_file, "\n");
        }

        fclose(potential_file);
        fclose(density_file);
        free(global_cells);
        free(displs);
        free(recvcounts);

        printf("Simulation finished. Output written to potential.dat and density.dat\n");
    } else {
        MPI_Gatherv(grid->cells + NGHOSTS * grid->local_size_y * grid->local_size_z, 
                    grid->local_size_x * grid->local_size_y * grid->local_size_z, 
                    MPI_Cell, 
                    NULL, NULL, NULL, MPI_Cell, 0, comm);
    }

    free_grid(grid);
    free_tree(tree);

    MPI_Type_free(&MPI_Cell);
    MPI_Finalize();
    return 0;
}
