#include "GalCos_variables.h"
#include "hdf5.h"
#include<malloc.h>

//This function would take care of partitioning the unclustered particles (the oclustered ones are the ones in the lines of the files
//between 0 and Npart_clustered, then the unclustered would be the ones between Npart_clustered and Npart_Total)
//this function should take two integers (a,b) to cut the unclustered particles in a subgroup of unclustered particle between a and b

int GalCos_PartitionParticles(int a, int b)
{
    a += Npart_clustered;
    b += Npart_clustered;

    printf("\nPartitioning from %d to %d\n", a, b); fflush(stdout);

    Particle = (struct part *) malloc((size_t) Npart_unclustered * sizeof(struct part));

    int start_idx[NumFiles], end_idx[NumFiles], offset = 0;

    // Precalculate global indices for each file
    for (int Fileid = 0; Fileid < NumFiles; Fileid++) {
        char infile[100];
        sprintf(infile, "%s.%d.hdf5", data_prefix, Fileid);

        hid_t file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) continue;

        hid_t header_group = H5Gopen(file_id, "/Header", H5P_DEFAULT);
        if (header_group < 0) {
            H5Fclose(file_id);
            continue;
        }

        hid_t attr_id = H5Aopen(header_group, "NumPart_ThisFile", H5P_DEFAULT);
        int npart[6] = {0};
        H5Aread(attr_id, H5T_STD_U32LE, npart);
        H5Aclose(attr_id);
        H5Gclose(header_group);
        H5Fclose(file_id);

        start_idx[Fileid] = offset;
        end_idx[Fileid] = offset + npart[1];
        offset = end_idx[Fileid];
    }

    // Loop files and read particles
    int local_idx = 0;
    for (int Fileid = 0; Fileid < NumFiles; Fileid++) {
        // If the file does not contain particles in the range [a, b), skip it
        if (end_idx[Fileid] <= a-Npart_clustered || start_idx[Fileid] >= b-Npart_clustered)
            continue;

        char infile[100];
        sprintf(infile, "%s.%d.hdf5", data_prefix, Fileid);
        hid_t file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) continue;

        hid_t group_id = H5Gopen(file_id, "/PartType1", H5P_DEFAULT);
        hid_t dataset_id = H5Dopen(group_id, "Coordinates", H5P_DEFAULT);
        hid_t filespace_id = H5Dget_space(dataset_id);
        hid_t memspace_id = H5Screate_simple(2, (hsize_t[]){1, 3}, NULL);

        for (int i = start_idx[Fileid]; i < end_idx[Fileid]; i++) {
            if (i < a-Npart_clustered || i >= b-Npart_clustered)
                continue;

            hsize_t offset_h[2] = {i - start_idx[Fileid], 0}; // Offset in the hyperslab, meaning the particle index in this file
            hsize_t count[2] = {1, 3};
            H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_h, NULL, count, NULL);

            float coord[3];
            herr_t status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace_id, filespace_id, H5P_DEFAULT, coord);

            if (status >= 0) {
                Particle[local_idx].pos[0] = coord[0];
                Particle[local_idx].pos[1] = coord[1];
                Particle[local_idx].pos[2] = coord[2];

                local_idx++;
            } else {
                fprintf(stderr, "Error leyendo partícula %d\n", i);
            }
        }

        H5Sclose(memspace_id);
        H5Sclose(filespace_id);
        H5Dclose(dataset_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
    }
    printf("Partition completed. \n");
    fflush(stdout);
    return 0;
}
