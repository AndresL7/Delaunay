#include <hdf5.h>
#include <stdlib.h>
#include <stdio.h>

#define GROUP_NAME "/PartType1"
#define DATASET_NAME "Coordinates2"
#define SPACE_DIM 3

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Uso: %s <archivo.h5>\n", argv[0]);
        return 1;
    }

    const char *file_name = argv[1];
    hid_t file_id, group_id, dataset_id, dataspace_id;
    hsize_t dims[2];
    herr_t status;
    
    int n = 512;
    dims[0] = (hsize_t)n * n * n;
    dims[1] = SPACE_DIM;

    dataspace_id = H5Screate_simple(2, dims, NULL);

    file_id = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
        printf("Error al abrir el archivo: %s\n", file_name);
        return 1;
    }

    group_id = H5Gopen(file_id, GROUP_NAME, H5P_DEFAULT);
    if (group_id < 0) {
        printf("Error al abrir el grupo: %s\n", GROUP_NAME);
        H5Fclose(file_id);
        return 1;
    }

    if (H5Lexists(group_id, DATASET_NAME, H5P_DEFAULT) > 0) {
        status = H5Ldelete(group_id, DATASET_NAME, H5P_DEFAULT);
        if (status < 0) {
            printf("Error al eliminar el dataset existente: %s\n", DATASET_NAME);
            H5Gclose(group_id);
            H5Fclose(file_id);
            return 1;
        }
    }

    dataset_id = H5Dcreate(group_id, DATASET_NAME, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        printf("Error al crear el dataset: %s\n", DATASET_NAME);
        H5Sclose(dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return 1;
    }

    // Reservar memoria en el heap
    float (*data)[SPACE_DIM] = malloc(dims[0] * sizeof(*data));
    if (data == NULL) {
        printf("Error al reservar memoria\n");
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return 1;
    }

    int index = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                data[index][0] = i * (100.0 / (n - 1));
                data[index][1] = j * (100.0 / (n - 1));
                data[index][2] = k * (100.0 / (n - 1));
                index++;
            }
        }
    }

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0) {
        printf("Error al escribir los datos en el dataset\n");
    }

    // Liberar memoria
    free(data);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Gclose(group_id);
    H5Fclose(file_id);

    printf("Dataset creado y datos escritos exitosamente en el archivo: %s\n", file_name);
    return 0;
}
