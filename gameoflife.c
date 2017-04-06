#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "endian.h"
#include <mpi.h>

#define calcIndex(width, x, y)  ((y)*(width) + (x))

void show(unsigned *currentfield, int w, int h) {
    printf("\033[H");
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) printf(currentfield[calcIndex(w, x, y)] ? "\033[07m  \033[m" : "  ");
        printf("\033[E");
    }
    fflush(stdout);
}


float convert2BigEndian(const float inFloat) {
    float retVal;
    char *floatToConvert = (char *) &inFloat;
    char *returnFloat = (char *) &retVal;

    // swap the bytes into a temporary buffer
    returnFloat[0] = floatToConvert[3];
    returnFloat[1] = floatToConvert[2];
    returnFloat[2] = floatToConvert[1];
    returnFloat[3] = floatToConvert[0];

    return retVal;
}

void writeVTK(unsigned *currentfield, int w, int h, int t, char *prefix) {
    char name[1024] = "\0";
    sprintf(name, "%s_%d.vtk", prefix, t);
    FILE *outfile = fopen(name, "w");

    /*Write vtk header */
    fprintf(outfile, "# vtk DataFile Version 3.0\n");
    fprintf(outfile, "frame %d\n", t);
    fprintf(outfile, "BINARY\n");
    fprintf(outfile, "DATASET STRUCTURED_POINTS\n");
    fprintf(outfile, "DIMENSIONS %d %d %d \n", w, h, 1);
    fprintf(outfile, "SPACING 1.0 1.0 1.0\n");//or ASPECT_RATIO
    fprintf(outfile, "ORIGIN 0 0 0\n");
    fprintf(outfile, "POINT_DATA %d\n", h * w);
    fprintf(outfile, "SCALARS data float 1\n");
    fprintf(outfile, "LOOKUP_TABLE default\n");

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            float value = currentfield[calcIndex(w, x, y)]; // != 0.0 ? 1.0:0.0;
            value = convert2BigEndian(value);
            fwrite(&value, 1, sizeof(float), outfile);
        }
    }
    fclose(outfile);
}

int coutLifingsPeriodic(unsigned *currentfield, int x, int y, int w, int h) {
    int n = 0;
    for (int y1 = y - 1; y1 <= y + 1; y1++) {
        for (int x1 = x - 1; x1 <= x + 1; x1++) {
            if (currentfield[calcIndex(w, (x1 + w) % w, (y1 + h) % h)]) {
                n++;
            }
        }
    }
    return n;
}

void evolve(int *currentfield, int *new_part_field, int w, int h, int xStart, int xEnd, int yStart, int yEnd) {
    int changes = 0;
    for (int y = yStart; y < yEnd; y++) {
        for (int x = xStart; x < xEnd; x++) {
            //new_part_field[]
            int n = coutLifingsPeriodic(currentfield, x, y, w, h);
            if (currentfield[calcIndex(w, x, y)]) n--;
            new_part_field[calcIndex(w, x, y)] = (n == 3 || (n == 2 && currentfield[calcIndex(w, x, y)]));

            if (new_part_field[calcIndex(w, x, y)] - currentfield[calcIndex(w, x, y)] != 0) {
                changes++;
            }
        }
    }
}

void filling(unsigned *currentfield, int w, int h) {
    // for (int i = 0; i < h*w; i++) {
    //   currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    // }

    currentfield[calcIndex(w, 0, 0)] = 1;

    currentfield[calcIndex(w, 10, 4)] = 1;
    currentfield[calcIndex(w, 10, 5)] = 1;
    currentfield[calcIndex(w, 10, 6)] = 1;

    currentfield[calcIndex(w, 15, 4)] = 1;
    currentfield[calcIndex(w, 15, 5)] = 1;
    currentfield[calcIndex(w, 15, 6)] = 1;

    currentfield[calcIndex(w, 28, 4)] = 1;
    currentfield[calcIndex(w, 28, 5)] = 1;
    currentfield[calcIndex(w, 28, 6)] = 1;

    currentfield[calcIndex(w, 1, 28)] = 1;
    currentfield[calcIndex(w, 2, 28)] = 1;
    currentfield[calcIndex(w, 3, 28)] = 1;

    currentfield[calcIndex(w, 3, 15)] = 1;
    currentfield[calcIndex(w, 4, 15)] = 1;
    currentfield[calcIndex(w, 5, 15)] = 1;

}

void debug_print(int *field, int w, int h) {
    for (int y = 1; y < 15; y++) {
        for (int x = 1; x < 15; x++) {
            printf("%d", field[calcIndex(w, x, y)]);
        }
        printf("\n");
    }
}

void game(int w, int h, int timesteps) {
    unsigned *currentfield = calloc(w * h, sizeof(unsigned));
    unsigned *newfield = calloc(w * h, sizeof(unsigned));

    filling(currentfield, w, h);
    int rank, size;
    MPI_Status statusOben;
    MPI_Status statusUnten;
    MPI_Status statusLinks;
    MPI_Comm card_comm;
    int dim[1], periodic[1];
    periodic[0] = 0;
    int reorder = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // ProzessID
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    dim[0] = size;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dim, periodic, reorder, &card_comm);

    int *current_part_field = calloc(16 * 16, sizeof(double));
    int *new_part_field = calloc(14 * 14, sizeof(double));

    for (int t = 0; t < timesteps; t++) {
        int xStart, xEnd, yStart, yEnd;
        int obenGhost[14];
        int untenGhost[14];
        int linksGhost[14];
        int rechtsGhost[14];

        if (rank == 0) {
            xStart = 1;
            xEnd = 15;
            yStart = 1;
            yEnd = 15;
        } else if (rank == 1) {
            xStart = 15;
            xEnd = 29;
            yStart = 1;
            yEnd = 15;
        } else if (rank == 2) {
            xStart = 1;
            xEnd = 15;
            yStart = 15;
            yEnd = 29;
        } else if (rank == 3) {
            xStart = 15;
            xEnd = 29;
            yStart = 15;
            yEnd = 29;
        }

        // Mapping auf Teilfeld vom großem
        for (int x = 0; x < 14; x++) {
            for (int y = 0; y < 14; y++) {
                current_part_field[calcIndex(w, x + 1, y + 1)] = currentfield[calcIndex(w, x + xStart, y + yStart)];
            }
        }
        int sendenRechtsGhost[14];
        int sendenLinksGhost[14];
        int sendenObenGhost[14];
        int sendenUntenGhost[14];

        MPI_Barrier(MPI_COMM_WORLD);

        for (int y = 1; y < 15; y++) {
            sendenLinksGhost[y - 1] = current_part_field[calcIndex(w, 1, y)];
            sendenRechtsGhost[y - 1] = current_part_field[calcIndex(w, 14, y)];
        }

        for (int x = 1; x < 15; x++) {
            sendenObenGhost[x - 1] = current_part_field[calcIndex(w, x, 1)];
            sendenUntenGhost[x - 1] = current_part_field[calcIndex(w, x, 14)];
        }
        
        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0) {
            printf("rank rec:%d\n", rank);
            MPI_Recv(&obenGhost, 14, MPI_INT, 2, 96, MPI_COMM_WORLD, &statusOben);
            MPI_Recv(&rechtsGhost, 14, MPI_INT, 1, 97, MPI_COMM_WORLD, &statusLinks);
            MPI_Recv(&untenGhost, 14, MPI_INT, 2, 98, MPI_COMM_WORLD, &statusUnten);
            MPI_Recv(&linksGhost, 14, MPI_INT, 1, 99, MPI_COMM_WORLD, &statusLinks);

            printf("rank send:%d\n", rank);
            MPI_Send(&sendenUntenGhost, 14, MPI_INT, 2, 96, MPI_COMM_WORLD);
            MPI_Send(&sendenRechtsGhost, 14, MPI_INT, 1, 97, MPI_COMM_WORLD);
            MPI_Send(&sendenObenGhost, 14, MPI_INT, 2, 98, MPI_COMM_WORLD);
            MPI_Send(&sendenLinksGhost, 14, MPI_INT, 1, 99, MPI_COMM_WORLD);

            printf("rank ready:%d\n", rank);

        } else if (rank == 1) {
            printf("rank send:%d\n", rank);
            MPI_Send(&sendenUntenGhost, 14, MPI_INT, 3, 96, MPI_COMM_WORLD);
            MPI_Send(&sendenRechtsGhost, 14, MPI_INT, 0, 97, MPI_COMM_WORLD);
            MPI_Send(&sendenObenGhost, 14, MPI_INT, 3, 98, MPI_COMM_WORLD);
            MPI_Send(&sendenLinksGhost, 14, MPI_INT, 0, 99, MPI_COMM_WORLD);

            printf("rank rec:%d\n", rank);
            MPI_Recv(&obenGhost, 14, MPI_INT, 3, 96, MPI_COMM_WORLD, &statusOben);
            MPI_Recv(&rechtsGhost, 14, MPI_INT, 0, 97, MPI_COMM_WORLD, &statusLinks);
            MPI_Recv(&untenGhost, 14, MPI_INT, 3, 98, MPI_COMM_WORLD, &statusUnten);
            MPI_Recv(&linksGhost, 14, MPI_INT, 0, 99, MPI_COMM_WORLD, &statusLinks);
            printf("rank ready:%d\n", rank);
        } else if (rank == 2) {
            printf("rank send :%d\n", rank);
            MPI_Send(&sendenUntenGhost, 14, MPI_INT, 0, 96, MPI_COMM_WORLD);
            MPI_Send(&sendenRechtsGhost, 14, MPI_INT, 3, 97, MPI_COMM_WORLD);
            MPI_Send(&sendenObenGhost, 14, MPI_INT, 0, 98, MPI_COMM_WORLD);
            MPI_Send(&sendenLinksGhost, 14, MPI_INT, 3, 99, MPI_COMM_WORLD);

            printf("rank rec:%d\n", rank);
            MPI_Recv(&obenGhost, 14, MPI_INT, 0, 96, MPI_COMM_WORLD, &statusOben);
            MPI_Recv(&rechtsGhost, 14, MPI_INT, 3, 97, MPI_COMM_WORLD, &statusLinks);
            MPI_Recv(&untenGhost, 14, MPI_INT, 0, 98, MPI_COMM_WORLD, &statusUnten);
            MPI_Recv(&linksGhost, 14, MPI_INT, 3, 99, MPI_COMM_WORLD, &statusLinks);
            printf("rank ready:%d\n", rank);
        } else if (rank == 3) {
            printf("rank rec:%d\n", rank);
            MPI_Recv(&obenGhost, 14, MPI_INT, 1, 96, MPI_COMM_WORLD, &statusOben);
            MPI_Recv(&rechtsGhost, 14, MPI_INT, 2, 97, MPI_COMM_WORLD, &statusLinks);
            MPI_Recv(&untenGhost, 14, MPI_INT, 1, 98, MPI_COMM_WORLD, &statusUnten);
            MPI_Recv(&linksGhost, 14, MPI_INT, 2, 99, MPI_COMM_WORLD, &statusLinks);

            printf("rank send:%d\n", rank);
            MPI_Send(&sendenUntenGhost, 14, MPI_INT, 1, 96, MPI_COMM_WORLD);
            MPI_Send(&sendenRechtsGhost, 14, MPI_INT, 2, 97, MPI_COMM_WORLD);
            MPI_Send(&sendenObenGhost, 14, MPI_INT, 1, 98, MPI_COMM_WORLD);
            MPI_Send(&sendenLinksGhost, 14, MPI_INT, 2, 99, MPI_COMM_WORLD);
            printf("rank ready:%d\n", rank);

        }


        MPI_Barrier(MPI_COMM_WORLD);


        if (rank == 0) {
            for (int i = 0; i < 14; i++) {
                printf("%d %d %d %d %d\n", i, obenGhost[i], rechtsGhost[i], untenGhost[i], linksGhost[i]);
                current_part_field[calcIndex(w, xStart - 1, i + 1)] = linksGhost[i];
                current_part_field[calcIndex(w, xEnd, i + 1)] = rechtsGhost[i];
                current_part_field[calcIndex(w, i + 1, xStart - 1)] = obenGhost[i];
                current_part_field[calcIndex(w, i + 1, xEnd)] = untenGhost[i];
            }
            debug_print(current_part_field, w, h);
            evolve(current_part_field, new_part_field, w, h, 1, 15, 1, 15);
            exit(1);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        /*else if (rank == 1) {
            int sendenRechtsGhost[13];
            int sendenLinksGhost[13];
            for (int x = 15; x < 29; x++) {
                for (int y = 1; y < 15; y++) {
                    if (x == 15) {
                        sendenLinksGhost[y - 1] = currentfield[calcIndex(w, x, y)];
                    }
                    if (x == 28) {
                        sendenRechtsGhost[y - 1] = currentfield[calcIndex(w, x, y)];
                    }
                    current_part_field[calcIndex(w, x + 1, y)] = currentfield[calcIndex(w, x, y)];
                }
            }
            //MPI_Send(&sendenLinksGhost, 14, MPI_INT, 0, 99, MPI_COMM_WORLD);
            //MPI_Send(&sendenRechtsGhost, 14, MPI_INT, 0, 97, MPI_COMM_WORLD);
        }
        if (rank == 2) {
            // Von Unten nach Oben
            int sendenUntenGhost[13];

            // Von Oben nach Unten
            int sendenObenGhost[13];

            for (int x = 1; x < 15; x++) {
                for (int y = 15; y < 29; y++) {
                    current_part_field[calcIndex(w, x, y + 1)] = currentfield[calcIndex(w, x, y)];
                    if (y == 28) {
                        //printf("%d, %d, %d\n", x, y, currentfield[calcIndex(w, x, y)]);
                        sendenUntenGhost[x - 1] = currentfield[calcIndex(w, x, y)];
                    }
                    if (y == 15) {
                        sendenObenGhost[x - 1] = currentfield[calcIndex(w, x, y)];
                    }
                    current_part_field[calcIndex(w, x + 1, y + 1)] = currentfield[calcIndex(w, x, y)];
                }
            }

            MPI_Send(&sendenObenGhost, 14, MPI_INT, 0, 98, MPI_COMM_WORLD);
            MPI_Send(&sendenUntenGhost, 14, MPI_INT, 0, 96, MPI_COMM_WORLD);
        }

*/

/*
        if (rank == 0) {
            for (int i = 0; i < 14; i++) {
                printf("%d, %d \n", linksGhost[i], rechtsGhost[i]);
                //current_part_field[calcIndex(w, xStart-1, i+1)] = linksGhost[i];
                //current_part_field[calcIndex(w, xEnd, i+1)] = rechtsGhost[i];
                current_part_field[calcIndex(w, i+1, xStart-1)] = obenGhost[i];
                current_part_field[calcIndex(w, i+1, xEnd)] = untenGhost[i];
            }
            printf("\n Nachaustausch\n");
            debug_print(current_part_field, 0, 16, 0, 16, w, h);

            //MPI_Recv(&untenGhost, 16, MPI_INT, 2, 99, MPI_COMM_WORLD, &statusUnten);
            //MPI_Recv(&untenGhost, 16, MPI_INT, 2, 99, MPI_COMM_WORLD, &statusUnten);

            evolve(current_part_field, new_part_field, w, h, 1, 15, 1, 15);
            printf("\n Nach evolven");
            printf("\n");
&&debug_print(new_part_field, 0, 16, 0, 16, w, h);
        }



        if (rank == 3) {
            int oben2Ghost[13];
            for (int x = 15; x < 29; x++) {
                for (int y = 15; y < 29; y++) {
                    current_part_field[calcIndex(w, x + 1, y + 1)] = currentfield[calcIndex(w, x, y)];
                }
            }
            //MPI_Send(&oben2Ghost, 14, MPI_INT, 0, 99, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        exit(1);
        int field[w * h];
        MPI_Gather(&new_part_field, 4, MPI_INT, &field, 14 * 14, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            // recv buffer
            //debug_print(current_part_field, 1, 15, 1, 15, w, h);
            printf("\n");
            printf("\n");
            //debug_print(field, 1, 29, 1, 29, w, h);
            exit(1);
            show(field, w, h);
            //exit(1);
        }*/
    }

    free(currentfield);
    free(newfield);
}

int main(int c, char **v) {
    int w = 30, h = 30, timesteps = 50;
    if (c > 1) w = atoi(v[1]); ///< read width
    if (c > 2) h = atoi(v[2]); ///< read height
    if (c > 3) timesteps = atoi(v[3]);
    MPI_Init(&c, &v);

    // Numma f/ Prozess
    // 1D Gebietszerlegung b)
    game(w, h, timesteps);
    MPI_Finalize();

    //MPI_Finalize();
    return 0;
}
