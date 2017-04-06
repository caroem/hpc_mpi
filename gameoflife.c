#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <endian.h>
#include <mpi.h>
 
#define calcIndex(width, x,y)  ((y)*(width) + (x))

void show(unsigned* currentfield, int w, int h) {
  printf("\033[H");
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) printf(currentfield[calcIndex(w, x,y)] ? "\033[07m  \033[m" : "  ");
    printf("\033[E");
  }
  fflush(stdout);
}
 

float convert2BigEndian( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat    = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

void writeVTK(unsigned* currentfield, int w, int h, int t, char* prefix) {
  char name[1024] = "\0";
  sprintf(name, "%s_%d.vtk", prefix, t);
  FILE* outfile = fopen(name, "w");

  /*Write vtk header */                                                           
  fprintf(outfile,"# vtk DataFile Version 3.0\n");       
  fprintf(outfile,"frame %d\n", t);     
  fprintf(outfile,"BINARY\n");     
  fprintf(outfile,"DATASET STRUCTURED_POINTS\n");     
  fprintf(outfile,"DIMENSIONS %d %d %d \n", w, h, 1);        
  fprintf(outfile,"SPACING 1.0 1.0 1.0\n");//or ASPECT_RATIO                            
  fprintf(outfile,"ORIGIN 0 0 0\n");                                              
  fprintf(outfile,"POINT_DATA %d\n", h*w);                                    
  fprintf(outfile,"SCALARS data float 1\n");                              
  fprintf(outfile,"LOOKUP_TABLE default\n");         
 
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      float value = currentfield[calcIndex(w, x,y)]; // != 0.0 ? 1.0:0.0;
      value = convert2BigEndian(value);
      fwrite(&value, 1, sizeof(float), outfile);
    }
  }
  fclose(outfile);
}

int coutLifingsPeriodic(unsigned* currentfield, int x , int y, int w, int h) {
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
 
unsigned* evolve(unsigned* currentfield, unsigned* array, int w, int h, int xStart, int xEnd, int yStart, int yEnd) {
int changes = 0;
int new_part_field[14*14];
  int * copy = malloc(sizeof(unsigned) * (14*14));
  mempcpy(copy, array, (14*14)*sizeof(unsigned));
  for (int y = yStart; y < yEnd; y++) {
    for (int x = xStart; x < xEnd; x++) {
      //new_part_field[]
      int n = coutLifingsPeriodic(currentfield, x , y, w, h);
      if (currentfield[calcIndex(w, x,y)]) n--;
      new_part_field[calcIndex(w, x,y)] = (n == 3 || (n == 2 && currentfield[calcIndex(w, x,y)]));
      
      if (new_part_field[calcIndex(w, x,y)] - currentfield[calcIndex(w, x,y)] != 0) {
        changes++;
      }
    }
  }
  //return new_part_field*;
  return copy;
}
 
void filling(unsigned* currentfield, int w, int h) {
  // for (int i = 0; i < h*w; i++) {
  //   currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
  // }
  currentfield[calcIndex(w, 10,4)] = 1;
  currentfield[calcIndex(w, 10,5)] = 1;
  currentfield[calcIndex(w, 10,6)] = 1;
  currentfield[calcIndex(w, 10,28)] = 1;

}
 
void game(int w, int h, int timesteps) {
  printf("Asdf");
  unsigned *currentfield = calloc(w*h, sizeof(unsigned));
  unsigned *newfield     = calloc(w*h, sizeof(unsigned));

  
  filling(currentfield, w, h);
  int rank, size;
  MPI_Status statusMPI;
  MPI_Status statusMPI2;
  MPI_Comm card_comm;
  int dim[1], periodic [1];
  periodic[0] = 0;
  int reorder = 1;
   
  MPI_Comm_size(MPI_COMM_WORLD, &size); // ProzessID
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  dim[0] = size;
  MPI_Cart_create(MPI_COMM_WORLD, 1, dim, periodic, reorder, &card_comm);


  int current_part_field[16*16];
  unsigned *part_field = calloc(14*14, sizeof(unsigned));
  //int part_field[14*14];
  for (int t = 0; t < timesteps; t++) {
  
    if(rank == 0) {
      for (int x = 1; x < 15; x++) {
        for(int y = 1; y< 15; y++) {
          current_part_field[calcIndex(w, x, y)] = currentfield[calcIndex(w,x,y)];
        }
      }

      // Erhaelt von 1 und 3 Ghost Layer, sendet 4
      int obenGhost[13];
      int untenGhost[13];
      MPI_Recv(&obenGhost, 16, MPI_INT, 2, 99, MPI_COMM_WORLD, &statusMPI);
      MPI_Recv(&untenGhost, 16, MPI_INT, 2, 99, MPI_COMM_WORLD, &statusMPI2);
      for(int test = 0; test < 14; test++) {
        //printf("%d\n", obenGhost[test]);
        //printf("%d\n", untenGhost[test]);
      }

      part_field = evolve(current_part_field, *part_field, w, h, 1, 15, 1, 15);

    }
    if(rank == 1) {
      for (int x = 15; x < 29; x++) {
        for(int y = 1; y < 15; y++) {
          current_part_field[calcIndex(w, x+1, y)] = currentfield[calcIndex(w,x,y)];
        }
      }
    }
    if(rank == 2) {
      int oben2Ghost[13];
      int unten2Ghost[13];
      for (int x = 1; x < 15; x++) {
        for(int y = 15; y < 29; y++) {
          current_part_field[calcIndex(w, x, y+1)] = currentfield[calcIndex(w,x,y)];
          if (y == 28) {
            oben2Ghost[x-1] = currentfield[calcIndex(w,x,y)];
            //printf("%d \n", x);
            //printf('%d', oben2Ghost[x-15]);
          }
          if (y == 15) {
            unten2Ghost[x-1] = currentfield[calcIndex(w,x,y)];
          }
          current_part_field[calcIndex(w, x+1, y+1)] = currentfield[calcIndex(w,x,y)];
        }
      }
      MPI_Send(&oben2Ghost, 14, MPI_INT, 0, 99, MPI_COMM_WORLD);
      MPI_Send(&unten2Ghost, 14, MPI_INT, 0, 99, MPI_COMM_WORLD);
    }

    if(rank == 3) {
      int oben2Ghost[13];
      for (int x = 15; x < 29; x++) {
        for(int y = 15; y < 29; y++) {
          current_part_field[calcIndex(w, x+1, y+1)] = currentfield[calcIndex(w,x,y)];
        }
      }
      //MPI_Send(&oben2Ghost, 14, MPI_INT, 0, 99, MPI_COMM_WORLD);
    }
    if(rank == 0) {
      // recv buffer
        int field[w*h];

    }
    
   //show(currentfield, w, h);
   exit(1);

   //int changes = 
   int xStart = 1;
   int xEnd = 15;
   int yStart = 1;
   int yEnd = 15;
   //evolve(currentfield, w, h, xStart, xEnd, yStart, yEnd);
   int allchanges;
   //MPI_Allreduce(&changes, &allchanges, size, MPI_INT, MPI_SUM, card_comm);
   //if (allchanges == 0) break;
   //if (changes == 0) break;
 
    usleep(200000);

    //SWAP
    unsigned *temp = currentfield;
    currentfield = newfield;
    newfield = temp;
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
