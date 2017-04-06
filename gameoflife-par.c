#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

#define calcIndex(width, x,y)  ((y)*(width) + (x))
long TimeSteps = 100;

void writeVTK2(long timestep, double *data, char prefix[1024], long w, long h) {
  int thread_num = omp_get_thread_num();
  int xStart = (thread_num % 2) * w / 2;
  int xEnd = (thread_num % 2 + 1) * w / 2;
  int yStart = (thread_num / 2) * h / 2;
  int yEnd = (thread_num / 2 + 1) * h / 2;
  if (xEnd > w) {
    xEnd = w;
  }
  if (yEnd > h) {
    yEnd = h;
  }
  
  char filename[2048];  
  int x,y; 
  
  long offsetX=0;
  long offsetY=0;
  float deltax=1.0;
  float deltay=1.0;
  long  nxy = w * h * sizeof(float);  

  snprintf(filename, sizeof(filename), "%s-%05ld%s", prefix, timestep, ".vti");
  FILE* fp = fopen(filename, "w");

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n", offsetX, offsetX + w-1, offsetY, offsetY + h-1, 0, 0, deltax, deltax, 0.0);
  fprintf(fp, "<CellData Scalars=\"%s\">\n", prefix);
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n", prefix);
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</ImageData>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fwrite((unsigned char*)&nxy, sizeof(long), 1, fp);

  // for (y = 0; y < h; y++) {
  //   for (x = 0; x < w; x++) {
  //     float value = data[calcIndex(h, x,y)];
  //     fwrite((unsigned char*)&value, sizeof(float), 1, fp);
  //   }
  // }
  // for (y = yStart; y < yEnd; y++) {
  //   for (x = xStart; x < xEnd; x++) {
  //     float value = data[calcIndex(h, x,y)];
  //     fwrite((unsigned char*)&value, sizeof(float), 1, fp);
  //   }
  // }
   for (y = yStart; y <= yEnd; y++) {
    for (x = xStart; x <= xEnd; x++) {
      float value = data[calcIndex(w,x,y)];
      fwrite((unsigned char*)&value, sizeof(float), 1, fp);
    }
}
  
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}


void show(double* currentfield, int w, int h) {
  printf("\033[H");
  int x,y;
  for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) printf(currentfield[calcIndex(w, x,y)] ? "\033[07m  \033[m" : "  ");
    //printf("\033[E");
    printf("\n");
  }
  fflush(stdout);
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

int evolve(unsigned* currentfield, unsigned* newfield, int w, int h) {
int changes = 0;
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      int n = coutLifingsPeriodic(currentfield, x , y, w, h);
      if (currentfield[calcIndex(w, x,y)]) n--;
      newfield[calcIndex(w, x,y)] = (n == 3 || (n == 2 && currentfield[calcIndex(w, x,y)]));

      if (newfield[calcIndex(w, x,y)] - currentfield[calcIndex(w, x,y)] != 0) {
        changes++;
      }
    }
  }
  return changes;
}

void filling(unsigned* currentfield, int w, int h) {
  for (int i = 0; i < h*w; i++) {
    currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
  }
}
 
void game2(int w, int h) {
  double *currentfield = calloc(w*h, sizeof(double));
  double *newfield     = calloc(w*h, sizeof(double));
  // reduction(+:changed)
  //printf("size unsigned %d, size long %d\n",sizeof(float), sizeof(long));
  filling(currentfield, w, h);
  long t;
  for (t=0;t<TimeSteps;t++) {
    bool changed;
    changed = false;
    show(currentfield, w, h);

    #pragma omp parallel num_threads(4) reduction(+:changed)
    {
      //writeVTK2(t,currentfield,"gol", xStart, xEnd, yStart, yEnd, w, h, omp_get_thread_num());
      //changed = evolve(currentfield, newfield, w, h);
    }
      if (changed == false) {
        sleep(3);
        break;
      }

    // evolve(currentfield, newfield, w, h);
    
    // printf("%ld timestep\n",t);
    // writeVTK2(t,currentfield,"gol", w, h);
    
    usleep(200000);

    //SWAP
    double *temp = currentfield;
    currentfield = newfield;
    newfield = temp;
  }
  free(currentfield);
  free(newfield);
}

void game(int w, int h) {
  double *currentfield = calloc(w*h, sizeof(double));
  double *newfield     = calloc(w*h, sizeof(double));
  int xStart, xEnd, yStart, yEnd;
  //printf("size unsigned %d, size long %d\n",sizeof(float), sizeof(long));
  
  filling(currentfield, w, h);
  long t;
  int changed;
  for (t=0;t<TimeSteps;t++) {
    show(currentfield, w, h);
    // reduction(+:changed)
    #pragma omp parallel num_threads(4) private(xStart, xEnd, yStart, yEnd) firstprivate(w,h)
    {
      int thread_num = omp_get_thread_num();

      int xStart = (thread_num % 2) * w;
      int xEnd = (((thread_num % 2) + 1) * w) - 1;
      int yStart = (thread_num / 2) * h;
      int yEnd = (((thread_num / 2) + 1) * h) -1;
      if (thread_num % 2 == 1) {
        xEnd = w - 1;
      }
      if (thread_num / 2 == 1) {
        yEnd = h - 1;
      }
      writeVTK2(t,currentfield,"gol", w, h);
      changed = evolve(currentfield, newfield, w, h, xStart, xEnd, yStart, yEnd);
      //#pragma omp barrier
    }
    
    printf("%ld timestep\n",t);
    //writeVTK2(t,currentfield,"gol", w, h);
    
    usleep(200000);

    //SWAP
    double *temp = currentfield;
    currentfield = newfield;
    newfield = temp;
  }
  
  free(currentfield);
  free(newfield);
  
}
 
int main(int c, char **v) {
  int w = 0, h = 0;
  if (c > 1) w = atoi(v[1]); ///< read width
  if (c > 2) h = atoi(v[2]); ///< read height
  if (w <= 0) w = 30; ///< default width
  if (h <= 0) h = 30; ///< default height
  game(w, h);
}
