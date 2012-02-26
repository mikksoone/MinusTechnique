// MinusTechnique.cpp : Defines the entry point for the console application.
//

// If Visual Studio (windows) then...else linux
// Could maybe use _WIN32 instead
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <climits>
// Neat profiler that currently isn't available for everyone
// Therefore it's commented out for the git
// #include "../ProfilingTimer/ProfilingTimer.h"
// #define DEBUG_TIMER

// Verbal version for debugging
// #define PRINT_DEBUG
#define PRINT_STATUS

#define SORT          // Could be undef if one want's to compare time differences
#define INT int       // Maybe a double is needed?

#define HARDCODED_DATA // If we don't want to specify input/output files 

// The timers in windows and linux are different
#ifdef _WIN32
#define TIMER_TYPE __int64
#define int64_t __int64
#else
#define TIMER_TYPE double
#endif


// Global data
int FILE_err=0;
INT nEdges=0, nNodes=0;
INT *g_conform=0;
double g_quadPart;
// Global time data
static TIMER_TYPE g_sort_time = 0;
static TIMER_TYPE g_find_time = 0;
static TIMER_TYPE g_timeWhenSorted = 0;
double g_timePerLine = LLONG_MAX;
bool g_bSort = 0;
static int loops_after_sort = 0;
static int elems_removed = 0;

typedef struct {
   INT *elem_buf;   // transactions, its buffer, and their sizes
   INT **elem;
   INT *elem_count;
   INT *frq;
   INT *conform;
   INT *seq; //sequence of rows eliminated
   INT seq_count;
   std::vector<INT> rows_left;
} TRSACT;


static inline void sort(TRSACT * T);

/* Helper function for debugging INT arrays 
   TODO: change %d to %f if needed */
void print_int_arr(INT *arr, INT len, const char *name){
   INT i;
   //return;
   for( i=0 ; i<len ; i++){
      printf("%d ", arr[i]);
   }
   printf ("==%s==========\n", name);
}

/* allocate memory with error routine (Takeaki Uno)*/
char *alloc_memory (size_t siz){
  char *p = (char *)calloc (siz, 1);
  if ( p == NULL ){
    printf ("out of memory\n");
    exit (1);
  }
  return (p);
}

/* re-allocate memory and fill the newly allocated part by 0, with error routine (Takeaki Uno)*/
char *realloc_memory (char *p, int sizof, size_t siz, size_t *org_siz){
  size_t i;
  p = (char *)realloc(p, siz*sizof);
  if ( p == NULL ){
    printf ("out of memory\n");
    exit (1);
  }
  for ( i=(*org_siz)*sizof ; i<siz*sizof ; i++ ) p[i] = 0;
  *org_siz = siz;
  return (p);
}

/* Compare routine for partial sort, currently not used */
int psort_cmp_conform (const INT x, const INT y){
  if ( g_conform[x] <= g_conform[y] ) return (-1);
  else return ( g_conform[x] > g_conform[y] );
}
/* Compare routine for qsort, reorders row items based on their conform */
int qsort_cmp_conform (const void *x, const void *y){
  if ( g_conform[*((INT *)x)] <= g_conform[*((INT *)y)] ) return (-1);
  else return ( g_conform[*((INT *)x)] > g_conform[*((INT *)y)] );
}

/* read an integer from the file (Takeaki Uno)*/
INT FILE_read_int (FILE *fp){
  INT item;
  int flag =1;
  int ch;
  FILE_err = 0;
  do {
    ch = fgetc (fp);

    if ( ch == '\n' ){ FILE_err = 5; return (0); }
    if ( ch < 0 ){ FILE_err = 6; return (0); }
    if ( ch=='-' ) flag = -1;
  } while ( ch<'0' || ch>'9' );
  for ( item=(int)(ch-'0') ; 1 ; item=item*10 +(int)(ch-'0') ){
    ch = fgetc (fp);
    if ( ch == '\n' ){ FILE_err = 1; return (flag*item); }
    if ( ch < 0 ){ FILE_err = 2; return (flag*item); }
    if ( (ch < '0') || (ch > '9')) return (flag*item);
  }
}

/* load a transaction from the input file to memory (Takeaki Uno style)*/
void TRSACT_file_load (TRSACT *T, const char *fname)
{
#ifdef DEBUG_TIMER   
   TIMER("load file");
#endif
   INT item, i, j;
   FILE *fp = fopen (fname,"r");

   if ( !fp ){ printf ("file open error\n"); exit (1); }

   nNodes = FILE_read_int (fp);
   if ( (FILE_err&4) != 0){ printf ("file structure error 1\n"); exit (1); }
   nEdges = FILE_read_int (fp);
   if ( (FILE_err&4) != 0){ printf ("file structure error 2\n"); exit (1); }
   
   T->elem_buf = (INT *)alloc_memory ( sizeof(INT) * (nEdges + nNodes) );
   T->elem = (INT **)alloc_memory( sizeof(INT) * nNodes );
   T->elem_count = (INT *)alloc_memory( sizeof(INT) * nNodes );
   i = 0; //counter
   int old_node = -1;
   do 
   {
      INT node = FILE_read_int (fp);
      INT item = FILE_read_int (fp);
      //FILE_read_int(fp); //get EOL
      if( old_node != node)
      {
         T->elem[node] = &T->elem_buf[i]; //save the starting point of a node
         T->elem_buf[i++] = node;
         ++T->elem_count[node];
      }
      T->elem_buf[i] = item;
      ++T->elem_count[node];
      i++;
      old_node = node;
  } while ( (FILE_err&2)==0);

   fclose(fp);
}

/* Prints the reordered table to a file 
void TRSACT_output_result(TRSACT * T, const char *fname)
{
#ifdef DEBUG_TIMER   
   TIMER("output file");
#endif
   INT i, j;
   FILE *fp = fopen (fname,"w");

   if ( !fp ){ printf ("file open error\n"); exit (1); }

   for( i=0 ; i<nCol ; i++ )
   {
      for( j=0 ; j<nRow ; j++) 
      {
         fprintf( fp, "%d ", T->buf[ T->seq[j]*nCol+i ] + MIN_ITEM );
      }
      fprintf(fp, "\n");
   }
   fclose(fp);
}
*/
/* Frees the transaction and global data */
void TRSACT_free( TRSACT * T )
{
#ifdef DEBUG_TIMER   
   TIMER("free");
#endif
   free( T->elem_buf );
   free( T->elem );
   free( T->elem_count );

   free( T->conform );
   free( T->seq );

   free( g_conform );

   free( T->frq );
}

/* Alloc the neccessary memory for transactions and do the initial calculations */
void TRSACT_init( TRSACT * T )
{
#ifdef DEBUG_TIMER   
   TIMER("init");
#endif
   INT i,j;

   T->frq = (INT*) alloc_memory ( sizeof(INT) * nNodes );
   
   // Fill the frequency table
   for( j=0 ; j<nNodes ; j++)
   {
      for( i=0 ; i<T->elem_count[j] ; i++ )
         ++T->frq[T->elem[j][i]];
   }


   T->conform      = (INT*) alloc_memory( sizeof(INT) * nNodes );
   T->seq         = (INT*) alloc_memory( sizeof(INT) * nNodes );
   T->seq_count   = 0;

   
   g_conform = (INT*) alloc_memory( sizeof(INT) * nNodes );

   // We add number from 0 to nRow to the rows_left vector
   // These number, of course, represent the rows that have not been kicked out just yet
   T->rows_left.reserve(nNodes);
   for( i=0 ; i<nNodes ; i++ )
   {
      if( T->elem_count[i] > 0 ) T->rows_left.push_back(i);
   }

   sort(T);
}

/* Here we transform the file buffer from vertical to horizontal 
void TRSACT_switch(TRSACT * T)
{
#ifdef DEBUG_TIMER   
   TIMER("switch");
#endif
   printf("Switching..\n");
   INT i, j, cnt = 0;
   // Fill the tmp container with values from ordered buf 
   INT * tmp = (INT*) alloc_memory( sizeof(INT) * nRow*nCol );
   for( i=0 ; i<nCol ; i++ )
      for( j=0 ; j<nRow ; j++)   
         tmp[cnt++] = T->buf[ T->seq[j]*nCol +i ]; 
   
#ifdef PRINT_DEBUG
   for( j=0; j<nRow ; j++ )
      print_int_arr(&T->buf[T->seq[j]*nCol], nCol, "buf ordered");
#endif
   // Free the currently used transaction container   
   TRSACT_free(T);

   T->buf = tmp;
   INT k = nRow;
   nRow = nCol;
   nCol = k;

#ifdef PRINT_DEBUG
   for( j=0 ; j<nRow ; j++ )
      print_int_arr(&T->buf[j*nCol], nCol, "buf switched");
#endif

   // Init again..
   TRSACT_init(T);
   printf("Switching done\n");
}
*/

/* Multiplatform function for getting the time in one big number */
static inline TIMER_TYPE get_time()
{
#ifdef _WIN32
   static LARGE_INTEGER li;
   QueryPerformanceCounter(&li);
   return li.QuadPart;
#else
   static timeval t;
   gettimeofday(&t, NULL);
   return (t.tv_sec * 1000.0) + t.tv_usec/1000.0;
#endif
}

/* Function for calculating fresh global conforms and then sorting the rows
   according to the just calculated conforms. 
   So, global conforms are the conforms at the time of the last sort 
*/
static inline void sort(TRSACT * T)
{
#ifdef DEBUG_TIMER   
   TIMER("sort");
#endif 
   // We measure the time it takes to calculate all coforms and sort the rows
   static TIMER_TYPE start_time = 0;
   start_time = get_time();
   g_timeWhenSorted = get_time();
   memset( &g_conform[0], 0, sizeof(INT) * nNodes );
   // Calcualte conforms
   for( auto it = T->rows_left.begin() ; it != T->rows_left.end() ; ++it )
   {
      for( int i=0 ; i<T->elem_count[*it] ; i++ )
         g_conform[*it] += T->frq[ T->elem[*it][i] ];
   }

   // And here we sort the rows according to their conform values
   std::qsort(&T->rows_left[0], T->rows_left.size(), sizeof(INT), qsort_cmp_conform);
   g_sort_time = get_time() - start_time;
   g_bSort = false;
   g_timePerLine = LLONG_MAX;
   loops_after_sort=0;
   elems_removed = 0;
   g_find_time = 0;
}

/* Routine for finding the row with minimum conform */
static inline bool find_min(TRSACT * T)
{
#ifdef DEBUG_TIMER   
   TIMER("find_min");
#endif
   static TIMER_TYPE start_time = 0;
   start_time = get_time();
   static auto it_remove = T->rows_left.begin();
   static int stat_max = nEdges; // Just some big enough number
   INT min = stat_max;
   INT min_row;
   INT i;
   INT loops_to_do = 0;
#ifdef PRINT_DEBUG
   print_int_arr(g_conform, nNodes, "g_conform");
#endif
   int loop_count = 0;
   for( auto it = T->rows_left.begin(); it != T->rows_left.end() ; ++it)
   {
#ifdef SORT
      // Here is the part that makes this implementation quick:
      // If the conform at time of the latest sort minus..
      // ..loops_after_sort*nCol (which is the max possible change in conform)..
      // ..is bigger than the already found minimum confom, there can't be a lower min..
      // ..therefore we can quit the loop
      if(g_conform[*it] - elems_removed >= min )
         break;
      
#endif
      // Calculate the conform for the current row..
      for( i=0 ; i<T->elem_count[*it] ; i++ )
         T->conform[*it] += T->frq[ T->elem[(*it)][i] ];
      // If the current conform is the minimum, mark this row to be removed
      if( T->conform[*it] < min ) { min = T->conform[*it]; it_remove = it; } 
         ++loops_to_do; 
   }

   min_row = *it_remove;
   ++loops_after_sort;
#ifdef PRINT_DEBUG
   print_int_arr(T->conform, nNodes, "conform");
   printf( "row=%d conform=%d\n", min_row, T->conform[min_row] );
   print_int_arr(&T->elem_buf[T->elem[min_row][0]], T->elem_count[min_row], "eliminated row");
#endif
   if( loops_after_sort % 10 == 0 )
   {
      TIMER_TYPE timeSinceSort = get_time() - g_timeWhenSorted;
      double secondsSinceSort = (get_time() - g_timeWhenSorted ) / (double) g_quadPart;
      double timePerLine =  (get_time() - g_timeWhenSorted ) / (double) (loops_after_sort);
      if ( timePerLine > g_timePerLine )
         g_bSort = true;

      g_timePerLine = timePerLine;
      static TIMER_TYPE interval_time = get_time();
      TIMER_TYPE total_time = get_time() - interval_time;

      double total_seconds = (double)total_time/(double)g_quadPart;
      interval_time = get_time();
      printf( "rows_left=%d seconds=%4.2f conform=%d\n",T->rows_left.size(), total_seconds, T->conform[min_row] );
   }
   // Remove the min row from the vector of rows left      
   T->rows_left.erase( it_remove );

   if ( T->rows_left.empty() )
      return false; //The end

   // Remember the row kicked out
   T->seq[ T->seq_count++ ] = min_row;
   elems_removed += T->elem_count[min_row];
   // Update the frequency table
   for( i=0 ; i<T->elem_count[min_row] ; i++ )
      --T->frq[ T->elem[min_row][i] ];
   
   g_find_time += get_time() - start_time;
   return true;
}

/* Iteration to order the table */
void minus(TRSACT * T)
{
#ifdef DEBUG_TIMER   
   TIMER("minus");
#endif
   loops_after_sort = 0;

   do
   {
#ifdef PRINT_DEBUG   
      print_int_arr(T->frq, nNodes, "frq");
#endif

      if ( ! find_min( T) ) 
         return; // The end
#ifdef SORT      
      if ( g_bSort )
      {
         sort(T);
#ifdef PRINT_DEBUG
         printf("sort/find: %4.2f, threshold=%d\n", double(g_sort_time)/double(g_find_time), SORT_THRESHOLD);
#endif
         //calculate_sort_threshold();
         printf("rows_left=%d sort_time=%4.2f\n", T->rows_left.size(),(double)g_sort_time/(double)g_quadPart);
      }
#endif
       //printf("%d\n", T->rows_left.size() );
      // We always find new conforms, and don't update the old ones using -- etc
      memset( T->conform, 0, sizeof(INT) * nNodes );
   }while(1);

   // Iteration is faster than recursion
   // minus(T); 
}

/* main main*/
int main(int argc, char* argv[])
{
   TRSACT T;
   TIMER_TYPE start_time = 0;
   start_time = get_time();

#ifdef _WIN32
   LARGE_INTEGER timerFreq;
   QueryPerformanceFrequency(&timerFreq);
   g_quadPart = (double)timerFreq.QuadPart;
#endif

#ifndef HARDCODED_DATA
   if( argc < 3 )
   {
      printf("Specify input and output files\n");
      return 1;
   }
   const char * inFileName = argv[1];
   const char * outFileName = argv[2];
#else
   //const char * inFileName = "C:\\Users\\Mikk\\data\\test.txt"; //"C:\\Users\\Mikk\\Dropbox\\git\\MinusTechnique\\data\\chess.dat";
   const char * inFileName = "C:\\Users\\Mikk\\data\\Amazon0302.txt";
   const char * outFileName = "C:\\Users\\Mikk\\data\\amazon06.out";

   // const char * inFileName = "..\\..\\data\\4ta2.txt";
   // const char * outFileName = "..\\..\\data\\4ta2.out";
#endif

   TRSACT_file_load(&T, inFileName);

   TRSACT_init(&T);

   minus(&T);

   // Because didn't want to program a separate minus function for doing the horizontal removal..
   // ..we have this switch that will fake the data a bit
   // TRSACT_switch(&T);

 //  minus(&T);

   // TRSACT_output_result(&T, outFileName);

   TRSACT_free(&T);

   TIMER_TYPE total_time = get_time() - start_time;
   double total_seconds = 0;
   // For getting the time in human-readable form (seconds)
#ifdef _WIN32
   total_seconds = (double)total_time/(double)g_quadPart;
#else
   total_seconds = total_time/1000.0;
#endif
   printf("Finished in about %4.2f seconds \n", total_seconds);
#ifdef DEBUG_TIMER
   TimerContainer::dump("minus_timer.txt");
#endif
   return 0;
}


