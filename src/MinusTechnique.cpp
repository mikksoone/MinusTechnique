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
#include <atomic>
#include <thread>
#include <future>
#include <climits>
// Neat profiler that currently isn't available for everyone
// Therefore it's commented out for the git
// #include "../ProfilingTimer/ProfilingTimer.h"
// #define DEBUG_TIMER

// Verbal version for debugging
// #define PRINT_DEBUG
//#define DEBUG_STATUS_TO_FILE

#define SORT          // Could be undef if one want's to compare time differences
#define MIN_ITEM 0    // If lowest item is 1, lower every item by 1
#define INT int       // Maybe a double is needed?
#define MAX_THREADS 4
//#define HARDCODED_DATA // If we don't want to specify input/output files 

// The timers in windows and linux are different
#ifdef _WIN32
#define TIMER_TYPE __int64
#define int64_t __int64
#else
#define TIMER_TYPE double
#endif

typedef struct {
   INT *buf;   // transactions, its buffer, and their sizes
   INT **frq;
   INT *frq_buf;
   INT *col_max;
   INT *conform;
   INT *seq; //sequence of rows eliminated
   INT seq_count;
   std::vector<INT> rows_left;
} TRSACT;

// Anonymous namespace for global data
namespace
{
   // General data
   int      FILE_err          = 0;
   int      loops_after_sort  = 0;
   int      elems_removed     = 0;
   int      nRows             = 0;
   int      nCols             = 0;
   INT      *g_conform        = 0;
   double   g_quadPart        = 0;
   double   g_timeForConform  = 0;
   double   g_timeForTotalConform=0;
   double   g_totalSortTime    = 0;
   double   g_timePerLine     = LLONG_MAX;
   bool     g_bSort           = 0;
   int g_average_break_point  = 0;
   FILE     *debug            = 0;
   // Timer data
   TIMER_TYPE g_sort_time     = 0;
   TIMER_TYPE g_timeWhenSorted= 0;
   TIMER_TYPE start_time      = 0;
}

static inline void sort(TRSACT * T);

/* Helper function for debugging INT arrays 
   TODO: change %d to %f if needed */
void print_int_arr(INT *arr, INT len, const char *name){
   INT i;
   //return;
   for( i=0 ; i<len ; ++i){
      printf("%d ", arr[i]);
   }
   printf ("==%s==========\n", name);
}

/* allocate memory with error routine (Takeaki Uno)*/
char *alloc_memory (size_t siz){
  char *p = (char *)calloc (siz, 1);
  if ( p == NULL ){
    printf ("out of memory %d\n", GetLastError());
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
  for ( i=(*org_siz)*sizof ; i<siz*sizof ; ++i ) p[i] = 0;
  *org_siz = siz;
  return (p);
}

/* Compare routine for partial sort, currently not used */
int psort_cmp_conform (const INT x, const INT y){
  return (g_conform[x] < g_conform[y]);
}
/* Compare routine for qsort, reorders row items based on their conform */
int qsort_cmp_conform (const void *x, const void *y){
  if ( g_conform[*((INT *)x)] <= g_conform[*((INT *)y)] ) return (-1);
  else return ( g_conform[*((INT *)x)] > g_conform[*((INT *)y)] );
}

class SortFunc
{
public:
   bool operator()(int left, int right) const
   {
      return (g_conform[left] < g_conform[right]);
   }
};
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
   
   if ( !fp ){ printf ("file open error %d\n ", GetLastError() ); exit (1); }
   
   nRows = FILE_read_int (fp);
   if ( (FILE_err&4) != 0){ printf ("file structure error 1\n"); exit (1); }
   nCols = FILE_read_int (fp);
   if ( (FILE_err&4) != 0){ printf ("file structure error 2\n"); exit (1); }
   
   T->buf = (INT *)alloc_memory ( sizeof(INT) * (nRows*nCols) );
   
   j = 0; //row id
   do 
   {
      i = 0; //col id
      do 
      {
         item = FILE_read_int (fp);
         if ( (FILE_err&4) == 0)
         {  // got an item from the file before reaching to line-end
            T->buf[j*nCols+i] = item-MIN_ITEM; //to get the elements started from 0
            i++;
         }
      } while ((FILE_err&3)==0);
      j++;
  } while ( (FILE_err&2)==0);

   fclose(fp);
}

/* Prints the reordered table to a file */
void TRSACT_output_result(TRSACT * T, const char *fname)
{
#ifdef DEBUG_TIMER   
   TIMER("output file");
#endif
   INT i, j;
   FILE *fp = fopen (fname,"w");

   if ( !fp ){ printf ("file open error\n"); exit (1); }

   for( i=0 ; i<nCols ; i++ )
   {
      for( j=0 ; j<nRows ; j++) 
      {
         fprintf( fp, "%d ", T->buf[ T->seq[j]*nCols+i ] + MIN_ITEM );
      }
      fprintf(fp, "\n");
   }
   fclose(fp);
}

/* Frees the transaction and global data */
void TRSACT_free( TRSACT * T )
{
#ifdef DEBUG_TIMER   
   TIMER("free");
#endif
   free( T->buf );
   free( T->col_max );
   free( T->conform );
   free( T->seq );

   free( g_conform );

   free( T->frq_buf );
   free( T->frq );
}

/* Alloc the neccessary memory for transactions and do the initial calculations */
void TRSACT_init( TRSACT * T )
{
#ifdef DEBUG_TIMER   
   TIMER("init");
#endif
   INT i,j;
   T->col_max = (INT*) alloc_memory( sizeof( INT ) * nCols );

   for( j=0 ; j<nRows ; j++)
      for( i=0 ; i<nCols ; i++ )
         if ( T->buf[j*nCols + i] + 1 > T->col_max[i] )
            T->col_max[i] = T->buf[j*nCols + i]+1; // Update the maximum item of col i

   T->frq_buf = (INT*) alloc_memory ( sizeof(INT) * nCols*nRows );
   T->frq = (INT**) alloc_memory( sizeof(INT) * nCols );
   
   size_t cnt = 0;
   for( i=0; i<nCols ; i++ )
   {
      // We just keep track and don't allocate for every col: (INT*)alloc_memory( sizeof(INT) * T->col_max[i] );
      T->frq[i] = &T->frq_buf[cnt]; 
      cnt += T->col_max[i];
   }
   // Fill the frequency table
   for( j=0 ; j<nRows ; j++)
      for( i=0 ; i<nCols ; i++ )
         ++T->frq[i][ T->buf[j*nCols + i] ];


   T->conform      = (INT*) alloc_memory( sizeof(INT) * nRows );
   T->seq         = (INT*) alloc_memory( sizeof(INT) * nRows );
   T->seq_count   = 0;


   g_conform = (INT*) alloc_memory( sizeof(INT) * nRows );

   
   // We add number from 0 to nRow to the rows_left vector
   // These number, of course, represent the rows that have not been kicked out just yet
   T->rows_left.reserve(nRows);
   for( i=0 ; i<nRows ; i++ )
      T->rows_left.push_back(i);

   sort(T);
}

/* Here we transform the file buffer from vertical to horizontal */
void TRSACT_switch(TRSACT * T)
{
#ifdef DEBUG_TIMER   
   TIMER("switch");
#endif
   printf("Switching..\n");
   INT i, j, cnt = 0;
   // Fill the tmp container with values from ordered buf 
   INT * tmp = (INT*) alloc_memory( sizeof(INT) * nRows*nCols );
   for( i=0 ; i<nCols ; i++ )
      for( j=0 ; j<nRows ; j++)   
         tmp[cnt++] = T->buf[ T->seq[j]*nCols +i ]; 
   
#ifdef PRINT_DEBUG
   for( j=0; j<nRow ; j++ )
      print_int_arr(&T->buf[T->seq[j]*nCol], nCol, "buf ordered");
#endif
   // Free the currently used transaction container   
   TRSACT_free(T);

   T->buf = tmp;
   INT k = nRows;
   nRows = nCols;
   nCols = k;

#ifdef PRINT_DEBUG
   for( j=0 ; j<nRow ; j++ )
      print_int_arr(&T->buf[j*nCol], nCol, "buf switched");
#endif

   // Init again..
   TRSACT_init(T);
   printf("Switching done\n");
}


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
   memset( &g_conform[0], 0, sizeof(INT) * nRows );
   // Calcualte conforms
#ifdef MAX_THREADS
   Concurrency::parallel_for_each(T->rows_left.begin(), T->rows_left.end(), [&](int row)
#else
   std::_For_each(T->rows_left.begin(), T->rows_left.end(), [&](int row)
#endif
	{
      for( INT i=0 ; i<nCols ; ++i )
         g_conform[row] += T->frq[i][ T->buf[(row)*nCols +i] ]; 
   });
   // Sort the rows
   // And here we sort the rows according to their conform values
   bool bSorted = false;
#ifdef _WIN32
   if( nRows < 1000000 )
   {
      Concurrency::parallel_sort(T->rows_left.begin(), T->rows_left.end(), SortFunc() );
      bSorted = true;
   }
#endif
   //std::qsort(&T->rows_left[0], T->rows_left.size(), sizeof(INT), qsort_cmp_conform);
   if( !bSorted )
      std::_Insertion_sort( T->rows_left.begin(), T->rows_left.end(), psort_cmp_conform );
   g_bSort = false;
   g_timePerLine = LLONG_MAX;
   loops_after_sort=0;
   elems_removed = 0;
   g_timeForConform=0;
   g_sort_time = get_time() - start_time;
   g_totalSortTime += (double)g_sort_time/(double)g_quadPart;
}

inline int calculate_conform(int j, int increment, int size, TRSACT * T, int min)
{
   if( j >= size )
      return min;
   int delta_check = j+increment*10;
   for( ; j < size ; j+=increment )
   {
      auto row = T->rows_left[j];
#ifdef SORT
      // Here is the part that makes this implementation quick:
      // If the conform at time of the latest sort minus..
      // ..loops_after_sort*nCol (which is the max possible change in conform)..
      // ..is bigger than the already found minimum confom, there can't be a lower min..
      // ..therefore we can quit the loop
      if( j % delta_check == 0 )
      {
         if(g_conform[row] - loops_after_sort*nCols > min )
         {
            //fprintf(debug, "%d\n", j+1);
            return min;
         }
      }
#endif
      // Calculate the conform for the current row..
      for( int i=0 ; i<nCols ; i++ )
         T->conform[row] += T->frq[i][ T->buf[(row)*nCols +i] ];

      // If the current conform is the minimum, mark this row to be removed
      if( T->conform[row] < min )
      {
         min = T->conform[row]; 
      }
   }
   return min;
}

inline 
/* Routine for finding the row with minimum conform */
static inline bool find_min(TRSACT * T)
{
#ifdef DEBUG_TIMER   
   TIMER("find_min");
#endif
   int min = LONG_MAX;
   INT min_row;
   INT i;
   INT loops_to_do = 0;
#ifdef PRINT_DEBUG
   print_int_arr(g_conform, nRow, "g_conform");
#endif
   auto size = T->rows_left.size();
   
   TIMER_TYPE time = get_time();
   auto it_remove = T->rows_left.begin();

#ifndef MAX_THREADS
   min = calculate_conform(0,1,size, T, min);
#else
#ifdef SORT   
   std::vector<std::future<int>> futures;
   for(i=0; i<MAX_THREADS; ++i)
   {
      futures.push_back( std::async(std::launch::async, calculate_conform, i, MAX_THREADS, size, T, min) );
   }
   //Wait for the threads to complete
   for(i=0; i<MAX_THREADS; ++i)
   {
      futures[i].wait();
   }
   //Get the results
   for(i=0; i<MAX_THREADS; ++i)
   {
      int tmp = futures[i].get();
      min = min < tmp ? min : tmp; 
   }
#else
   Concurrency::parallel_for_each(  T->rows_left.begin() , T->rows_left.end() , [T](int row)
   {
      // Calculate the conform for the current row..
      for( int i=0 ; i<nCols ; i++ )
         T->conform[row] += T->frq[i][ T->buf[(row)*nCols +i] ];
   });

   auto last = T->rows_left.end();
   auto first = T->rows_left.begin();

   while (++first!=last)
      if( T->conform[*first] <  T->conform[*it_remove] )
         it_remove=first;
#endif
#endif

#ifdef DEBUG_STATUS_TO_FILE
   fprintf(debug, "g_timeForConform=%.8f\n", g_timeForConform);
#endif

#ifndef MAX_THREADS
   for ( ;it_remove!=T->rows_left.end() ; it_remove++)
   {
      if ( T->conform[*it_remove]==min ) break;
   }
#endif

   g_timeForConform = ((double)(get_time()-time)/(double)g_quadPart); //-g_nThreads*g_threadCreateTime;
   g_timeForTotalConform += g_timeForConform; //+g_nThreads*g_threadCreateTime;

   min_row = *it_remove;
//   printf("min_row=%d conform=%d\n", min_row, T->conform[min_row]);
   ++loops_after_sort;
#ifdef PRINT_DEBUG
   print_int_arr(T->conform, nRow, "conform");
   printf( "row=%d conform=%d\n", min_row, T->conform[min_row] );
   print_int_arr(&T->buf[min_row*nCol], nCol, "eliminated row");
#endif
   //if( loops_after_sort == 0 )
   {
      double timePerLine =  (get_time() - g_timeWhenSorted ) / (double) (loops_after_sort);
      if ( timePerLine > g_timePerLine )
         g_bSort = true;
      //if( g_break_count % 100 == 0)
      //   g_bSort = true;
      g_timePerLine = timePerLine;

      static int cnt = 0;
	   if( ++cnt % 10000 == 0 )
      {
         double total_time = (get_time() - start_time ) / (double) g_quadPart;
         printf( "rows_left=%d seconds=%4.2f conform=%d, sort_time=%4.2f\n",T->rows_left.size(), total_time, T->conform[min_row], (double)g_sort_time/(double)g_quadPart );
      }
   }
   // Remove the min row from the vector of rows left      
   T->rows_left.erase( it_remove );

   // Remember the row kicked out
   T->seq[ T->seq_count++ ] = min_row;

   if ( T->rows_left.empty() )
      return false; //The end

   // Update the frequency table
   for( i=0 ; i<nCols ; i++ )
      --T->frq[i][ T->buf[min_row*nCols + i] ];
   fflush (stdout);
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
   for( INT i=0 ; i<nCols ; i++ )
      print_int_arr(T->frq[i], T->col_max[i], "frq");
#endif

      if ( ! find_min( T) ) 
         return; // The end
#ifdef SORT      
      if ( g_bSort )
      {
	      int tmp = loops_after_sort;
         sort(T);
#ifdef PRINT_DEBUG
         printf("sort/find: %4.2f, threshold=%d\n", double(g_sort_time)/double(g_find_time), SORT_THRESHOLD);
#endif
#ifdef DEBUG_STATUS_TO_FILE
		   fprintf(debug, "%d\t%d\t%4.2f\n", T->rows_left.size(), tmp, total_seconds );
         fflush(debug);
#endif
      }
#endif
      // We always find new conforms, and don't update the old ones using -- etc
      memset( T->conform, 0, sizeof(INT) * nRows );
   }while(1);

   // Iteration is faster than recursion
   // minus(T); 
}

//dummy thread to measure thread execution time
void test_thread()
{

}

void global_init()
{
   TIMER_TYPE start_time = get_time();
#ifdef _WIN32
   LARGE_INTEGER timerFreq;
   QueryPerformanceFrequency(&timerFreq);
   g_quadPart = (double)timerFreq.QuadPart;
#endif
#ifdef DEBUG_STATUS_TO_FILE
   debug = fopen ("debug.out","w");
   if( !debug ) { printf("file open err\n"); exit(1); }
	fprintf(debug, "start..\n");
	fflush(debug);
#endif
  
}

/* main main*/
int main(int argc, char* argv[])
{
   TRSACT T;
   start_time = get_time();
   global_init();

#ifndef HARDCODED_DATA
   if( argc < 3 )
   {
      printf("Specify input and output files\n");
      return 1;
   }
   const char * inFileName = argv[1];
   const char * outFileName = argv[2];
#else
   const char * inFileName = "C:\\Users\\soonem\\Dropbox\\data\\connect_monsa.dat";
   const char * outFileName = "C:\\Users\\soonem\\Dropbox\\data\\connect_monsa_2.out";
#endif
   printf("\nStarting: %s\n", inFileName);
   TRSACT_file_load(&T, inFileName);

   TRSACT_init(&T);

   minus(&T);

   // Because didn't want to program a separate minus function for doing the horizontal removal..
   // ..we have this switch that will fake the data a bit
   TRSACT_switch(&T);

   minus(&T);

   TRSACT_output_result(&T, outFileName);

   TRSACT_free(&T);

   TIMER_TYPE total_time = get_time() - start_time;
   double total_seconds = 0;
   // For getting the time in human-readable form (seconds)
#ifdef _WIN32
   total_seconds = (double)total_time/(double)g_quadPart;
#else
   total_seconds = total_time/1000.0;
#endif
   printf("\n%s\t%6.4f\t%6.4f\t%6.4f\n",inFileName,total_seconds, g_timeForTotalConform, g_totalSortTime);
#ifdef DEBUG_TIMER
   TimerContainer::dump("minus_timer.txt");
#endif
#ifdef DEBUG_STATUS_TO_FILE
   fclose(debug);
#endif
   return 0;
}


