/**
 * @brief Problem 1, Laboratory of Algorithms and Data Structures.
 * @version Student
 */

// ----- INCLUDED LIBRARIES ----- //

// Standard input-output library (e.g., fprintf).
#include <stdio.h>
// Time library (e.g., time, clock()).
#include <time.h>
// Standard library (e.g., rand, srand).
#include <stdlib.h>
// Boolean library (e.g., bool).
#include <stdbool.h>
// String library (e.g., memcpy, strcmp)
#include <string.h>

// ----- End INCLUDED LIBRARIES ----- //

// ----- AUXILIARY DATA STRUCTURES ----- //

/**
 * @brief Enumeration data type for the output.
 */
typedef enum {
    ONCONSOLE,  // On console.
    ONFILE      // On file.
} outputEnumType;

/**
 * @brief Pair data structure.
 */
typedef struct {
    clock_t time;   // Time needed to sort.
    bool isSorted;  // Flag representing if the algorithm performed correctly its task.
} pairType;

// ----- End AUXILIARY DATA STRUCTURES ----- //

// ----- GLOBAL VARIABLES ----- //

// IMPORTANT: these constants are only for didactic purposes: you must find the correct values!

// Seed (important for reproducibility).
time_t SEED = 20;
// Minimum size of the array.
const int minSize = 100;
// Maximum size of the aray.
const int maxSize = 1000;
// Number of experiments.
const int numExperiments = 150;
// Granularity of the experiment.
const int granularity = 250;
// Maximum random integer allowed when generating random numbers.
const int maxRandInt = 1000000;
// Thereshold parameter for the base case of HybridSort.            IMPORTANT: this is the result of the first part of the experiment!
const int threshold = 75;
// Output type.
const outputEnumType outputType = ONFILE;
// Output pointer (for printing).
FILE* outputPointer;

// ----- End GLOBAL VARIABLES ----- //

// ----- AUXILIARY FUNCTIONS ----- //

/**
 * @brief Generate a collection of random numbers into an array A of size n.
 * @param A Array of random numbers.
 * @param n Size of the array.
 */
void generateRandomArray(int* A, const int n) {
    // For each i in 0..n-1, generate a random number in the interval [1,maxRandInt].
    for (int i = 0; i < n; i++) A[i] = rand() % maxRandInt + 1;
}

/**
 * @brief Display the array A of size n.
 * @param A Array to be displayed.
 * @param n Size of the array.
 */
void displayArray(const int* A, const int n) {
    // For each i in 0..n-1, print A[i].
    for (int i = 0; i < n; i++) printf("%d ", A[i]);
    printf("\n");
}

// ----- End AUXILIARY FUNCTIONS ----- //

// ----- ANTAGONISTIC FUNCTIONS ----- //

/**
 * @brief Unit test: check if the input array A of size n is sorted.
 * @param A Array to be checked if sorted.
 * @param n Size of the array.
 * @return true if it is sorted; otherwise, false
 */
bool isSorted(const int* A, const int n) {
    // For each i in 0..n-2, if the current element is greater than the next one,
    // then it is unsorted.
    for (int i = 0; i < n-1; i++) if (A[i] > A[i+1]) return false;
    // Otherwise it is.
    return true;
}

// ----- End ANTAGONISTIC FUNCTIONS ----- //

// ----- CORE FUNCTIONS ----- //

/**
 * @brief InsertionSort algorithm.
 * @param A Array of random numbers to be sorted.
 * @param low Left-end index of the array.
 * @param high Right-end index of the array.
 * @property It takes time O(n^2), where n=high-low+1 is the size of the input array.
 */
void insertionSort(int* A, const int low, const int high) {

  int key,k;
  for(int i=low+1;i<=high;i++){
     key=A[i];
     k=i-1;
     while(k>low-1 && A[k]>key){
         A[k+1]=A[k];
         k--;
     }
     A[k+1]=key;
  }
}

/**
 * @brief Merge algorithm.
 * @param A Array to be merged.
 * @param low Left-end index of the array.
 * @param mid Mid index of the array.
 * @param high Right-end index of the array.
 * @property It takes O(n), where n=high-low+1 is the size of the input array.
 */
void merge(int* A, const int low, const int mid, const int high) {

     int n1 = mid - low +1;
     int n2 = high - mid;

     int L[n1], R[n2];

     for(int i=0;i<n1;i++){
        L[i] = A[low + i -1];
     }

     for(int j=0;j<n2;j++){
        R[j] = A[mid+j];
     }

     int i=0, j=0;

     for(int k=low;k<=high;k++){

        if(i<=n1){
            if(j<=n2){
                if(L[i]<=R[j]){
                    A[k] = L[i];
                    i++;
                }else{
                    A[k] = R[j];
                    j++;
                }
            }else{
                 A[k] = L[i];
                 i++;
            }
        }else{
            A[k] = R[j];
            j++;
        }

    }

}


/**
 * @brief MergeSort algorithm.
 * @param A Array of random numbers to be sorted.
 * @param low Left-end index of the array.
 * @param high Right-end index of the array.
 * @property It takes O(n*logn), where n=high-low+1 is the size of the input array.
 */
void mergeSort(int* A, const int low, const int high) {
    if(low<high){
         int middle=(low+high)/2;
         mergeSort(A,low,middle);
         mergeSort(A,middle+1,high);
         merge(A,low,middle,high);
    }

}

/**
 * @brief HybridSort algorithm.
 * @param A Array of random numbers to be sorted.
 * @param low Left-end index of the array.
 * @param high Right-end index of the array.
 * @property It takes O(n*logn), where n=high-low+1 is the size of the input array.
 */
void hybridSort(int* A, const int low, const int high) {

     if(high-low +1 <threshold)//if the array is small then use insertion sort
        insertionSort(A,low,high);
     else{
            int middle=(low+high)/2;
            hybridSort(A,low,middle);
            hybridSort(A,middle+1,high);
            merge(A,low,middle,high);
     }

}

/**
 * @brief Polymorphic function that calls different sorting algorithms.
 * @param randomArray Array of random numbers to be sorted.
 * @param dim Dimension of the input array.
 * @param algo Sorting algorithm to be called. The possible values are: insertionSort,
 *             mergeSort and hybridSort.
 * @return Pair containing the total time needed to sort and the isSorted flag.
 */
pairType sortArray(const int* randomArray, const int dim, const char* algo) {
    // Initiliazation of a pairType with values time = 0 and isSorted = true.
    pairType pair = {0, true};

    // Start and end time.
    clock_t startTime, endTime = 0;

    // Allocate memory for dim integers.
    int* sliceRandomArray = malloc(dim * sizeof(int));
    // Put every i-th element of randomArray into sliceRandomArray for each i in [0, ..., dim-1].
    for (int i = 0; i < dim; i++) sliceRandomArray[i] = randomArray[i];

    // Start the clock.
    startTime = clock();
    // Use InsertionSort.
    if (strcmp(algo, "insertionSort") == 0) insertionSort(sliceRandomArray, 0, dim-1);
    // Use MergeSort.
    else if (strcmp(algo, "mergeSort") == 0) mergeSort(sliceRandomArray, 0, dim-1);
    // Use HybridSort.
    else if (strcmp(algo, "hybridSort") == 0) hybridSort(sliceRandomArray, 0, dim-1);
    // Error
    else {
        fprintf(stderr, "ERROR: There is no such sorting algorithm called %s\n", algo);
        exit(1);
    }
    // Stop the clock.
    endTime = clock();

    // Total time needed to sort.
    pair.time = endTime - startTime;
    // Have we sorted the instance? If not, then the flag is set to false; otherwise, it remains true.
    if (!isSorted(sliceRandomArray, dim)) pair.isSorted = false;

    // Free sliceRandomArray.
    free(sliceRandomArray);

    return pair;
}

// ----- End CORE FUNCTIONS ----- //

// ----- MAIN FUNCTION ----- //

/**
 * @brief Main function.
 * @return Exit code 0.
 */
int main() {
    // Initialize the random seed only once.
    srand(SEED);

    // Accumulated times.
    clock_t timeIS = 0;         // InsertionSort
    clock_t timeMS = 0;         // MergeSort
    clock_t timeHS = 0;         // HybridSort
    // Flags saying either that the output of the execution is correct or not.
    bool isSortedIS = true;     // InsertionSort
    bool isSortedMS = true;     // MergeSort
    bool isSortedHS = true;     // HybridSort
    // Pairs for each of the sorting algorithms.
    pairType pairIS;            // InsertionSort
    pairType pairMS;            // MergeSort
    pairType pairHS;            // HybridSort

    // Allocate an array of maxSize*sizeof(int) cells on the heap.
    // We use this array as a container.
    int* randomArray = malloc(maxSize * sizeof(int));

    // What is the outputPointer?
    if (outputType == ONCONSOLE || outputType == ONFILE) {
        // On console.
        if (outputType== ONCONSOLE) outputPointer = stdout;
        // On file.
        else {
            // Open file.
            outputPointer = fopen("results.txt", "w");
            // Have we opened the file?
            if (outputPointer == NULL) {
                fprintf(stderr, "Error: The outputPointer has not been created\n");
                exit(1);
            }
        }
    }
    // Error
    else {
        fprintf(stderr, "Error: The outputType can be only ONCONSOLE or ONFILE\n");
        exit(1);
    }

    // // Print the header, only if it is on console.
    if (outputType == ONCONSOLE) {
        fprintf(outputPointer, "+-----------+-------------------------------+-------------------------------+-------------------------------+\n");
        fprintf(outputPointer, "| ######### | InsertionSort                 | MergeSort                     | HybridSort                    |\n");
        fprintf(outputPointer, "+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+\n");
        fprintf(outputPointer, "| Dimension | Time              | isSorted? | Time              | isSorted? | Time              | isSorted? |\n");
        fprintf(outputPointer, "+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+\n");
    }

    // Going from minSize to maxSize with step equal to granularity.
    for (int dim = minSize; dim <= maxSize; dim+=granularity) {
        // Reset the accumulated times from one experiment to another.
        timeIS = 0; // InsertionSort
        timeMS = 0; // MergeSort
        timeHS = 0; // HybridSort
        // Reset the isSorted flag for InsertionSort from one experiment to another.
        // We set this flag to true at the beginning, then if at least one time
        // the InsertionSort algorithm fails to sort the input such flag will be
        // set to false; otherwise, it remains true. Similarly for the other algorithms.
        isSortedIS = true;  // InsertionSort
        isSortedMS = true;  // MergeSort
        isSortedHS = true;  // HybridSort

        // Repeat the experiment a numExperiments times for the fixed size (dim).
        for (int exper = 0; exper < numExperiments; exper++) {
            // Fill the array with (pseudo-) random numbers. That is, initialize only the
            // prefix of size dim (<= maxSize) with random numbers.
            generateRandomArray(randomArray, dim);

            // InsertionSort.
            pairIS = sortArray(randomArray, dim, "insertionSort");
            timeIS += pairIS.time;
            isSortedIS = pairIS.isSorted;

            //MergeSort
            pairMS= sortArray(randomArray, dim, "mergeSort");
            timeMS += pairMS.time;
            isSortedMS = pairMS.isSorted;

            //HybridSort
            pairHS = sortArray(randomArray, dim, "hybridSort");
            timeHS += pairHS.time;
            isSortedHS = pairHS.isSorted;


        }
        // Printing the (sample mean as) result. Use TAB (\t) on file.
        if (outputType == ONCONSOLE)
            fprintf(outputPointer, "| %9d | %17f | %9s | %17f | %9s | %17f | %9s |\n",
                dim,
                (double) timeIS/numExperiments, isSortedIS ? "true" : "false",      // InsertionSort
                (double) timeMS/numExperiments, isSortedMS ? "true" : "false",      // MergeSort
                (double) timeHS/numExperiments, isSortedHS ? "true" : "false");     // HybridSort
        else
            fprintf(outputPointer, "%9d\t %17f\t %9s\t %17f\t %9s\t %17f\t %9s\n",
                dim,
                (double) timeIS/numExperiments, isSortedIS ? "true" : "false",      // InsertionSort
                (double) timeMS/numExperiments, isSortedMS ? "true" : "false",      // MergeSort
                (double) timeHS/numExperiments, isSortedHS ? "true" : "false");     // HybridSort
    }

    // Print the ending part, only if it is on console.
    if (outputType == ONCONSOLE)
        fprintf(outputPointer, "+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+\n");

    // Free the allocated memory.
    free(randomArray);

    // If the output is on file, we need to close the file.
    if (outputType == ONFILE) fclose(outputPointer);

    // We are done.
    return 0;
}

// ----- End MAIN FUNCTION ----- //
