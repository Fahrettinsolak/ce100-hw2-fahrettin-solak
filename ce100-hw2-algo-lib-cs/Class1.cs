using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ce100_hw2_algo_lib_cs
{

    /**
  * @file CE100-hw2-algo-lib-cs
  * @author Fahrettin SOLAK
  * @date 31 March 2022
  *
  * @brief <b> HW-2 Functions </b>
  *
  * HW-2 Sample Lib Functions
  *
  * @see http://bilgisayar.mmf.erdogan.edu.tr/en/
  *
  */



    public class Class1
    {
        //HEAP SORT ALGORITHM

        /**
        *
        *	  @name   Heap Sort ()
        *
        *	  @brief Heap Sort Algorithm
        *
        *  	 Heap sort is a comparison-based sorting technique based on Binary Heap data structure.
        *  	 It is similar to selection sort where we first find the minimum element and place the minimum element at the beginning.
        *     
        *     @param  [in] mass [\arr,n int]   Heap Sort in the serie
        *     @retval [\arr,n int] MATRİX MULTIPLICATION RECURSIVE
        *	  
        **/

        public static void heapSort(int[] arr, int n)
        {
            for (int i = n / 2 - 1; i >= 0; i--)
                heapify(arr, n, i);
            for (int i = n - 1; i >= 0; i--)
            {
                int temp = arr[0];
                arr[0] = arr[i];
                arr[i] = temp;
                heapify(arr, i, 0);
            }
        }
        static void heapify(int[] arr, int n, int i)
        {
            int largest = i;
            int left = 2 * i + 1;
            int right = 2 * i + 2;
            if (left < n && arr[left] > arr[largest])
                largest = left;
            if (right < n && arr[right] > arr[largest])
                largest = right;
            if (largest != i)
            {
                int swap = arr[i];
                arr[i] = arr[largest];
                arr[largest] = swap;
                heapify(arr, n, largest);
            }
        }


        //RADIX SORT
        /**
        *
        *	  @name   Radix Sort ()
        *
        *	  @brief Radix Sort Algorithm
        *
        *	  In computer science, radix sort is a non-comparative sorting algorithm. 
        *	  It avoids comparison by creating and distributing elements into buckets according to their radix. 
        *	  For elements with more than one significant digit, this bucketing process is repeated for each digit, while preserving the ordering of the prior step, until all digits have been considered.
        *
        *
        *       @param  [in] mass [\arr,n int]   Radix Sort in the serie
                @retval [\arr,n int] Radix Sort
        *	  
        **/

        public static int getMax(int[] arr, int n)
        {
            int mx = arr[0];
            for (int i = 1; i < n; i++)
                if (arr[i] > mx)
                    mx = arr[i];
            return mx;
        }
        public static void countSort(int[] arr, int n, int exp)
        {
            int[] output = new int[n]; 
            int i;
            int[] count = new int[10];

            for (i = 0; i < 10; i++)
                count[i] = 0;

            for (i = 0; i < n; i++)
                count[(arr[i] / exp) % 10]++;

            for (i = 1; i < 10; i++)
                count[i] += count[i - 1];

            for (i = n - 1; i >= 0; i--)
            {
                output[count[(arr[i] / exp) % 10] - 1] = arr[i];
                count[(arr[i] / exp) % 10]--;
            }

            for (i = 0; i < n; i++)
                arr[i] = output[i];
        }
        public static void radixsort(int[] arr, int n)
        {
            int m = getMax(arr, n);
            for (int exp = 1; m / exp > 0; exp *= 10)
                countSort(arr, n, exp);
        }



        //COUNTING SORT

        /**
        *
        *	  @name   Counting Sort ()
        *
        *	  @brief Counting Sort Algorithm
        *
        *	  Counting sort is a sorting technique based on keys between a specific range. 
        *	  It works by counting the number of objects having distinct key values (kind of hashing). 
        *	  Then doing some arithmetic to calculate the position of each object in the output sequence.
        *     
        *     @param  [in] mass [\array int]   Counting Sort in the serie
             @retval [\array int ] Counting Sort
        *	  
        **/

        public static void countingsort(int[] Array)
        {
            int n = Array.Length;
            int max = 0;

            for (int i = 0; i < n; i++)
            {
                if (max < Array[i])
                {
                    max = Array[i];
                }
            }

            int[] freq = new int[max + 1];
            for (int i = 0; i < max + 1; i++)
            {
                freq[i] = 0;
            }
            for (int i = 0; i < n; i++)
            {
                freq[Array[i]]++;
            }

            for (int i = 0, j = 0; i <= max; i++)
            {
                while (freq[i] > 0)
                {
                    Array[j] = i;
                    j++;
                    freq[i]--;
                }
            }
        }





        //QUICK SORT LOMUTO

        /**
        *
        *	  @name   Quick Sort Lomuto
        *
        *	  @brief Quick Sort Lomuto Algorithm
        *
        *	  This algorithm works by assuming the pivot element as the last element. 
        *	  If any other element is given as a pivot element then swap it first with the last element.
        *	  Now initialize two variables i as low and j also low, swap whenever arr[I] <= arr[j], and increment i, otherwise only increment j.
        *	  After coming out from the loop swap arr[i] with arr[hi]. This i stores the pivot element.
        *     
        *     @param  [in] mass [\arr,low,high int]   Quick Sort Lomuto in the serie
             @retval [\arr,low,high int] Quick Sort Lomuto
        *	  
        **/

        static void Swap(int[] array,
                 int position1,
                 int position2)
        {
            int temp = array[position1];
            array[position1] = array[position2];
            array[position2] = temp;
        }

        /* This function takes last element as
        pivot, places the pivot element at its
        correct position in sorted array, and
        places all smaller (smaller than pivot)
        to left of pivot and all greater elements
        to right of pivot */
        static int partition(int[] arr, int low,
                                        int high)
        {
            int pivot = arr[high];
            int i = (low - 1);

            for (int j = low; j <= high - 1; j++)
            {
                if (arr[j] <= pivot)
                {
                    i++; 
                    Swap(arr, i, j);
                }
            }
            Swap(arr, i + 1, high);
            return (i + 1);
        }
        public static void quickSortLomuto(int[] arr, int low,
                                         int high)
        {
            if (low < high)
            {
                int pi = partition(arr, low, high);
                quickSortLomuto(arr, low, pi - 1);
                quickSortLomuto(arr, pi + 1, high);
            }
        }



        //QUICK SORT HOARE'S

        /**
        *
        *	  @name   Quick Sort Hoare's
        *
        *	  @brief Quick Sort Hoare's Algorithm
        *
        *	 Hoare's scheme is more efficient than Lomuto's partition scheme because it does three times fewer swaps on average.
        *	 Also, as mentioned, the implementation given creates a balanced partition even when all values are equal.[10][self-published source?], which Lomuto's scheme does not.
        *	 Like Lomuto's partition scheme, Hoare's partitioning also would cause Quicksort to degrade to O(n2) for already sorted input, if the pivot was chosen as the first or the last element.
        *     
        *     @param  [in] mass [\array,position1,position2 int]   Quick Sort Hoare's in the serie
             @retval [\array,position1,position2 int] Quick Sort Hoare's
        *	  
        **/

        static void SwapHoares(int[] array,
                 int position1,
                 int position2)
        {

            int temp = array[position1];
            array[position1] = array[position2];
            array[position2] = temp;
        }
        static int partitionHoares(int[] arr, int low,
                                        int high)
        {
            int pivot = arr[high];
            int i = (low - 1);

            for (int j = low; j <= high - 1; j++)
            {
                if (arr[j] <= pivot)
                {
                    i++;       
                    SwapHoares(arr, i, j);
                }
            }
            SwapHoares(arr, i + 1, high);
            return (i + 1);
        }

        /* The main function that
           implements QuickSort
        arr[] --> Array to be sorted,
        low --> Starting index,
        high --> Ending index */
        public static void quickSortHoares(int[] arr, int low,
                                         int high)
        {
            if (low < high)
            {
                /* pi is partitioning index,
                arr[p] is now at right place */
                int pi = partition(arr, low, high);
                quickSortHoares(arr, low, pi - 1);
                quickSortHoares(arr, pi + 1, high);
            }
        }







        //RANDOM QUICK SORT LOMUTO


        /**
        *
        *	  @name   Random Quick Sort Lomuto
        *
        *	  @brief Random Quick Sort Lomuto Algorithm
        *
        *	 In this article, we will discuss how to implement QuickSort using random pivoting. 
        *	 In QuickSort we first partition the array in place such that all elements to the left of the pivot element are smaller, while all elements to the right of the pivot are greater than the pivot.
        *     
        *     @param  [in] arr [\arr,low,high int]   Quick Sort Hoare's in the serie
             @retval [\arr,low,high int] Quick Sort Hoare's
        *	  
        **/

        /* This function takes last element as pivot,
    places the pivot element at its correct
    position in sorted array, and places all
    smaller (smaller than pivot) to left of
    pivot and all greater elements to right
    of pivot */
        static int partitionRL(int[] arr, int low, int high)
        {
            randomQSL(arr, low, high);
            int pivot = arr[high];

            int i = (low - 1); 
            for (int j = low; j < high; j++)
            {
                if (arr[j] < pivot)
                {
                    i++;
                    int tempp = arr[i];
                    arr[i] = arr[j];
                    arr[j] = tempp;
                }
            }
            int tempp2 = arr[i + 1];
            arr[i + 1] = arr[high];
            arr[high] = tempp2;

            return i + 1;
        }
        static int randomQSL(int[] arr, int low, int high)
        {

            Random rand = new Random();
            int pivot = rand.Next() % (high - low) + low;

            int tempp1 = arr[pivot];
            arr[pivot] = arr[high];
            arr[high] = tempp1;

            return partition(arr, low, high);
        }

        /* The main function that implements Quicksort()
          arr[] --> Array to be sorted,
          low --> Starting index,
          high --> Ending index */
        public static void randomQuickSortLomuto(int[] arr, int low, int high)
        {
            if (low < high)
            {
                int pi = partition(arr, low, high);
                randomQuickSortLomuto(arr, low, pi - 1);
                randomQuickSortLomuto(arr, pi + 1, high);
            }
        }





        //RANDOMIZE QUICK SORT HORAE'S

        /**
        *
        *	     @name  Randomize quick Sort Hoare's
        *   
        *	     @brief Randomize quick Sort Hoare's Algorithm
        *
        *	     In QuickSort we first partition the array in place such that all elements to the left of the pivot element are smaller,while all elements to the right of the pivot are greater than the pivot.
        *     
        *        @param  [in] arr [\array,low,high int]   Randomize quick Sort Hoare's in the serie
        *        @retval [\array,low,high int] Randomize quick Sort Hoare's
        *	  
        **/
        public static int randomhoarepartition(int[] arr, int low, int high)
        {
            int pivot = arr[low];
            int i = low - 1, j = high + 1;

            while (true)
            {
                do
                {
                    i++;
                } while (arr[i] < pivot);

                do
                {
                    j--;
                } while (arr[j] > pivot);

                if (i >= j)
                    return j;

                int tempp = arr[i];
                arr[i] = arr[j];
                arr[j] = tempp;
            }
        }
        public static int Random(int[] arr, int low, int high)
        {
            Random rand = new Random();
            int pivot = rand.Next() % (high - low) + low;

            int tempp1 = arr[pivot];
            arr[pivot] = arr[high];
            arr[high] = tempp1;

            return randomhoarepartition(arr, low, high);
        }
        public static void randomQuicksortHoare(int[] array, int lw, int high)
        {
            if (lw < high)
            {
                int pi = Random(array, lw, high);
                randomQuicksortHoare(array, lw, pi);
                randomQuicksortHoare(array, pi + 1, high);
            }
        }






        //MATRİX MULTIPLICATION RECURSIVE


        /*
        *
        *      @name MATRİX MULTIPLICATION RECURSIVE
        *
        *      @brief MATRİX MULTIPLICATION RECURSIVE Function
        *
        *      First check if multiplication between matrices is possible or not.
        *      For this, check if number of columns of first matrix is equal to number of rows of second matrix or not. 
        *      If both are equal than proceed further otherwise generate output “Not Possible”.
        *
        *      @param[in] mass [\row1,col1,A,row2,col2,B,C int]   MATRİX MULTIPLICATION RECURSIVE in the serie
        *
        *      @retval[\row1,col1,A,row2,col2,B,C  int] MATRİX MULTIPLICATION RECURSIVE
       */

        public static int MAX = 100;
        public static int i = 0, j = 0, k = 0;

        public static void multiplyMatrixRec(int row1, int col1,
                                  int[,] A, int row2,
                                  int col2, int[,] B,
                                  int[,] C)
        {
            if (i >= row1)
                return;

            if (j < col2)
            {
                if (k < col1)
                {
                    C[i, j] += A[i, k] * B[k, j];
                    k++;

                    multiplyMatrixRec(row1, col1, A,
                                      row2, col2, B, C);
                }

                k = 0;
                j++;
                multiplyMatrixRec(row1, col1, A,
                                  row2, col2, B, C);
            }

            j = 0;
            i++;
            multiplyMatrixRec(row1, col1, A,
                              row2, col2, B, C);
        }



        // MATRİX MULTIPLICATION ITERATIVE

        /*
        *
        *      @name MATRİX MULTIPLICATION ITERATIVE
        *
        *      @brief MATRİX MULTIPLICATION ITERATIVE Function
        *      This program can multiply any two square or rectangular matrices.
        *      The below program multiplies two square matrices of size 4 * 4.
        *      There is also an example of a rectangular matrix for the same code (commented below).
        *      We can change the Matrix value with the number of rows and columns (from MACROs) for Matrix-1 and Matrix-2 for different dimensions.
        *      @param[in] mass[\mat1, mat2, res int]   MATRİX MULTIPLICATION ITERATIVE in the serie
        *
        *      @retval[\mat1, mat2, rec int] MATRİX MULTIPLICATION ITERATIVE
        */

        static int N = 4;
        public static void multiply(int[,] mat1,
                         int[,] mat2, int[,] res)
        {
            int i, j, k;
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    res[i, j] = 0;
                    for (k = 0; k < N; k++)
                        res[i, j] += mat1[i, k]
                                     * mat2[k, j];
                }
            }
        }



        //MATRİX MULTIPLICATION STRASSEN

        /*
        *
        *      @name MATRİX MULTIPLICATION STRASSEN
        *
        *      @brief MATRİX MULTIPLICATION STRASSEN Function
        *
        *       In the above divide and conquer method, the main component for high time complexity is 8 recursive calls. 
        *       The idea of Strassen’s method is to reduce the number of recursive calls to 7. Strassen’s method is similar to above simple divide
        *       and conquer method in the sense that this method also divide matrices to sub-matrices of size N/2 x N/2 as shown in the above diagram, 
        *       but in Strassen’s method, the four sub-matrices of result are calculated using following formulae.
        *
        * @param[in] mass[\a,b,a11,a12,b21,b22,region,n int] MATRİX MULTIPLICATION STRASSEN in the serie
        *
        *      @retval[\a,b,a11,a12,b21,b22,region,n int] MATRİX MULTIPLICATION STRASSEN
        */


        static void calculate(int n)
        {
            Random rng = new Random();
            float[,] m1 = new float[n, n];
            float[,] m2 = new float[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    m1[i, j] = (float)rng.NextDouble();
                }
            }

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    m2[i, j] = (float)rng.NextDouble();
                }
            }

            float[,] m3 = strassen(m1, m2, n);
        }

        public static float[,] strassen(float[,] a, float[,] b, int n)
        {
            if (n == 2)
            {
                var m1 = (a[0, 0] + a[1, 1]) * (b[0, 0] + b[1, 1]);
                var m2 = (a[1, 0] + a[1, 1]) * b[0, 0];
                var m3 = a[0, 0] * (b[0, 1] - b[1, 1]);
                var m4 = a[1, 1] * (b[1, 0] - b[0, 0]);
                var m5 = (a[0, 0] + a[0, 1]) * b[1, 1];
                var m6 = (a[1, 0] - a[0, 0]) * (b[0, 0] + b[0, 1]);
                var m7 = (a[0, 1] - a[1, 1]) * (b[1, 0] + b[1, 1]);
                a[0, 0] = m1 + m4 - m5 + m7;
                a[0, 1] = m3 + m5;
                a[1, 0] = m2 + m4;
                a[1, 1] = m1 - m2 + m3 + m6;
                return a;
            }
            else
            {
                float[,] a11 = matrixDivide(a, n, 11);
                float[,] a12 = matrixDivide(a, n, 12);
                float[,] a21 = matrixDivide(a, n, 21);
                float[,] a22 = matrixDivide(a, n, 22);

                float[,] b11 = matrixDivide(b, n, 11);
                float[,] b12 = matrixDivide(b, n, 12);
                float[,] b21 = matrixDivide(b, n, 21);
                float[,] b22 = matrixDivide(b, n, 22);

                float[,] p1 = strassen(a11, matrixDiff(b12, b22, n / 2), n / 2);
                float[,] p2 = strassen(matrixSum(a11, a12, n / 2), b22, n / 2);
                float[,] p3 = strassen(matrixSum(a21, a22, n / 2), b11, n / 2);
                float[,] p4 = strassen(a22, matrixDiff(b21, b11, n / 2), n / 2);
                float[,] p5 = strassen(matrixSum(a11, a22, n / 2), matrixSum(b11, b22, n / 2), n / 2);
                float[,] p6 = strassen(matrixDiff(a12, a22, n / 2), matrixSum(b21, b22, n / 2), n / 2);
                float[,] p7 = strassen(matrixDiff(a11, a21, n / 2), matrixSum(b11, b12, n / 2), n / 2);

                float[,] c11 = matrixDiff(matrixSum(p5, p4, n / 2), matrixDiff(p2, p6, n / 2), n / 2);
                float[,] c12 = matrixSum(p1, p2, n / 2);
                float[,] c21 = matrixSum(p3, p4, n / 2);
                float[,] c22 = matrixDiff(matrixSum(p1, p5, n / 2), matrixSum(p3, p7, n / 2), n / 2);

                for (int i = 0; i < n / 2; i++)
                {
                    for (int j = 0; j < n / 2; j++)
                    {
                        a[i, j] = c11[i, j];
                        a[i, j + n / 2] = c12[i, j];
                        a[i + n / 2, j] = c21[i, j];
                        a[i + n / 2, j + n / 2] = c22[i, j];
                    }
                }
                return a;
            }
        }

        static float[,] matrixSum(float[,] a, float[,] b, int n)
        {
            float[,] c = new float[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    c[i, j] = a[i, j] + b[i, j];
            return c;
        }


        static float[,] matrixDiff(float[,] a, float[,] b, int n)
        {
            float[,] c = new float[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    c[i, j] = a[i, j] - b[i, j];
            return c;
        }

        static float[,] matrixCombine(float[,] a11, float[,] a12, float[,] a21, float[,] a22, int n)
        {
            float[,] a = new float[n, n];
            for (int i = 0; i < n / 2; i++)
            {
                for (int j = 0; j < n / 2; j++)
                {
                    a[i, j] = a11[i, j];
                    a[i, j + n / 2] = a12[i, j];
                    a[i + n / 2, j] = a21[i, j];
                    a[i + n / 2, j + n / 2] = a22[i, j];
                }
            }
            return a;
        }

        static float[,] matrixDivide(float[,] a, int n, int region)
        {
            float[,] c = new float[n / 2, n / 2];
            if (region == 11)
            {
                for (int i = 0, x = 0; x < n / 2; i++, x++)
                {
                    for (int j = 0, y = 0; y < n / 2; j++, y++)
                    {
                        c[i, j] = a[x, y];
                    }
                }
            }
            else if (region == 12)
            {
                for (int i = 0, x = 0; x < n / 2; i++, x++)
                {
                    for (int j = 0, y = n / 2; y < n; j++, y++)
                    {
                        c[i, j] = a[x, y];
                    }
                }
            }
            else if (region == 21)
            {
                for (int i = 0, x = n / 2; x < n; i++, x++)
                {
                    for (int j = 0, y = 0; y < n / 2; j++, y++)
                    {
                        c[i, j] = a[x, y];
                    }
                }
            }
            else if (region == 22)
            {
                for (int i = 0, x = n / 2; x < n; i++, x++)
                {
                    for (int j = 0, y = n / 2; y < n; j++, y++)
                    {
                        c[i, j] = a[x, y];
                    }
                }
            }
            return c;
        }



        //PRIORİTY QUEUE

        /**
        *
        *	  @name  Priority queue
        *
        *	  @brief Priority queue Algorithm
        *
        *	 Priority Queue is an abstract data type, which is similar to a queue, however, in the priority queue, every element has some priority. 
        *	 The priority of the elements in a priority queue determines the order in which elements are removed from the priority queue. 
        *	 Therefore all the elements are either arranged in an ascending or descending order.
        *     
        *     @param  [in] mass [\i,p,j int]   Priority queue in the serie
            @retval [\b n] Priority queue
        *	  
        **/


        static int[] H = new int[50];
        static int sized = -1;
        static int parent(int i)
        {
            return (i - 1) / 2;
        }
        static int leftChild(int i)
        {
            return ((2 * i) + 1);
        }
        static int rightChild(int i)
        {
            return ((2 * i) + 2);
        }
        static void shiftUp(int i)
        {
            while (i > 0 &&
                   H[parent(i)] < H[i])
            {
                swapheap(parent(i), i);

                i = parent(i);
            }
        }
        static void shiftDown(int i)
        {
            int maxIndex = i;

            int l = leftChild(i);

            if (l <= sized &&
                H[l] > H[maxIndex])
            {
                maxIndex = l;
            }

            int r = rightChild(i);

            if (r <= sized &&
                H[r] > H[maxIndex])
            {
                maxIndex = r;
            }

            if (i != maxIndex)
            {
                swapheap(i, maxIndex);
                shiftDown(maxIndex);
            }
        }
        public static void insert(int p)
        {
            sized = sized + 1;
            H[sized] = p;

            shiftUp(sized);
        }
        public static int extractMax()
        {
            int result = H[0];

            H[0] = H[sized];
            sized = sized - 1;

            shiftDown(0);
            return result;
        }
        static void changePriority(int i,
                                   int p)
        {
            int oldp = H[i];
            H[i] = p;

            if (p > oldp)
            {
                shiftUp(i);
            }
            else
            {
                shiftDown(i);
            }
        }
        static int getMaxheap()
        {
            return H[0];
        }
        static void Remove(int i)
        {
            H[i] = getMaxheap() + 1;

            shiftUp(i);
            extractMax();
        }
        static void swapheap(int i, int j)
        {
            int temp = H[i];
            H[i] = H[j];
            H[j] = temp;
        }







        //SELECTION PROBLEM 

        /*
        *
        *   @name SELECTION PROBLEM
        *
        *   @brief  SELECTION PROBLEM Function
        *
        *   The activity selection problem is a combinatorial optimization problem concerning the selection of non-conflicting activities to perform within a given time frame,
        *   given a set of activities each marked by a start time (si) and finish time (fi).
        *   The problem is to select the maximum number of activities that can be performed by a single person or machine, 
        *   assuming that a person can only work on a single activity at a time. 
        *
        *    @param[in] mass [\i,arr,n,l,r,k,j,x,a,b int]   SELECTION PROBLEM in the serie
        *
        *   @retval[\i,arr,n,l,r,k,j,x,a,b int] SELECTION PROBLEM
               */


        /// <summary getMax>
        /// 
        /// </summary Function to get value of
        /// the current maximum element>
        /// <returns></H[0]>
        static int getMax()
        {
            return H[0];
        }


        /// <summary Remove>
        /// 
        /// </summary Function to remove the element
        /// located at given index>
        /// <param name="i"></int>
        static void Removed(int i)
        {
            H[i] = getMax() + 1;

            shiftUp(i);

            extractMax();
        }

        static int findMedian(int[] arr, int i, int n)
        {
            if (i <= n)
                Array.Sort(arr, i, n); 
            else
                Array.Sort(arr, n, i);
            return arr[n / 2]; 
        }

        public static int SelectionProblem(int[] arr, int l,
                                    int r, int k)
        {

            if (k > 0 && k <= r - l + 1)
            {
                int n = r - l + 1; 

                int i;

                int[] median = new int[(n + 4) / 5];
                for (i = 0; i < n / 5; i++)
                    median[i] = findMedian(arr, l + i * 5, 5);

                if (i * 5 < n)
                {
                    median[i] = findMedian(arr, l + i * 5, n % 5);
                    i++;
                }

                int medOfMed = (i == 1) ? median[i - 1] :
                                        SelectionProblem(median, 0, i - 1, i / 2);
                int pos = partitionforselection(arr, l, r, medOfMed);
                if (pos - l == k - 1)
                    return arr[pos];
                if (pos - l > k - 1) 
                    return SelectionProblem(arr, l, pos - 1, k);

                return SelectionProblem(arr, pos + 1, r, k - pos + l - 1);
            }
            return int.MaxValue;
        }

        public static int[] Swapforselection(int[] arr, int i, int j)
        {
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
            return arr;
        }
        static int partitionforselection(int[] arr, int l,
                                int r, int x)
        {
            int i;
            for (i = l; i < r; i++)
                if (arr[i] == x)
                    break;

            Swapforselection(arr, i, r);
            i = l;
            for (int j = l; j <= r - 1; j++)
            {
                if (arr[j] <= x)
                {
                    Swapforselection(arr, i, j);
                    i++;
                }
            }
            Swapforselection(arr, i, r);
            return i;
        }
        int size = 0;
        static void Swap(int a, int b)
        {
            int temp = b;
            b = a;
            a = temp;
        }



    }
}
