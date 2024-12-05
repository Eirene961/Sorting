#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <string.h>
#include <algorithm>
#include <chrono>
using namespace std;

template <class T>
void Swap(T& a, T& b)
{
	T x = a;
	a = b;
	b = x;
}

//-------------------------------------------------

// This function generates a random array
void GenerateRandomData(int a[], int n)
{
	srand((unsigned int)time(NULL));

	for (int i = 0; i < n; i++)
	{
		a[i] = rand() % n;
	}
}

// This function generates a sorted array (ascending order)
void GenerateSortedData(int a[], int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = i;
	}
}

// This function generates a reverse-sorted array (descending order)
void GenerateReverseData(int a[], int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = n - 1 - i;
	}
}

// This function generates a nearly-sorted array
void GenerateNearlySortedData(int a[], int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = i;
	}

	srand((unsigned int)time(NULL));

	for (int i = 0; i < 10; i++)
	{
		int r1 = rand() % n;
		int r2 = rand() % n;
		Swap(a[r1], a[r2]);
	}
}

void GenerateData(int a[], int n, int dataType)
{
	switch (dataType)
	{
	case 0:
		GenerateRandomData(a, n);
		break;
	case 1:
		GenerateSortedData(a, n);
		break;
	case 2:
		GenerateReverseData(a, n);
		break;
	case 3:
		GenerateNearlySortedData(a, n);
		break;
	default:
		printf("Error: unknown data type!\n");
	}
}


unsigned long long count_comparison = 0;
chrono::duration<double> elapsed;
int secsToMilisecs = 1000;

// ---------------------- THE ANH ----------------------

void SelectionSort(int a[], int n) {
	count_comparison = 0;
	for (int i = 0; ++count_comparison && i < n - 1; i++)
	{

		int min = i;
		for (int j = i + 1; ++count_comparison && j < n; j++)
		{
			if (++count_comparison && a[j] < a[min]) {
				min = j;
			}
		}
		swap(a[i], a[min]);
	}
}

void Heapify(int a[], int n, int i)
{
	int largest = i;
	int left = 2 * i + 1;
	int right = 2 * i + 2;

	if (++count_comparison && left < n && a[left]>a[largest])
		largest = left;

	if (++count_comparison && right<n && a[right]>a[largest])
		largest = right;

	if (++count_comparison && largest != i)
	{
		swap(a[i], a[largest]);
		Heapify(a, n, largest);
	}
}

void HeapSort(int a[], int n)
{
	count_comparison = 0;
	for (int i = n / 2 - 1; ++count_comparison && i >= 0; i--)
		Heapify(a, n, i);

	for (int i = n - 1; +count_comparison && i > 0; i--)
	{
		swap(a[0], a[i]);
		Heapify(a, i, 0);
	}
}

void Merge(int a[], int left, int mid, int right)
{
	int n1 = mid - left + 1;
	int n2 = right - mid;
	int* L = new int[n1];
	int* R = new int[n2];
	for (int i = 0; ++count_comparison && i < n1; i++)
		L[i] = a[left + i];
	for (int j = 0; ++count_comparison && j < n2; j++)
		R[j] = a[mid + 1 + j];
	int i = 0;
	int j = 0;
	int k = left;
	while (++count_comparison && i < n1 && j < n2)
	{
		if (++count_comparison && L[i] <= R[j])
		{
			a[k] = L[i];
			i++;
		}
		else
		{
			a[k] = R[j];
			j++;
		}
		k++;
	}
	while (++count_comparison && i < n1)
	{
		a[k] = L[i];
		i++;
		k++;
	}
	while (++count_comparison && j < n2)
	{
		a[k] = R[j];
		j++;
		k++;
	}
	delete[] L;
	delete[] R;
}

void MergeSort(int a[], int left, int right)
{
	if (++count_comparison && left < right)
	{
		int mid = left + (right - left) / 2;
		MergeSort(a, left, mid);
		MergeSort(a, mid + 1, right);
		Merge(a, left, mid, right);
	}
}


// ---------------------- NHAT ANH ----------------------

void InsertionSort(int a[], int n) {
	count_comparison = 0;

	for (int i = 1; ++count_comparison && i < n; i++) {
		int x = a[i];
		int j = i;

		while (++count_comparison && j >= 1 && ++count_comparison && a[j - 1] > x) {
			a[j] = a[j - 1];
			j--;
		}
		a[j] = x;
	}
}


void ShellSort(int a[], int n) {
	count_comparison = 0;

	vector<int> jumps;
	int h = 0, i = 0;
	while (1) {
		h = 3 * i + 1;
		if (h >= n)
			break;
		jumps.push_back(h);
		i++;
	}
	int m = jumps.size();

	for (i = m - 1; ++count_comparison && i >= 0; i--) {
		int jump = jumps[i];
		for (int j = jump; ++count_comparison && j < n; j += jump) {
			int x = a[j];
			int t = j;
			while (++count_comparison && t >= jump && ++count_comparison && a[t - jump] > x) {
				a[t] = a[t - jump];
				t -= jump;
			}
			a[t] = x;
		}
	}
}


int BinarySearchForPosition(int a[], int n, int x) {
	int from = 0, to = n - 1;
	int mid = 0;
	while (++count_comparison && from <= to) {
		mid = (from + to) / 2;
		if (++count_comparison && a[mid] <= x)
			from = mid + 1;
		else
			to = mid - 1;
	}
	return from;
}


void BinaryInsertionSort(int a[], int n) {
	count_comparison = 0;

	for (int i = 1; ++count_comparison && i < n; i++) {
		int x = a[i], k = BinarySearchForPosition(a, i, x);
		for (int j = i; ++count_comparison && j > k; j--)
			a[j] = a[j - 1];
		a[k] = x;
	}
}

// ---------------------- BINH ----------------------
void BubbleSort(int a[], int n) {
	count_comparison = 0;

	for (int i = n - 1; i >= 1; i--) {
		count_comparison++;
		for (int j = 0; j < i; j++) {
			count_comparison++;
			if (++count_comparison && a[j] > a[j + 1])
				swap(a[j], a[j + 1]);
		}
	}
}


void ShakerSort(int a[], int n) {
	count_comparison = 0;

	int h = 0;
	int k = n - 1;
	while (++count_comparison && h < k) {
		count_comparison += 2;

		for (int i = h; i < k; i++) {
			count_comparison++;
			if (a[i] > a[i + 1])
				swap(a[i], a[i + 1]);
		}
		k--;


		for (int j = k; j > h; j--) {
			count_comparison++;
			if (a[j] < a[j - 1])
				swap(a[j], a[j - 1]);
		}
		h++;
	}
}


void QuickSort(int a[], int low, int high) {
	if (++count_comparison && low < high) {
		int mid = low + (high - low) / 2;
		swap(a[mid], a[high]);
		int pivot = a[high];

		int i = low - 1;

		for (int j = low; j < high; j++) {
			count_comparison++;
			if (a[j] < pivot) {
				i++;
				std::swap(a[i], a[j]);
			}
		}
		std::swap(a[i + 1], a[high]);
		int pi = i + 1;

		QuickSort(a, low, pi - 1);
		QuickSort(a, pi + 1, high);
	}
}


// ---------------------- BAU ----------------------

void CountingSort(int a[], int n) {
	count_comparison = 0;

	// Step 1: Find the range of the array
	int max_val = *max_element(a, a + n);
	int min_val = *min_element(a, a + n);
	int range_of_elements = max_val - min_val + 1;

	// Step 2: Initialize count array and output array
	vector<int> output(n);
	vector<int> count(range_of_elements, 0);

	// Step 3: Count occurrences of each element
	for (int i = 0; ++count_comparison && i < n; i++) {
		count[a[i] - min_val]++;
	}

	// Step 4: Update count array with cumulative counts
	for (int i = 1; ++count_comparison && i < range_of_elements; i++) {
		count[i] += count[i - 1];
	}

	// Step 5: Build the output array
	int i = n - 1;
	while (++count_comparison && i >= 0) {
		count[a[i] - min_val]--;
		output[count[a[i] - min_val]] = a[i];
		i--;
	}

	// Step 6: Copy the sorted output array back to the original array
	for (int i = 0; ++count_comparison && i < n; i++) {
		a[i] = output[i];
	}
}

void CountingSortForRadix(int a[], int n, int exp) {
	vector<int> output(n); // Output array
	int count[10]; // Count array for digits 0-9
	memset(count, 0, sizeof(count));

	// Count occurrences of each digit in the current place
	for (int i = 0; ++count_comparison && i < n; i++) {
		int index = (a[i] / exp) % 10;
		count[index]++;
	}

	// Update count[i] to hold the position of the next digit in output
	for (int i = 1; ++count_comparison && i < 10; i++) {
		count[i] += count[i - 1];
	}

	// Build the output array
	int i = n - 1;
	while (++count_comparison && i >= 0) {
		int index = (a[i] / exp) % 10;
		count[index]--;
		output[count[index]] = a[i];
		i--;
	}

	// Copy the output array to the original array
	for (int i = 0; ++count_comparison && i < n; i++) {
		a[i] = output[i];
	}
}

void RadixSort(int a[], int n) {
	count_comparison = 0;

	// Find the maximum number to determine the number of digits
	int max_num = *max_element(a, a + n);
	int exp = 1; // Start with the least significant digit

	// Perform counting sort for each digit
	while (++count_comparison && max_num / exp > 0) {
		CountingSortForRadix(a, n, exp);
		exp *= 10;
	}
}

void FlashSort(int a[], int n) {
	count_comparison = 0;

	if (n <= 1)
		return;
	
	// Step 1: Find the minimum and maximum elements
	int min_val = *min_element(a, a + n);
	int max_val = *max_element(a, a + n);
	if (min_val == max_val)
		return; // Array is already sorted

	// Step 2: Initialize the class boundaries
	int m = (int)(0.45 * n); // Number of classes
	vector<int> class_counts(m, 0);
	double scaling_factor = double(m - 1) / (max_val - min_val);

	// Step 3: Classify elements
	for (int i = 0; ++count_comparison && i < n; i++) {
		int class_index = int(scaling_factor * (a[i] - min_val));
		class_counts[class_index]++;
	}

	// Step 4: Compute cumulative counts (class boundaries)
	for (int i = 1; ++count_comparison && i < m; i++) {
		class_counts[i] += class_counts[i - 1];
	}

	// Step 5: Permute elements into their respective classes
	int count = 0;
	int i = 0;
	while (++count_comparison && count < n) {
		int class_index = int(scaling_factor * (a[i] - min_val));
		if (++count_comparison && i >= class_counts[class_index]) {
			i++;
		}
		else {
			class_counts[class_index]--;
			swap(a[i], a[class_counts[class_index]]);
			count++;
		}
	}

	// Step 6: Sort within each class
	int start = 0;
	for (int idx = 0; ++count_comparison && idx < m; idx++) {
		int end = class_counts[idx];
		if (++count_comparison && end - start > 1) {
			for (int i = start; ++count_comparison && i < end; i++) {
				int key = a[i];
				int j = i - 1;
				while (++count_comparison && j >= start && ++count_comparison && a[j] > key) {
					a[j + 1] = a[j];
					j--;
				}
				a[j + 1] = key;
			}
		}
		start = end;
	}
}


bool isanumber(string s) {
	for (int i = 0; i < s.length(); i++) {
		if (!isdigit(s[i]))
			return false;
	}
	return true;
}

void ExecuteSortAlgorithm(int a[], int n, string algorithm_name) {
	auto start = chrono::high_resolution_clock::now();
	count_comparison = 0;

	if (algorithm_name == "selection-sort") {
		SelectionSort(a, n);
	}
	else if (algorithm_name == "counting-sort") {
		CountingSort(a, n);
	}
	else if (algorithm_name == "radix-sort") {
		RadixSort(a, n);
	}
	else if (algorithm_name == "flash-sort") {
		FlashSort(a, n);
	}
	else if (algorithm_name == "heap-sort") {
		HeapSort(a, n);
	}
	else if (algorithm_name == "merge-sort") {
		MergeSort(a, 0, n - 1);
	}
	else if (algorithm_name == "insertion-sort") {
		InsertionSort(a, n);
	}
	else if (algorithm_name == "shell-sort") {
		ShellSort(a, n);
	}
	else if (algorithm_name == "binary-insertion-sort") {
		BinaryInsertionSort(a, n);
	}
	else if (algorithm_name == "bubble-sort") {
		BubbleSort(a, n);
	}
	else if (algorithm_name == "quick-sort") {
		QuickSort(a, 0, n - 1);
	}
	else if (algorithm_name == "shaker-sort") {
		ShakerSort(a, n);
	}
	else {
		cout << "Algorithm name is invalid!" << endl;
	}

	auto end = chrono::high_resolution_clock::now();
	elapsed = end - start;
}

void HandleInputOrder(int a[], int n, string input_order) {
	if (input_order == "-rand") {
		GenerateRandomData(a, n);
	}
	else if (input_order == "-nsorted") {
		GenerateNearlySortedData(a, n);
	}
	else if (input_order == "-sorted") {
		GenerateSortedData(a, n);
	}
	else if (input_order == "-rev") {
		GenerateReverseData(a, n);
	}
	else {
		cout << "Input order is invalid!" << endl;
	}
}

string InputOrder(string input_order) {
	string answer = "";
	if (input_order == "-rand") {
		answer = "Randomize";
	}
	else if (input_order == "-nsorted") {
		answer = "Nearly Sorted";
	}
	else if (input_order == "-sorted") {
		answer = "Sorted";
	}
	else if (input_order == "-rev") {
		answer = "Reversed";
	}
	else {
		answer = "Input order is invalid";
	}
	return answer;
}

void HandleOutputParameter(string output_parameter) {
	if (output_parameter == "-time") {
		cout << "Running time: " << elapsed.count() * secsToMilisecs << endl;
	}
	else if (output_parameter == "-comp") {
		cout << "Comparisons: " << count_comparison << endl;
	}
	else if (output_parameter == "-both") {
		cout << "Running time: " << elapsed.count() * secsToMilisecs << endl;
		cout << "Comparisons: " << count_comparison << endl;
	}
	else {
		cout << "Output parameter is invalid!" << endl;
	}
}

void WriteFile(int a[], int n, string file_name) {
	ofstream fout(file_name);
	if (!fout.is_open()) {
		cout << "Cannot open output.txt to write the sorted array" << endl;
		return;
	}

	fout << n << endl;
	for (int i = 0; i < n; i++) {
		fout << a[i] << ' ';
	}

	fout.close();
}


const int N = 1000000;
int a[N];
int b[N];

int main(int argc, char* argv[])
{
	string mode = argv[1];
	int n = 0; // Size of array
	const string input_file = "input.txt";
	const string output_file = "output.txt";

	if (mode == "-a") // Algorithm mode
	{
		cout << "ALGORITHM MODE" << endl;
		if (argc == 5) // Command 1 or command 3
		{
			if (!isanumber(argv[3])) // Command 1: [Execution file] -a [Algorithm] [Input filename] [Output parameter(s)]
			{
				string algorithm_name = argv[2];
				string input_file = argv[3];
				string output_parameter = argv[4];
				cout << "Algorithm: " << algorithm_name << endl;
				cout << "Input file: " << input_file << endl;

				ifstream fin(input_file);
				if (!fin.is_open()) {
					cout << "Cannot open input file" << endl;
					return 1;
				}

				fin >> n;
				cout << "Input size: " << n << endl;
				cout << "-------------------------\n";

				for (int i = 0; i < n; i++) {
					fin >> a[i];
				}

				fin.close();

				ExecuteSortAlgorithm(a, n, algorithm_name);
				WriteFile(a, n, output_file);
				HandleOutputParameter(output_parameter);
			}
			else // Command 3: [Execution file] -a [Algorithm] [Input size] [Output parameter(s)]
			{
				string algorithm_name = argv[2];
				string input_size = argv[3];
				string output_parameter = argv[4];
				cout << "Algorithm: " << algorithm_name << endl;
				cout << "Input size: " << input_size << endl;
				cout << endl;

				n = stoi(input_size);

				cout << "Input order: Randomize" << endl;
				cout << "-------------------------\n";
				GenerateRandomData(a, n);
				WriteFile(a, n, "input_1.txt");
				ExecuteSortAlgorithm(a, n, algorithm_name);
				HandleOutputParameter(output_parameter);
				cout << endl;

				cout << "Input order: Nearly Sorted" << endl;
				cout << "-------------------------\n";
				GenerateNearlySortedData(a, n);
				WriteFile(a, n, "input_2.txt");
				ExecuteSortAlgorithm(a, n, algorithm_name);
				HandleOutputParameter(output_parameter);
				cout << endl;

				cout << "Input order: Sorted" << endl;
				cout << "-------------------------\n";
				GenerateSortedData(a, n);
				WriteFile(a, n, "input_3.txt");
				ExecuteSortAlgorithm(a, n, algorithm_name);
				HandleOutputParameter(output_parameter);
				cout << endl;

				cout << "Input order: Reversed" << endl;
				cout << "-------------------------\n";
				GenerateReverseData(a, n);
				WriteFile(a, n, "input_4.txt");
				ExecuteSortAlgorithm(a, n, algorithm_name);
				HandleOutputParameter(output_parameter);
				cout << endl;
			}
		}
		else if (argc == 6) // Command 2: [Execution file] -a [Algorithm] [Input size] [Input order] [Output parameter(s)]
		{
			string algorithm_name = argv[2];
			string input_size = argv[3];
			string input_order = argv[4];
			string output_parameter = argv[5];
			cout << "Algorithm: " << algorithm_name << endl;
			cout << "Input size: " << input_size << endl;
			cout << "Input order: " << InputOrder(input_order) << endl;
			cout << "-------------------------\n";

			n = stoi(input_size);
			HandleInputOrder(a, n, input_order);
			WriteFile(a, n, input_file);
			ExecuteSortAlgorithm(a, n, algorithm_name);
			WriteFile(a, n, output_file);
			HandleOutputParameter(output_parameter);
		}
		else
		{
			cout << "Number of arguments is invalid!" << endl;
		}
	}
	else if (mode == "-c") // Comparison mode
	{
		cout << "COMPARE MODE" << endl;
		if (argc == 5) // Command 4: [Execution file] -c [Algorithm 1] [Algorithm 2] [Input filename]
		{
			string algorithm_name_1 = argv[2];
			string algorithm_name_2 = argv[3];
			string input_file = argv[4];
			cout << "Algorithm: " << algorithm_name_1 << " | " << algorithm_name_2 << endl;
			cout << "Input file: " << input_file << endl;

			ifstream fin(input_file);
			if (!fin.is_open()) {
				cout << "Cannot open input file" << endl;
				return 1;
			}

			fin >> n;
			cout << "Input size: " << n << endl;
			cout << "-------------------------\n";

			for (int i = 0; i < n; i++) {
				fin >> a[i];
				b[i] = a[i];
			}

			fin.close();

			ExecuteSortAlgorithm(a, n, algorithm_name_1);
			double running_time_1 = elapsed.count() * secsToMilisecs;
			int comparison_1 = count_comparison;
			ExecuteSortAlgorithm(b, n, algorithm_name_2);
			double running_time_2 = elapsed.count() * secsToMilisecs;
			int comparison_2 = count_comparison;

			cout << "Running time: " << running_time_1 << " | " << running_time_2 << endl;
			cout << "Comparisons: " << comparison_1 << " | " << comparison_2 << endl;
		}
		else if (argc == 6) // Command 5: [Execution file] -c [Algorithm 1] [Algorithm 2] [Input size] [Input order]
		{
			string algorithm_name_1 = argv[2];
			string algorithm_name_2 = argv[3];
			string input_size = argv[4];
			string input_order = argv[5];
			cout << "Algorithm: " << algorithm_name_1 << " | " << algorithm_name_2 << endl;
			cout << "Input file: " << input_size << endl;
			cout << "Input order: " << InputOrder(input_order) << endl;
			cout << "-------------------------\n";

			n = stoi(input_size);
			HandleInputOrder(a, n, input_order);
			for (int i = 0; i < n; i++) {
				b[i] = a[i];
			}
			WriteFile(a, n, input_file);

			ExecuteSortAlgorithm(a, n, algorithm_name_1);
			double running_time_1 = elapsed.count() * secsToMilisecs;
			int comparison_1 = count_comparison;
			ExecuteSortAlgorithm(b, n, algorithm_name_2);
			double running_time_2 = elapsed.count() * secsToMilisecs;
			int comparison_2 = count_comparison;

			cout << "Running time: " << running_time_1 << " | " << running_time_2 << endl;
			cout << "Comparisons: " << comparison_1 << " | " << comparison_2 << endl;
		}
		else
		{
			cout << "Number of arguments is invalid!" << endl;
		}
	}
	else
	{
		cout << "Mode is invalid!" << endl;
		return 1;
	}
	return 0;
}