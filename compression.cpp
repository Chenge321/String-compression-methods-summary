/**
 * @file compression.cpp
 * @author Chenge Liu (liu1217@mcmaster.ca)
 * @brief A program to encode the input file using four different compression methods and rank the methods based on their compression rate, 
 * and decode a file generated by the LZ77 method.
 * @version 1.0
 * @date 2021-12-19
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <fstream>
#include "classes.h"
using namespace std;
/**
 * @brief The main function to run encode or decode process
 * 
 * @param argc The number of input arguments from command line
 * @param argv The inputarguments from command line
 * @return int
 */
int main(int argc, char *argv[])
{
	// Check the invalid input command line arguments
	if (argc != 4 and argc != 5)
	{
		cout << "Please enter correct arguments." << '\n';
		cout
			<< "To encode, 1. E 2. Show result (S) or not (N) 3. input file name 4. output filename" << '\n';
		cout << "To decode, 1. D 2. encoded file name 3. decoded file name" << '\n';
		return -1;
	}
	if (*argv[1] != 'E' and *argv[1] != 'D')
	{
		cout << "Please enter correct arguments." << '\n';
		cout
			<< "To encode, 1. E 2. Show result (S) or not (N) 3. input file name 4. output filename" << '\n';
		cout << "To decode, 1. D 2. encoded file name 3. decoded file name" << '\n';
		return -1;
	}
	if (*argv[1] == 'E' and *argv[2] != 'S' and *argv[2] != 'N')
	{
		cout << "Please enter correct arguments." << '\n';
		cout
			<< "To encode, 1. E 2. Show result (S) or not (N) 3. input file name 4. output filename" << '\n';
		cout << "To decode, 1. D 2. encoded file name 3. decoded file name" << '\n';
		return -1;
	}
	if ((*argv[1] == 'E' and argc != 5) || (*argv[1] == 'D' and argc != 4))
	{
		cout << "Please enter correct arguments." << '\n';
		cout
			<< "To encode, 1. E 2. Show result (S) or not (N) 3. input file name 4. output filename" << '\n';
		cout << "To decode, 1. D 2. encoded file name 3. decoded file name" << '\n';
		return -1;
	}

	// This part does the encoding process
	if (*argv[1] == 'E')
	{
		lz77 lz77_coding;			  // lz77 encoding
		Arithmetic Arithmetic_coding; // Arithmetic encoding
		string filename1 = argv[3];	  // Input file name
		string filename2 = argv[4];
		fstream file;			  // File object
		deque<char_A> char_deque; // deque to store character
		char_A char_A_temp;
		string str; //Input string
		ifstream input(filename1);
		bool show_result = false; // Use this variable to decide to show the result or not
		if (*argv[2] == 'S')	  // Show the result
		{
			show_result = true;
		}
		double ari_rate, lzw_rate, lr_rate, lz77_rate; // Compression rates for methods
		// Read string from the input file
		file.open(filename1, ios::in | ios::binary);
		while (file.read((char *)&char_A_temp, A))
		{
			char_deque.push_back(char_A_temp);
		}
		char_deque.push_back(file.gcount() == 0 ? 0 : (char_A_temp & ((1 << 8 * (A - 1)) - 1)) | char_A((uint64_t)file.gcount() << 8 * (A - 1)));
		file.close();
		for (uint64_t i = 0, j = uint64_t(char_deque.size() - 1); i < j; i++)
		{
			if (char_deque[i] == 0)
			{
				show_result = false;
			}

			str.push_back(char_deque[i]);
		}
		// If the input file contains too many letters, do not show the result of each methods (to avoid redundant message)
		if (str.length() == 0)
		{
			cout << "The input file is empty" << endl;
			return -1;
		}
		// Doing the encoding process
		lzw_rate = lzw_encode(str, show_result); // LZW encoding
		ari_rate = Arithmetic_coding.runArithmetic(str);
		//cout << str.length() << '\n';
		lr_rate = RLencode(str, show_result); // RL encoding
		lz77_rate = lz77_coding.file_encode(filename1, true, filename2);
		// Rank four methods based on the compression rate
		map<double, string> method_name = {{
											   ari_rate,
											   "Arithmetic",
										   },
										   {
											   lzw_rate,
											   "LZW",
										   },
										   {
											   lr_rate,
											   "RLencode",
										   },
										   {
											   lz77_rate,
											   "lz77",
										   }};
		// Print the result
		cout << "The rank of methods (in increasing compression rate): " << '\n';
		for (auto x : method_name)
		{
			cout << x.first << " " << x.second << endl;
		}
	}
	// This part does the decoding process
	if (*argv[1] == 'D')
	{
		lz77_decode decode;
		decode.file_decode(argv[2], true, argv[3]);
	}
	return 0;
}
