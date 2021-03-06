



<h1 align="center">String Compression Methods Summary</h1>


<p align="center">
  </a>
  <a href="https://www.linkedin.com/in/chenge-liu/">
    <img src="https://img.shields.io/badge/-LinkedIn-black.svg?style=plastic-square&logo=linkedin&colorB=555"
      alt="linkedin" />
  </a>
  <a href="https://github.com/Chenge321/String-compression-methods-summary/issues">
    <img src="https://img.shields.io/github/issues-raw/e-hengirmen/huffman-coding"
      alt="open issues" />
  </a>
</p>


C++ program to rank string compression methods based on their compression rate and decode files generated by the lz77 compression method
## Table of Contents
- [Table of Contents](#table-of-contents)
- [How does it work?](#how-does-it-work)
  - [Encoder:](#encoder)
  - [Decoder:](#decoder)
- [How to use it?](#how-to-use-it)
- [Note](#note)
- [Versions](#versions)
- [References](#references)


## How does it work?
_This program has two functions: The first is to encode the input file using four different compression methods (RL encoding, Arithmetic Encoding, LZW encoding, LZ77) and rank the methods based on their compression rate.
The second function is to decode a file generated by the LZ77 method, which is one of the methods adopted in the first part, then give an output file after the decoding process. You can check out documentation inside [classes.h](https://github.com/Chenge321/String-compression-methods-summary/blob/main/classes.h) and [compression.cpp](https://github.com/Chenge321/String-compression-methods-summary/blob/main/compression.cpp) files to help you understand how does each method work._


### Encoder:
The first function of this program is encoding the input file by four different methods:


1. Run Length Encoding (RL encoding)
   
   This is a very simple algorithm for lossless compression. It replaces repeated bytes with a simple description of the repeated characters and the number of times they are repeated; for example, if the input string is “aaabba”, the output of RL encoding is “3a2b1a”. Consequently, the method will have a bad compression rate when dealing with a non-repeating string. Although simple and very inefficient for normal compression, it can sometimes be very useful (eg. JPEG).


2. Arithmetic Encoding
   
   This is an encoding method based on the probabilities of letter occurrences in the string. It uses those probabilities to calculate a range for every letter, then uses the range of the last letter in the string to represent the encoded string; for example, if the input string is “ab”, then the arithmetic coding process is 0.5*0.5 = 0.25. The main idea behind this method is the same as in Huffman coding, and the reason I choose this method is that it has a better compression rate than Huffman coding and it uses more mathematical techniques. However, it consumes too much time to encode and decode for a long string.


3. LZW Encoding (Lempel-Ziv-Welch)
   
   LZW creates a dictionary (code table) to add the unknown string sequence to the dictionary. When encountering such a string sequence again, LZW replaces the sequence with the dictionary index number. As the number of bytes occupied by an index number is often much smaller than by the replaced string, the purpose of file compression is achieved. For example, if the input string contains “abc”, the dictionary will first add “abc” to its contents, and use the corresponding symbol to represent “abc” the next time it meets “abc” at the string. For the implementation, I choose to use the map containers as the dictionary and use different numbers to represent substrings in the encoding process. This method has a better compression rate when dealing with a highly repetitive string.


4. LZ77
   
    LZ77 is a compression method evolved by Abraham Lempel and Jacob Ziv in 1977. The main idea of LZ77 is the same as the LZW; this method also uses a dictionary and symbols to represent the encoded string. The difference is that LZ77 uses a different method to create and construct the dictionary. LZ77 uses a forward buffer and a sliding window (represented as Tree_mark class in the code). LZ77 first loads a portion of the data into a forward buffer. Once the phrase in the data passes through the forward buffer, it moves to the sliding window and becomes part of the dictionary. The main idea of the LZ77 algorithm is to constantly find the longest phrase in the forward buffer that can match the phrase in the dictionary. As you decode each tag, encode the tag as a character and copy it to the sliding window. Whenever a phrase marker is encountered, the corresponding offset is looked up in the sliding window, along with the phrase of the specified length found there. Each time a symbol tag is encountered, a symbol saved in the tag is generated. The method has a better compression rate and processing time than the basic LZW method, so it is still widely used today.




After the encoding process, the program will show the following things:


1. The encode table of LZW
2. compression rate of LZW
3. result code and compression rate of Arithmetic Encoding and RL Encoding
4. compression rate of lz77
5. Rank of four methods based on the compression rate(Increasing order)






### Decoder:
The second function of this program is to decode the encoded file generated by the first part of the program. The reason I choose an LZ77 type file to decode is that the LZ77 method is the most practical and commonly used of those four methods.




## How to use it?
If using GCC, please compile with "g++ compression.cpp -std=c++20 -o compression"
1. To use the encoding function, the first argument is "E", the second argument is "S" if you want to show the result of methods, or "N" to only show the compression rate and the rank. (Note: for any not text files (eg. image), the program will only show the compression rate and the rank automatically because the result of other kinds of files is meaningless.)The third argument is the name of the input text file, and the fourth argument is the name of the output file. In the following example, the program will encode the "test.txt" and show the result, then output an encoded file "test.lz77".
   
   For example:
   
   .\compression E S inefficient_input.txt test.lz77


   The processing table of LZW method
  
    a 97


    b 98


    c 99


    d 100
    
    e 101
    
    f 102
    
    g 103
    
    h 104
    
    i 105
    
    j 106
    
    k 107
    
    L 76


    M 77
    
    N 78
    
    1 49
    
    2 50
    
    3 51


    Compression rate of LZW: 1
    
    Resulting code of Arithmetic coding: 0.378906
    
    Compression rate of Arithmetic coding: 0.477941
    
    The result of RLencode: 
    1a1b1c1d1e1f1g1h1i1j1k1L1M1N111213
    
    Compression rate of RLencode: 2
    
    Compression rate of lz77 is 0.705882
    
    The rank of methods (in increasing compression rate):
    
    0.477941 Arithmetic
    
    0.705882 lz77
    
    1 LZW
    
    2 RLencode


    Example2:


    .\compression E N efficient_input.txt test.bin


    Compression rate of LZW: 0.6
    
    Resulting code of Arithmetic coding: 0.166677
    
    Compression rate of Arithmetic coding: 0.235714
    
    Compression rate of RLencode: 0.4
    
    Compression rate of lz77 is 0.171429
    
    The rank of methods (in increasing compression rate):
    
    0.171429 lz77
    
    0.235714 Arithmetic
    
    0.4 RLencode
    
    0.6 LZW
   
   
 
 
2. To use decoding:
    The first argument is "D", the second argument is the name of the encoded file, the third argument is the name of the output file that contains the decoded result. In the following example, the program will decode “test.bin” and write the result in “decode.txt”
 
    For example:
 
    ./compression D test.bin decode.txt
3. The program will check the invalid input, for example:


    ./compression




    Please enter correct arguments.
    
    To encode, 1. E 2. Show result (S) or not (N) 3. file name 4. output filename


    To decode, 1. D 2. encoded file name 3. decoded file name


## Note
1. This program can accept any format file as the input file, but the main purpose of this program is ranking string compression methods, for another format file, it may have a bad result, the output file size may be larger than the original file size.
2. Due to the limit of the double type of c++. The arithmetic coding will be not accurate when dealing with a big file. The results are for reference only.
3. The compression rates used to rank methods are calculated theoretically, the compression rate in practice will be worse. For example, the arithmetic encoding needs to store the occurrence probability of every letter. The LZW encoding needs to store the original dictionary. This compression rate only can be used for research purposes.
4. There are two sample inputs, in the “inefficient_input.txt”, which contains a string of non-repeating letters, methods will have a bad result of this sample input. In the “efficient_input.txt”, which contains a highly repetitive string. Methods will have a good result for this sample input. Also, there are two sample image inputs, they can be used as different types of input files.




## Versions
* [Version 1.0](https://github.com/Chenge321/String-compression-methods-summary) 
  * The final version of the program
  * Update the implementation of the LZW method
  * Fix the bug when compressing the file that contains whitespace characters
  * Update comments and README file
  * This program can accept any type of input files
  * The total number of codes is 881 lines
* [Version 0.1](https://github.com/Chenge321/Compression-methods-summary) 
  * The draft version of the program
  * Useful for educative purposes
  * The total number of codes is 918 lines
  




## References
1. https://en.wikipedia.org/wiki/Arithmetic_coding#:~:text=Arithmetic%20coding%20(AC)%20is%20a,as%20in%20the%20ASCII%20code.&text=It%20represents%20the%20current%20information,range%2C%20defined%20by%20two%20numbers.




2. https://www.geeksforgeeks.org/lzw-lempel-ziv-welch-compression-technique/
3. https://en.wikipedia.org/wiki/LZ77_and_LZ78
4. https://blog.csdn.net/luoshixian099/article/details/50331883
5. The slides from the CAS722 course (by instructor Neerja Mhaskar), the lecture notes from CSE701 (by Dr. Barak Shoshany).





