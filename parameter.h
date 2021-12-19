/**
 * @file parameter.h
 * @author Chenge Liu (liu1217@mcmaster.ca)
 * @brief Contains all self-defined variables that are using in the classes.h
 * @version 1.0
 * @date 2021-12-19
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#pragma once
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
constexpr uint64_t LogN = 8;
constexpr uint64_t LogL = 3;
typedef unsigned char char_A;
constexpr uint64_t A = sizeof(char_A);
constexpr uint64_t N = ((uint64_t)1 << LogN) - 1;
constexpr uint64_t L = 1 << LogL;
constexpr uint64_t end_number = 1;
constexpr uint64_t charactereof = '#';
