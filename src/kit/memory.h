/*
 * Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2016-04-19 11:37:42
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-21 Thursday 10:00:34
 */

#ifndef MEMORY_H
#define MEMORY_H

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#if defined(_WIN32)
#include <Windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/types.h>
#include <sys/param.h>
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#else
#error "Unable to define getMemorySize( ) for an unknown OS."
#endif

size_t getMemorySize();
#endif

#if defined(__unix__)
void usedMemory(const string& lab="");
#endif