#pragma once

// ---------------------------------------------------------------------------------------------------------------
// --------------------- This file contains misc. utility methods for working on std::string ---------------------
// ---------------------------------------------------------------------------------------------------------------


#include <string>

/** Compares two strings without regard to case.
    @return true if the strings are identical */
bool EqualsIgnoringCase(const std::string& str1, const std::string& str2);
bool EqualsIgnoringCase(const char* str1, const char* str2);

/** Compares at most 'nCharacters' of two strings without regard to case.
    @return true if the strings are identical */
bool EqualsIgnoringCase(const std::string& str1, const std::string& str2, unsigned int nCharacters);

/** Trims out the provided characters (usually spaces) from the beginning and the end of the string */
void Trim(std::string& str, const char* characters = " ");

/** Call this member function to remove instances of ch from the string. 
    Comparisons for the character are case-sensitive.
    @return the number of characters removed. */
void Remove(std::string& str, char character);

/** Extracts the right-most or left-most characters in a std::string */
std::string Right(const std::string& input, size_t nChars);
std::string Left(const std::string& input, size_t nChars);

/** Finds the last occurrence of the given character in the provided string, -1 if the character is not found. */
int ReverseFind(const std::string& str, char ch);

/** Converts all letters to upper-case in the provided string. */
void MakeUpper(std::string& str);

/** Converts all letters to lower-case in the provided string. */
void MakeLower(std::string& str);

/** @return true if the provided string contains the given substring */
bool Contains(const std::string& string, const std::string& substring);

/** This function takes a string and removes any 'special' (ASCII code < 32) characters in it */
std::string CleanString(const std::string &in);

std::string SimplifyString(const std::string& str);