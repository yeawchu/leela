//
//  utils.hpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include <utility>
#include <type_traits>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <chrono>
#include <cstring>
#include "vec.h"

bool has_suffix(const std::string &str, const std::string &suffix);
void filename_rename_ifexist(const std::string &filename);

//
// class templates
//

// timer
template<typename T = std::chrono::high_resolution_clock>
class timer
{
    typename T::time_point t1, t2;

public:
    void start()
    {
        t1 = T::now();
    }
    
    void stop()
    {
        t2 = T::now();
    }

    std::chrono::hours elapse_hours()
    {
        return std::chrono::duration_cast<std::chrono::hours>(t2 - t1);
    }
    
    std::chrono::minutes elapse_minutes()
    {
        return std::chrono::duration_cast<std::chrono::minutes>(t2 - t1);
    }

    std::chrono::seconds elapse_seconds()
    {
        return std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    }

    std::chrono::milliseconds elapse_milliseconds()
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    }

    std::chrono::microseconds elapse_microseconds()
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    }

    std::chrono::nanoseconds elapse_nanoseconds()
    {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
    }

    friend std::ostream& operator<<(std::ostream &os, timer &t)
    {
        auto hours_to_mins = std::chrono::minutes(t.elapse_hours());
        auto mins_to_secs = std::chrono::seconds(t.elapse_minutes());
        auto secs_to_millis = std::chrono::milliseconds(t.elapse_seconds());
        auto millis_to_micros = std::chrono::microseconds(t.elapse_milliseconds());
        auto micros_to_nanos = std::chrono::nanoseconds(t.elapse_microseconds());
        
        return os << t.elapse_hours().count() << " hrs, " << t.elapse_minutes().count() - hours_to_mins.count() << " mins, " << t.elapse_seconds().count() - mins_to_secs.count() << " s, " << t.elapse_milliseconds().count() - secs_to_millis.count() << " ms, " << t.elapse_microseconds().count() - millis_to_micros.count() << " mus, " << t.elapse_nanoseconds().count() - micros_to_nanos.count() << " ns";
    }
};

//
// function templates
//

// variadic min
template<typename T>
T vmin(T &&t)
{
    return std::forward<T>(t);
}

template<typename T0, typename T1, typename... Ts>
typename std::common_type<T0, T1, Ts...>::type vmin(T0 &&val1, T1 &&val2, Ts &&... vs)
{
    if (val2 < val1)
        return vmin(val2, std::forward<Ts>(vs)...);
    else
        return vmin(val1, std::forward<Ts>(vs)...);
}

// variadic max
template<typename T>
T vmax(T &&t)
{
    return std::forward<T>(t);
}

template<typename T0, typename T1, typename... Ts>
typename std::common_type<T0, T1, Ts...>::type vmax(T0 &&val1, T1 &&val2, Ts &&... vs)
{
    if (val2 > val1)
        return vmin(val2, std::forward<Ts>(vs)...);
    else
        return vmin(val1, std::forward<Ts>(vs)...);
}

// convert T to a string
template <typename T>
std::string to_string(const T value)
{
    std::ostringstream out;
    out << value;
    return out.str();
}

// convert T to a string with a pre-defined precision, eg 1.0042400000
template <typename T>
std::string to_string_with_precision(const T value, const int n = 6)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(n) << value;
    return out.str();
}

// convert T to a string with preceeding zeros, eg 000023
template <typename T>
std::string to_string_with_padding(const T value, const char &c, const int n = 6)
{
    std::ostringstream out;
    out << std::setfill(c) << std::setw(n) << value;
    return out.str();
}

// near equal float/double comparison
template<typename T>
bool near_equal(T a, T b)
{
    return std::nextafter(a, std::numeric_limits<T>::lowest()) <= b
    && std::nextafter(a, std::numeric_limits<T>::max()) >= b;
}

template<typename T>
bool near_equal(T a, T b, long factor /* a factor of epsilon */)
{
    T min_a = a - (a - std::nextafter(a, std::numeric_limits<T>::lowest())) * factor;
    T max_a = a + (std::nextafter(a, std::numeric_limits<T>::max()) - a) * factor;
    
    return min_a <= b && max_a >= b;
}

// write variable to file
template<typename T>
void file_write_var(std::ofstream &ofs, const std::string &label, const T &var)
{
    ofs << label << " " << var << std::endl;
}

// write boolalpha variable to file
template<typename T>
void file_write_var(std::ofstream &ofs, const std::string &label, const bool &var)
{
    ofs << label << " " << std::boolalpha << var << std::endl;
}

// write vec3 to file
template<typename T>
void file_write_vec(std::ofstream &ofs, const std::string &label, vec3<T> &vec)
{
    ofs << label << " " << vec << std::endl;
}

// write 1D array to file
template<typename T>
void file_write_array1d(std::ofstream &ofs, const std::string &label, std::vector<T> &array1d)
{
    ofs << label << " " << array1d.size();                          // 1d array of vectors
    for (long i = 0, i_size = array1d.size(); i != i_size; ++i)
    {
        ofs << " " << array1d[i];
    }
    ofs << std::endl;
}

// write 1D array of vec3 to file
template<typename T>
void file_write_array1dvec(std::ofstream &ofs, const std::string &label, std::vector<vec3<T>> &array1dvec)
{
    ofs << label << " " << array1dvec.size();                          // 1d array of vectors
    for (long i = 0, i_size = array1dvec.size(); i != i_size; ++i)
    {
        ofs << " " << array1dvec[i];
    }
    ofs << std::endl;
}

// write 2D array to file
template<typename T>
void file_write_array2d(std::ofstream &ofs, const std::string &label, std::vector<std::vector<T>> &array2d)
{
    ofs << label << " " << array2d.size();                          // 2d array of ints
    for (long i = 0, i_size = array2d.size(); i != i_size; ++i)
    {
        ofs << " " << array2d[i].size();
        for (long j = 0, j_size = array2d[i].size(); j != j_size; ++j)
        {
            ofs << " " << array2d[i][j];
        }
    }
    ofs << std::endl;
}

// read variable from file
template<typename T>
void file_read_label_var(const std::vector<std::string> &v, T &var, long &i)
{
    std::stringstream is;
    is << v[++i];
    is >> var;
}

// read boolapha variable from file
template<typename T>
void file_read_label_var(const std::vector<std::string> &v, bool &var, long &i)
{
    std::stringstream is;
    is << v[++i];
    is >> std::boolalpha >> var;
}

// read vec3 from file
template<typename T>
void file_read_label_vec(const std::vector<std::string> &v, vec3<T> &vec, long &i)
{
    for (long j = 0, j_size = vec.data.size(); j != j_size; ++j)
    {
        std::stringstream is;
        is << v[++i];
        is >> vec.data[j];
    }
}

// read 1D array from file
template<typename T>
void file_read_label_array1d(const std::vector<std::string> &v, std::vector<T> &array1d, long &i)
{
    std::stringstream is;
    long a_size = 0;
    
    is << v[++i];
    is >> a_size;               // read array size
    is.clear();
    
    array1d.resize(a_size);  // resize array
    
    for (long a = 0; a != a_size; ++a)
    {
        is << v[++i];
        is >> array1d[a];
        is.clear();
    }
}

// read 1D array of vec3 from file
template<typename T>
void file_read_label_array1dvec(const std::vector<std::string> &v, std::vector<vec3<T>> &array1dvec, long &i)
{
    std::stringstream is;
    long a_size = 0;
    
    is << v[++i];
    is >> a_size;               // read array size
    is.clear();
    
    array1dvec.resize(a_size);  // resize array
    
    for (long a = 0; a != a_size; ++a)
    {
        for (long j = 0, j_size = array1dvec[a].data.size(); j != j_size; ++j)
        {
            is << v[++i];
            is >> array1dvec[a].data[j];
            is.clear();
        }
    }
}

// read 2D array from file
template<typename T>
void file_read_label_array2d(const std::vector<std::string> &v, std::vector<std::vector<T>> &array2d, long &i)
{
    std::stringstream is;
    long a1_size = 0;
    
    is << v[++i];
    is >> a1_size;           // read particle array size
    is.clear();
    
    array2d.resize(a1_size);  // resize array
    
    for (long a1 = 0; a1 != a1_size; ++a1)
    {
        long a2_size = 0;
        
        is << v[++i];
        is >> a2_size;           // read particle's NN list size
        is.clear();
        
        array2d[a1].resize(a2_size);  // resize array
        
        for (long a2 = 0; a2 != a2_size; ++a2)
        {
            is << v[++i];
            is >> array2d[a1][a2];
            is.clear();
        }
    }
}

// tokeniser
// example usage: std::vector<std::string> v = tokenise<std::string>("I am, alive;", " ,;");
// to tokenise " ", "," and ";"
template<typename T>
std::vector<T> tokenise(const T &str, const T &delimiter)
{
    std::vector<T> v;                                       // create temporary vector
    typename T::size_type start = 0;                        // set string position to 0
    auto pos = str.find_first_of(delimiter, start);         // find position of first delimiter
    while (pos != T::npos)                                  // check if delimiter is end-of-line
    {
        if (pos != start)                                   // ignore empty tokens
            v.emplace_back(str, start, pos - start);        // push tokens into v
        start = pos + 1;                                    // setup start of next token position
        pos = str.find_first_of(delimiter, start);          // find position of next delimiter
    }
    if (start < str.length())                               // ignore trailing delimiter
        v.emplace_back(str, start, str.length() - start);   // push last token into v
    return v;                                               // return v to caller
}

#endif /* utils_hpp */
