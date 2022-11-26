//
// Created by rolf on 26-11-22.
//

#ifndef PCOG_TESTS_GRAPHIO_HPP
#define PCOG_TESTS_GRAPHIO_HPP

///\file This file contains functions for writing and reading graphs

#include "DenseGraph.hpp"
#include <fstream>
#include <optional>
#include <string>
#include <utility>

class DimacsFileHeader {
 public:
   DimacsFileHeader(std::size_t t_vertices, std::size_t t_edges_expected,
                    long t_start_edges_pos, std::string t_header,
                    std::ifstream &t_file)
       : vertices{t_vertices}, edges_expected{t_edges_expected},
         start_edges_pos{t_start_edges_pos},header{std::move(t_header)},file{t_file} {         };
   std::optional<DenseGraph> ReadAsDense();

 private:
   std::size_t vertices = 0;
   std::size_t edges_expected = 0;
   long start_edges_pos = 0;
   std::string header;
   std::ifstream &file;
};
std::optional<DimacsFileHeader> readDimacsHeader(std::ifstream &fileStream);

bool writeToDimacsFile(const DenseGraph &graph, std::ofstream& file,
                       const std::string& name, const std::string &header);
#endif // PCOG_TESTS_GRAPHIO_HPP
