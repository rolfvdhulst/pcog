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
         start_edges_pos{t_start_edges_pos}, header{std::move(t_header)},
         file{t_file} {};
   /// Reads the graph as a dense graph.
   /// \return The dense graph, or std::nullopt if reading somehow failed.
   std::optional<DenseGraph> ReadAsDense();

 private:
   std::size_t vertices = 0;
   std::size_t edges_expected = 0;
   long start_edges_pos = 0;
   std::string header;
   std::ifstream &file;
};
/// Reads the header of a dimacs file. Returns std::nullopt if this somehow
/// failed. \param fileStream the stream to read from \return The data needed to
/// read the graph
std::optional<DimacsFileHeader> readDimacsHeader(std::ifstream &fileStream);

/// Writes a graph to a dimacs file
/// \param graph The graph to write
/// \param file the stream to write to
/// \param name The name of the graph
/// \param header The comment header for the dimacs file
/// \return true if writing was done succesfully, false otherwise
bool writeToDimacsFile(const DenseGraph &graph, std::ofstream &file,
                       const std::string &name, const std::string &header);

/// Writes the graph out the dot graph format
/// \param graph The graph to write out
/// \param stream The stream to write the dot file to
/// \return True if succesfull
bool graphToDot(const DenseGraph& graph, std::ostream& stream);
#endif // PCOG_TESTS_GRAPHIO_HPP
