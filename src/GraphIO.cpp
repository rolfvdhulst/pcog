//
// Created by rolf on 26-11-22.
//

#include "pcog/GraphIO.hpp"
namespace pcog {
std::optional<DimacsFileHeader> readDimacsHeader(std::ifstream &file) {
   if (!file.good()) {
      return std::nullopt;
   }
   std::string line;
   std::string header;
   std::string line_prefix = "p edge";

   while (std::getline(file, line)) {
      if (!line.compare(0, line_prefix.size(), line_prefix)) {
         std::string arguments = line.substr(6);
         const auto begin = arguments.find_first_not_of(' ');
         if (begin == std::string::npos) {
            break;
         }
         const auto end = arguments.find_last_not_of(' ');
         const auto substring = arguments.substr(begin, end - begin + 1);

         const auto space = substring.find_first_of(' ');

         const auto first_substring = substring.substr(0, space);
         const auto second_substring = substring.substr(space);

         try {
            std::size_t vertices = std::stoul(first_substring);
            std::size_t edges_expected = std::stoul(second_substring);
            long start_edges_pos = file.tellg();
            return DimacsFileHeader(vertices, edges_expected, start_edges_pos,
                                    header, file);
         } catch (std::invalid_argument &error) {
            return std::nullopt;
         }

      } else {
         header += line;
      }
   }
   return std::nullopt;
}

bool writeToDimacsFile(const DenseGraph &graph, std::ofstream &ofstream,
                       const std::string &name, const std::string &header) {
   if (!ofstream.good()) {
      return false;
   }
   // split header into multiple lines
   auto result = std::vector<std::string>{};
   auto ss = std::stringstream{header};
   for (std::string line; std::getline(ss, line, '\n');) {
      result.push_back(line);
   }
   const std::string comment_preamble = "c ";
   std::string firstLine = comment_preamble + "FILE: " + name + "\n";
   ofstream << firstLine;
   ofstream << comment_preamble << "\n";
   ofstream << comment_preamble << "DESCRIPTION: \n";
   for (const auto &line : result) {
      ofstream << comment_preamble << line << "\n";
   }
   std::size_t num_nodes = graph.numNodes();
   std::size_t num_edges = graph.numEdges();
   ofstream << "p edge " << num_nodes << " " << num_edges << "\n";
   std::string buffer;
   for (Node i = 0; i < num_nodes; ++i) {
      const auto &neighbourhood = graph.neighbourhood(i);
      Node neighbour = neighbourhood.find_next(i);
      while (neighbour != INVALID_NODE) {
         buffer += ("e " + std::to_string(i + 1) + " " +
                    std::to_string(neighbour + 1) +
                    "\n"); // all indices in the format start with 1...
         neighbour = neighbourhood.find_next(neighbour);
      }
      ofstream << buffer;
      buffer.clear();
   }
   return true;
}
std::optional<DenseGraph> DimacsFileHeader::ReadAsDense() {
   DenseGraph graph(vertices);
   file.seekg(start_edges_pos);
   std::string line;
   std::string line_prefix = "e ";
   while (std::getline(file, line)) {
      if (!line.compare(0, line_prefix.size(), line_prefix)) {
         auto substring = line.substr(2);
         std::size_t index = 0;
         try {
            std::size_t v1 = std::stoul(substring, &index);
            auto smallerSubstring = substring.substr(index);
            std::size_t v2 = std::stoul(smallerSubstring, &index);
            graph.addEdge(v1, v2);
         } catch (std::invalid_argument &error) {
            return std::nullopt;
         }
      }
   }
   return graph;
}
std::size_t DimacsFileHeader::numVertices() const { return vertices; }
std::size_t DimacsFileHeader::numExpectedEdges() const { return edges_expected; }
bool graphToDot(const DenseGraph &graph, std::ostream &stream) {
   assert(graph.isConsistent());
   std::size_t n = graph.numNodes();
   stream << "graph g{\n";

   for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = i; j < n; ++j) {
         if (graph.isEdge(i, j)) {
            stream << i << " -- " << j << ";\n";
         }
      }
   }
   stream << "}\n";
   stream.flush();
   return stream.good();
}
} // namespace pcog