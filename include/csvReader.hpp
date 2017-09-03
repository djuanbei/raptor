/**
 * @file   csvReader.hpp
 * @author Liyun Dai <dlyun2009@gmail.com>
 * @date   Sun Sep  3 11:43:18 2017
 *
 * @brief  csv util reader
 *
 *
 */
#include <iostream>
#include <locale>
#include <vector>
namespace raptor {
using namespace std;

struct csv_reader : std::ctype<char> {
  csv_reader() : std::ctype<char>(get_table()) {}
  static std::ctype_base::mask const* get_table() {
    static std::vector<std::ctype_base::mask> rc(table_size,
                                                 std::ctype_base::mask());

    rc[','] = std::ctype_base::space;
    rc['\n'] = std::ctype_base::space;
    rc[' '] = std::ctype_base::space;
    return &rc[0];
  }
};

struct csv_istream {
  std::istream& is_;
  csv_istream(std::istream& is) : is_(is) {
    is_.imbue(std::locale(std::locale(), new csv_reader()));
  }

  istream& operator>>(bool& val) { return is_ >> val; }
  istream& operator>>(short& val) { return is_ >> val; }
  istream& operator>>(unsigned short& val) { return is_ >> val; }
  istream& operator>>(int& val) { return is_ >> val; }
  istream& operator>>(unsigned int& val) { return is_ >> val; }
  istream& operator>>(long& val) { return is_ >> val; }
  istream& operator>>(unsigned long& val) { return is_ >> val; }

  istream& operator>>(float& val) { return is_ >> val; }
  istream& operator>>(double& val) { return is_ >> val; }
  istream& operator>>(long double& val) { return is_ >> val; }
  istream& operator>>(void*& val) { return is_ >> val; }

  istream& operator>>(string& str) { return is_ >> str; }

  istream& operator>>(streambuf* sb) { return is_ >> sb; }
};
}
