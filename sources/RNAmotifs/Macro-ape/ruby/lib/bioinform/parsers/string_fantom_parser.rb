require_relative '../support'
require_relative '../parsers/string_parser'

module Bioinform
  class StringFantomParser < StringParser
    def header_pat
      /NA (?<name>[\w.+:-]+)\n[\w\d]+ A C G T.*\n/
    end

    def row_pat
      /[\w\d]+ (?<row>(#{number_pat} )*#{number_pat})\n?/
    end

    def scan_splitter
      scanner.scan(/(\/\/\n)+/)
    end

    def parse_matrix
      matrix = []
      while row_string = scan_row
        matrix << split_row(row_string)[0,4]
      end
      matrix.transpose
    end

    def parse!
      scan_any_spaces
      scan_splitter
      name = parse_name
      matrix = parse_matrix
      Parser.parse!(matrix).tap{|result| result.name = name}
    end

  end
end