require_relative '../support'
require_relative '../parsers/string_parser'

module Bioinform
  class JasparParser < StringParser
    def header_pat
      /(?<name>)/
    end

    def row_pat
      /[ACGT]\s*\[\s*(?<row>(#{number_pat}\s+)*#{number_pat})\s*\]\n?/
    end

    def scan_splitter
      scanner.scan(/(\/\/\n)+/)
    end

    def parse_matrix
      matrix = []
      while row_string = scan_row
        matrix << split_row(row_string)
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